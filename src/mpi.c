/*
 * mpi.c — versão paralela com MPI, modelo de fases paralelas.
 *
 * Decomposição: o mapa do RS é dividido em faixas horizontais de linhas
 * consecutivas, uma por processo (processo 0 = norte, último = extremo sul). Cada processo
 * guarda suas linhas de Z e W e duas linhas fantasma
 * (gtop = última linha do vizinho de cima, gbot = primeira linha do vizinho de baixo).
 *
 * Laço de fases (repetido por todos os processos até a convergência global):
 *   Fase 1 — processamento local: varre o próprio bloco até ele ficar estável
 *            (ou no máximo -k varreduras), só nas linhas que ainda podem mudar.
 *   Fase 2 — verificação distribuída: MPI_Allreduce(OU lógico) de "meu bloco mudou?".
 *            Se ninguém mudou, o resultado global está pronto.
 *   Fase 3 — troca com os vizinhos: cada processo envia sua primeira linha para
 *            rank-1 e sua última para rank+1 (MPI_Sendrecv), recebendo as fantasmas.
 */
#include "common.h"
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

/*
 * Animação (opção -a): a cada rodada, cada processo amostra 1 de cada F linhas e colunas da
 * sua faixa, o processo 0 junta tudo (MPI_Gatherv) e grava um quadro. Valor de cada pixel:
 *   profundidade W - Z;  +inf = nível ainda desconhecido ("sem saída" até agora);  NaN = sem dado.
 * Arquivo: "FLOODANI" + int64 (linhas, colunas, F, R, C, np) + int64 início da faixa de cada
 * processo; depois, para cada quadro, int64 rodada + float32 linhas*colunas.
 */
typedef struct {
    long F, rs, cs;            /* passo da amostragem e tamanho do quadro */
    long first, n;             /* primeira linha global amostrada desta faixa e quantas */
    int *cnt, *dsp;            /* contagens do Gatherv (só no processo 0) */
    float *buf, *all;
    FILE *f;
} Anim;

static long sampled_rows(long a, long b, long F)   /* linhas múltiplas de F em [a, b) */
{
    long first = (a + F - 1) / F * F;
    return first >= b ? 0 : (b - 1 - first) / F + 1;
}

static void anim_open(Anim *A, const Params *p, int rank, int np, const long *inicio)
{
    long r0 = inicio[rank], nr = inicio[rank + 1] - inicio[rank];
    long R = p->rows, C = p->cols;
    A->F = ((R > C ? R : C) + 899) / 900;          /* quadro com no máximo ~900 px */
    A->rs = (R + A->F - 1) / A->F;
    A->cs = (C + A->F - 1) / A->F;
    A->first = (r0 + A->F - 1) / A->F * A->F;
    A->n = sampled_rows(r0, r0 + nr, A->F);
    A->buf = malloc(sizeof(float) * (A->n * A->cs + 1));
    A->all = NULL; A->cnt = A->dsp = NULL; A->f = NULL;
    if (rank == 0) {
        A->all = malloc(sizeof(float) * A->rs * A->cs);
        A->cnt = malloc(sizeof(int) * np); A->dsp = malloc(sizeof(int) * np);
        long off = 0;
        A->f = fopen(p->anim_path, "wb");
        if (!A->f) { fprintf(stderr, "erro abrindo %s\n", p->anim_path); MPI_Abort(MPI_COMM_WORLD, 1); }
        int64_t hdr[6] = { A->rs, A->cs, A->F, R, C, np };
        fwrite("FLOODANI", 1, 8, A->f); fwrite(hdr, 8, 6, A->f);
        for (int r = 0; r < np; r++) {
            long a = inicio[r], b = inicio[r + 1];
            int64_t start = a;
            fwrite(&start, 8, 1, A->f);
            A->cnt[r] = (int)(sampled_rows(a, b, A->F) * A->cs);
            A->dsp[r] = (int)off;
            off += A->cnt[r];
        }
    }
}

static void anim_frame(Anim *A, long round, const float *W, const float *Z, long r0, long C)
{
    for (long k = 0; k < A->n; k++) {
        long i = A->first + k * A->F - r0;
        for (long jj = 0; jj < A->cs; jj++) {
            long j = jj * A->F;
            float z = Z[i * C + j], w = W[i * C + j];
            A->buf[k * A->cs + jj] = (z <= NODATA_Z) ? NAN : (isinf(w) ? INFINITY : w - z);
        }
    }
    MPI_Gatherv(A->buf, (int)(A->n * A->cs), MPI_FLOAT, A->all, A->cnt, A->dsp, MPI_FLOAT,
                0, MPI_COMM_WORLD);
    if (A->f) {
        int64_t r = round;
        fwrite(&r, 8, 1, A->f);
        fwrite(A->all, sizeof(float), A->rs * A->cs, A->f);
    }
}

static void anim_close(Anim *A)
{
    if (A->f) fclose(A->f);
    free(A->buf); free(A->all); free(A->cnt); free(A->dsp);
}

/*
 * Divide as R linhas do mapa em np faixas consecutivas.
 *   igual (padrão): o mesmo número de linhas em cada faixa;
 *   balanceada (-b): o mesmo número de CÉLULAS COM DADO em cada faixa. Linhas com muito
 *   oceano ou país vizinho ("sem dado", que quase não custam nada) ficam em faixas mais
 *   altas, e linhas só com terra em faixas mais baixas.
 * Cada linha pesa (células com dado + 1), para linhas vazias também serem distribuídas.
 */
static void divide_faixas(const float *Z, long R, long C, int np, int balance, long *inicio)
{
    inicio[0] = 0; inicio[np] = R;
    if (!balance) {
        for (int r = 1; r < np; r++) inicio[r] = (long)r * R / np;
        return;
    }
    double total = 0;
    double *peso = malloc(sizeof(double) * R);
    for (long i = 0; i < R; i++) {
        long n = 0;
        for (long j = 0; j < C; j++) n += (Z[i * C + j] > NODATA_Z);
        peso[i] = (double)n + 1.0;
        total += peso[i];
    }
    double acum = 0;
    long i = 0;
    for (int r = 1; r < np; r++) {
        double alvo = total * r / np;              /* a faixa r começa quando o acumulado passa daqui */
        while (i < R && acum + peso[i] <= alvo) acum += peso[i++];
        long min_i = inicio[r - 1] + 1;            /* toda faixa tem pelo menos 1 linha */
        long max_i = R - (np - r);                 /* e sobra pelo menos 1 linha para cada faixa seguinte */
        inicio[r] = i < min_i ? min_i : (i > max_i ? max_i : i);
        while (i < inicio[r]) acum += peso[i++];
    }
    free(peso);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    Params p;
    parse_args(argc, argv, &p);
    long R = p.rows, C = p.cols;
    if (R < np) {
        if (rank == 0) fprintf(stderr, "erro: menos linhas (%ld) que processos (%d)\n", R, np);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* um "tipo linha" para as contagens de Scatterv/Gatherv não estourarem int */
    MPI_Datatype row_t;
    MPI_Type_contiguous((int)C, MPI_FLOAT, &row_t);
    MPI_Type_commit(&row_t);

    /* O processo 0 lê o mapa do RS (fora da medição de tempo) e decide onde cada faixa
     * começa: inicio[r] = primeira linha do processo r; inicio[np] = R. */
    float *all = NULL;
    long *inicio = malloc(sizeof(long) * (np + 1));
    if (rank == 0) {
        all = malloc(sizeof(float) * R * C);
        if (!all) MPI_Abort(MPI_COMM_WORLD, 2);
        dem_load(&p, all);                        /* com -r f, o mapa já vem reduzido */
        divide_faixas(all, R, C, np, p.balance, inicio);
    }
    MPI_Bcast(inicio, np + 1, MPI_LONG, 0, MPI_COMM_WORLD);

    /* faixa de linhas deste processo: [r0, r0+nr) */
    long r0 = inicio[rank];
    long nr = inicio[rank + 1] - r0;

    float *Z = malloc(sizeof(float) * nr * C);
    float *W = malloc(sizeof(float) * nr * C);
    float *gtop = malloc(sizeof(float) * C), *gbot = malloc(sizeof(float) * C);
    float *ntop = malloc(sizeof(float) * C), *nbot = malloc(sizeof(float) * C);
    if (!Z || !W || !gtop || !gbot || !ntop || !nbot) MPI_Abort(MPI_COMM_WORLD, 2);

    int *counts = malloc(sizeof(int) * np), *displs = malloc(sizeof(int) * np);
    for (int r = 0; r < np; r++) {
        displs[r] = (int)inicio[r];
        counts[r] = (int)(inicio[r + 1] - inicio[r]);
    }
    /* entrega a cada processo a sua faixa (norte -> sul) */
    MPI_Scatterv(all, counts, displs, row_t, Z, (int)nr, row_t, 0, MPI_COMM_WORLD);
    free(all);
    init_water(W, Z, r0, nr, &p);

    int up   = (rank > 0)      ? rank - 1 : MPI_PROC_NULL;
    int down = (rank < np - 1) ? rank + 1 : MPI_PROC_NULL;
    for (long j = 0; j < C; j++) gtop[j] = gbot[j] = INFINITY;

    /* linhas locais que podem ser atualizadas (a primeira/última linha global é borda fixa) */
    long i_first = (r0 == 0) ? 1 : 0;
    long i_last  = (r0 + nr == R) ? nr - 2 : nr - 1;

    Anim A;
    if (p.anim_path) {
        anim_open(&A, &p, rank, np, inicio);
        anim_frame(&A, 0, W, Z, r0, C);           /* quadro 0: estado inicial */
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime(), t_local = 0, t_check = 0, t_exch = 0, ta;

    /* troca inicial: todos começam conhecendo a borda real dos vizinhos */
    MPI_Sendrecv(W, 1, row_t, up, 0, gbot, 1, row_t, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(W + (nr - 1) * C, 1, row_t, down, 1, gtop, 1, row_t, up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    long rounds = 0, rows_visited = 0;
    int dir = 0;
    /* dirty[i] = 1 se a linha local i pode estar instável e precisa ser visitada */
    unsigned char *dirty = malloc(nr);
    memset(dirty, 1, nr);
#define MARK(i) do { if ((i) >= i_first && (i) <= i_last) dirty[i] = 1; } while (0)
    for (;;) {
        rounds++;

        /* ---------- Fase 1: processamento local ---------- */
        ta = MPI_Wtime();
        /* (se nenhuma linha está marcada, só percorre o vetor dirty e retorna 0) */
        long changed = relax_block(W, Z, nr, C, gtop, gbot, i_first, i_last,
                                   dirty, &dir, p.eps, p.max_sweeps, &rows_visited);
        t_local += MPI_Wtime() - ta;
        if (p.anim_path) anim_frame(&A, rounds, W, Z, r0, C);   /* (o tempo disso entra no t_total) */

        /* ---------- Fase 2: verificação distribuída ---------- */
        ta = MPI_Wtime();
        int any = (changed > 0);
        MPI_Allreduce(MPI_IN_PLACE, &any, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        t_check += MPI_Wtime() - ta;
        if (!any) break;                     /* ninguém mudou: convergiu */

        /* ---------- Fase 3: troca com os vizinhos ---------- */
        ta = MPI_Wtime();
        MPI_Sendrecv(W, 1, row_t, up, 0, nbot, 1, row_t, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(W + (nr - 1) * C, 1, row_t, down, 1, ntop, 1, row_t, up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /* fantasma nova => marca a linha da borda correspondente para revisitar */
        if (up != MPI_PROC_NULL && memcmp(ntop, gtop, sizeof(float) * C) != 0) {
            float *t = gtop; gtop = ntop; ntop = t; MARK(0);
        }
        if (down != MPI_PROC_NULL && memcmp(nbot, gbot, sizeof(float) * C) != 0) {
            float *t = gbot; gbot = nbot; nbot = t; MARK(nr - 1);
        }
        t_exch += MPI_Wtime() - ta;
    }
    double t_total = MPI_Wtime() - t0;

    /* ---------- resultados ---------- */
    Stats s, g;
    stats_rows(W, Z, r0, nr, C, &s);
    MPI_Reduce(&s.checksum, &g.checksum, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&s.flooded, &g.flooded, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&s.volume, &g.volume, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&s.max_depth, &g.max_depth, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double tm[4] = { t_total, t_local, t_check, t_exch }, tmax[4], soma_local;
    MPI_Reduce(tm, tmax, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_local, &soma_local, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* desbalanceamento (mesma métrica do T1/T2): quanto a faixa mais lenta excede a média */
    double media_local = soma_local / np;
    double desbal = media_local > 0 ? 100.0 * (tmax[1] - media_local) / media_local : 0.0;

    if (p.verbose) {                    /* -v: uma linha por processo, para o gráfico de carga */
        double eu[4] = { (double)r0, (double)nr, t_local, (double)rows_visited }, *todos = NULL;
        if (rank == 0) todos = malloc(sizeof(double) * 4 * np);
        MPI_Gather(eu, 4, MPI_DOUBLE, todos, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            printf("# rank,linha_inicial,linhas,t_local,linhas_visitadas\n");
            for (int r = 0; r < np; r++)
                printf("# %d,%.0f,%.0f,%.4f,%.0f\n", r, todos[4*r], todos[4*r+1], todos[4*r+2], todos[4*r+3]);
            free(todos);
        }
    }
    long sw_total;
    MPI_Reduce(&rows_visited, &sw_total, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        print_header();
        print_result("mpi", np, &p, rounds, sw_total,
                     tmax[0], tmax[1], tmax[2], tmax[3], desbal, &g);
    }

    if (p.out_path) {
        for (long k = 0; k < nr * C; k++)
            W[k] = (Z[k] <= NODATA_Z) ? 0.0f : W[k] - Z[k];
        float *all = (rank == 0) ? malloc(sizeof(float) * R * C) : NULL;
        MPI_Gatherv(W, (int)nr, row_t, all, counts, displs, row_t, 0, MPI_COMM_WORLD);
        if (rank == 0) { write_depth(p.out_path, all, R, C); free(all); }
    }

    if (p.anim_path) anim_close(&A);
    MPI_Type_free(&row_t);
    free(Z); free(W); free(gtop); free(gbot); free(ntop); free(nbot);
    free(counts); free(displs); free(dirty); free(inicio);
    MPI_Finalize();
    return 0;
}
