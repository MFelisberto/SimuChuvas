/*
 * hibrido.c — versão híbrida MPI + OpenMP, modelo de fases paralelas em dois níveis.
 *
 * Nível MPI (entre nós): igual ao mpi.c. O mapa é dividido em faixas horizontais, uma por
 * processo (na atlantica, 1 processo por nó), com uma linha fantasma em cima e outra embaixo.
 *
 * Nível OpenMP (dentro do nó): a faixa do processo é dividida de novo em sub-faixas, uma por
 * thread. Cada thread guarda a sua própria cópia das linhas vizinhas (fantasmas da thread),
 * então durante a Fase 1 nenhuma thread lê linhas que outra está escrevendo.
 *
 * Laço de fases (todos os processos e threads, até a convergência global):
 *   Fase 1 — processamento local: cada thread varre a sua sub-faixa (no máximo -k varreduras);
 *            o total de mudanças do processo sai de uma redução OpenMP.
 *   Fase 2 — verificação distribuída: MPI_Allreduce(OU lógico) entre os processos.
 *   Fase 3 — troca com os vizinhos: (a) a thread mestre troca as bordas da faixa com os
 *            processos vizinhos (MPI_Sendrecv); (b) cada thread atualiza as suas fantasmas a
 *            partir da memória compartilhada (ou das fantasmas MPI, nas pontas da faixa).
 *
 * Só a thread mestre chama MPI, sempre fora das regiões paralelas (MPI_THREAD_FUNNELED).
 */
#include "common.h"
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    int nivel;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &nivel);
    int rank, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (nivel < MPI_THREAD_FUNNELED) {
        if (rank == 0) fprintf(stderr, "erro: o MPI não suporta MPI_THREAD_FUNNELED\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    Params p;
    parse_args(argc, argv, &p);
    if (p.anim_path || p.out_path) {
        if (rank == 0) fprintf(stderr, "erro: -a e -o só no ./mpi\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    long R = p.rows, C = p.cols;
    int T = omp_get_max_threads();             /* threads por processo (OMP_NUM_THREADS) */

    MPI_Datatype row_t;
    MPI_Type_contiguous((int)C, MPI_FLOAT, &row_t);
    MPI_Type_commit(&row_t);

    /* ---------- distribuição entre processos (igual ao mpi.c) ---------- */
    float *all = NULL;
    long *inicio = malloc(sizeof(long) * (np + 1));
    if (rank == 0) {
        all = malloc(sizeof(float) * R * C);
        if (!all) MPI_Abort(MPI_COMM_WORLD, 2);
        dem_load(&p, all);
        divide_faixas(all, R, C, np, p.balance, inicio);
    }
    MPI_Bcast(inicio, np + 1, MPI_LONG, 0, MPI_COMM_WORLD);
    long r0 = inicio[rank], nr = inicio[rank + 1] - r0;
    if (nr < T) {
        if (rank == 0) fprintf(stderr, "erro: faixa com menos linhas (%ld) que threads (%d)\n", nr, T);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    float *Z = malloc(sizeof(float) * nr * C);
    float *W = malloc(sizeof(float) * nr * C);
    float *gtop = malloc(sizeof(float) * C), *gbot = malloc(sizeof(float) * C);
    unsigned char *dirty = malloc(nr);
    if (!Z || !W || !gtop || !gbot || !dirty) MPI_Abort(MPI_COMM_WORLD, 2);

    /* "primeiro toque" em paralelo: cada thread toca as linhas que deve processar, para as
     * páginas ficarem na memória do soquete dela (cada nó tem 2 soquetes) */
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < nr; i++) {
        memset(Z + i * C, 0, sizeof(float) * C);
        memset(W + i * C, 0, sizeof(float) * C);
    }

    int *counts = malloc(sizeof(int) * np), *displs = malloc(sizeof(int) * np);
    for (int r = 0; r < np; r++) {
        displs[r] = (int)inicio[r];
        counts[r] = (int)(inicio[r + 1] - inicio[r]);
    }
    MPI_Scatterv(all, counts, displs, row_t, Z, (int)nr, row_t, 0, MPI_COMM_WORLD);
    free(all);
    init_water(W, Z, r0, nr, &p);

    /* ---------- divisão da faixa entre as threads (mesmo critério, -b ou não) ---------- */
    long *sub = malloc(sizeof(long) * (T + 1));
    divide_faixas(Z, nr, C, T, p.balance, sub);

    int up   = (rank > 0)      ? rank - 1 : MPI_PROC_NULL;
    int down = (rank < np - 1) ? rank + 1 : MPI_PROC_NULL;
    for (long j = 0; j < C; j++) gtop[j] = gbot[j] = INFINITY;

    /* fantasmas de cada thread: cópia da linha logo acima e logo abaixo da sub-faixa */
    float **ttop = malloc(sizeof(float *) * T), **tbot = malloc(sizeof(float *) * T);
    double *t_thr = calloc(T, sizeof(double));            /* tempo de cálculo de cada thread */
    long rounds = 0, rows_visited = 0;
    double t_local = 0, t_check = 0, t_exch = 0, ta;

    memset(dirty, 1, nr);
    MPI_Barrier(MPI_COMM_WORLD);                          /* só para todos ligarem o relógio juntos */
    double t0 = MPI_Wtime();

    /* troca inicial entre processos */
    MPI_Sendrecv(W, 1, row_t, up, 0, gbot, 1, row_t, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(W + (nr - 1) * C, 1, row_t, down, 1, gtop, 1, row_t, up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    #pragma omp parallel
    {
        int t = omp_get_thread_num();
        ttop[t] = malloc(sizeof(float) * C);
        tbot[t] = malloc(sizeof(float) * C);
        long a = sub[t], b = sub[t + 1];
        memcpy(ttop[t], t == 0     ? gtop : W + (a - 1) * C, sizeof(float) * C);
        memcpy(tbot[t], t == T - 1 ? gbot : W + b * C,       sizeof(float) * C);
    }

    int *dir_thr = calloc(T, sizeof(int));                /* sentido da próxima varredura de cada thread */
    for (;;) {
        rounds++;

        /* ---------- Fase 1: processamento local (threads) ---------- */
        ta = MPI_Wtime();
        long changed = 0;
        #pragma omp parallel reduction(+ : changed, rows_visited)
        {
            int t = omp_get_thread_num();
            double tt = omp_get_wtime();
            long a = sub[t], b = sub[t + 1], n = b - a;
            long i_first = (r0 + a == 0) ? 1 : 0;         /* 1ª/última linha do mapa: borda fixa */
            long i_last  = (r0 + b == R) ? n - 2 : n - 1;
            changed += relax_block(W + a * C, Z + a * C, n, C, ttop[t], tbot[t], i_first, i_last,
                                   dirty + a, &dir_thr[t], p.eps, p.max_sweeps, &rows_visited);
            t_thr[t] += omp_get_wtime() - tt;
        }                                                 /* (barreira implícita do OpenMP) */
        t_local += MPI_Wtime() - ta;

        /* ---------- Fase 2: verificação distribuída (processos) ---------- */
        ta = MPI_Wtime();
        int any = (changed > 0);
        MPI_Allreduce(MPI_IN_PLACE, &any, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        t_check += MPI_Wtime() - ta;
        if (!any) break;

        /* ---------- Fase 3: troca com os vizinhos ---------- */
        ta = MPI_Wtime();
        /* (a) entre processos: bordas da faixa */
        MPI_Sendrecv(W, 1, row_t, up, 0, gbot, 1, row_t, down, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(W + (nr - 1) * C, 1, row_t, down, 1, gtop, 1, row_t, up, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /* (b) entre threads: cada uma copia as linhas vizinhas; se mudaram, revisita a borda */
        #pragma omp parallel
        {
            int t = omp_get_thread_num();
            long a = sub[t], b = sub[t + 1];
            long i_first = (r0 + a == 0) ? 1 : 0, i_last = (r0 + b == R) ? (b - a) - 2 : (b - a) - 1;
            const float *nt = t == 0     ? gtop : W + (a - 1) * C;
            const float *nb = t == T - 1 ? gbot : W + b * C;
            if (memcmp(nt, ttop[t], sizeof(float) * C) != 0) {
                memcpy(ttop[t], nt, sizeof(float) * C);
                if (i_first == 0 && i_last >= 0) dirty[a] = 1;
            }
            if (memcmp(nb, tbot[t], sizeof(float) * C) != 0) {
                memcpy(tbot[t], nb, sizeof(float) * C);
                if (i_last == (b - a) - 1 && i_last >= i_first) dirty[b - 1] = 1;
            }
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

    double tm[4] = { t_total, t_local, t_check, t_exch }, tmax[4];
    MPI_Reduce(tm, tmax, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    /* desbalanceamento entre TODOS os trabalhadores (np x T threads), como no mpi.c */
    double meu_max = 0, minha_soma = 0, glob_max, glob_soma;
    for (int t = 0; t < T; t++) { minha_soma += t_thr[t]; if (t_thr[t] > meu_max) meu_max = t_thr[t]; }
    MPI_Reduce(&meu_max, &glob_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&minha_soma, &glob_soma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    long sw_total;
    MPI_Reduce(&rows_visited, &sw_total, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double media = glob_soma / ((double)np * T);
        double desbal = media > 0 ? 100.0 * (glob_max - media) / media : 0.0;
        print_header();
        print_result("hibrido", np, T, &p, rounds, sw_total,
                     tmax[0], tmax[1], tmax[2], tmax[3], desbal, &g);
    }

    for (int t = 0; t < T; t++) { free(ttop[t]); free(tbot[t]); }
    free(ttop); free(tbot); free(t_thr); free(sub); free(dir_thr);
    MPI_Type_free(&row_t);
    free(Z); free(W); free(gtop); free(gbot); free(dirty);
    free(counts); free(displs); free(inicio);
    MPI_Finalize();
    return 0;
}
