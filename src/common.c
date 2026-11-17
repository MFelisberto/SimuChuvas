/* common.c — implementação das partes compartilhadas. Ver common.h. */
#define _POSIX_C_SOURCE 200809L
#include "common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static void usage(const char *prog)
{
    fprintf(stderr,
        "uso: %s [opções] [mapa.dem]\n"
        "  mapa.dem  relevo do RS gerado por tools/prepare_dem.py (padrão: dados/rs.dem)\n"
        "  -L h      nível d'água imposto nas bordas do mapa, em metros\n"
        "            (simula o Guaíba / a Lagoa dos Patos cheios; padrão: nenhum)\n"
        "  -e eps    desnível mínimo de escoamento (padrão 0)\n"
        "  -o f      grava as profundidades (float32) para o tools/render.py\n"
        "  -a f      grava um quadro por rodada para a animação (tools/animar.py)\n"
        "  -k n      no máximo n varreduras por rodada na Fase 1 (padrão: até estabilizar;\n"
        "            com -a, o padrão vira 1, para a animação mostrar o passo a passo)\n",
        prog);
    exit(1);
}

void parse_args(int argc, char **argv, Params *p)
{
    p->dem_path = "dados/rs.dem"; p->out_path = NULL; p->anim_path = NULL; p->max_sweeps = -1; p->level = -INFINITY; p->eps = 0.0f;
    /* Aceita as opções antes ou depois do nome do mapa. (O getopt do macOS para no
     * primeiro argumento que não é opção; aqui ele é guardado e a leitura continua.) */
    while (optind < argc) {
        int c = getopt(argc, argv, "L:e:o:a:k:h");
        if (c == -1) { p->dem_path = argv[optind++]; continue; }
        switch (c) {
        case 'L': p->level = (float)atof(optarg); break;
        case 'e': p->eps = (float)atof(optarg); break;
        case 'o': p->out_path = optarg; break;
        case 'a': p->anim_path = optarg; break;
        case 'k': p->max_sweeps = atoi(optarg); break;
        default:  usage(argv[0]);
        }
    }
    if (p->max_sweeps < 0) p->max_sweeps = p->anim_path ? 1 : 0;
    if (dem_read_header(p->dem_path, &p->rows, &p->cols) != 0) {
        fprintf(stderr, "erro: não consegui ler %s (gere com: python3 tools/prepare_dem.py)\n",
                p->dem_path);
        exit(1);
    }
    if (p->rows < 3 || p->cols < 3) usage(argv[0]);
}

static inline uint64_t mix64(uint64_t x)          /* finalizador splitmix64 (checksum) */
{
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27; x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;
    return x;
}

/* ------------------------------------------------------------------ */
/* Arquivos                                                            */
/* ------------------------------------------------------------------ */

int dem_read_header(const char *path, long *rows, long *cols)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    char magic[8]; int64_t r, c;
    int ok = fread(magic, 1, 8, f) == 8 && memcmp(magic, "FLOODDEM", 8) == 0 &&
             fread(&r, 8, 1, f) == 1 && fread(&c, 8, 1, f) == 1;
    fclose(f);
    if (!ok) return -1;
    *rows = (long)r; *cols = (long)c;
    return 0;
}

void dem_read_all(const char *path, float *Z, long rows, long cols)
{
    FILE *f = fopen(path, "rb");
    if (!f || fseek(f, 24, SEEK_SET) != 0 ||
        fread(Z, sizeof(float), (size_t)rows * cols, f) != (size_t)rows * cols) {
        fprintf(stderr, "erro lendo %s\n", path);
        exit(1);
    }
    fclose(f);
}

void write_depth(const char *path, const float *D, long rows, long cols)
{
    FILE *f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "erro abrindo %s\n", path); return; }
    int64_t r = rows, c = cols;
    fwrite("FLOODDEP", 1, 8, f); fwrite(&r, 8, 1, f); fwrite(&c, 8, 1, f);
    fwrite(D, sizeof(float), (size_t)rows * cols, f);
    fclose(f);
}

/* ------------------------------------------------------------------ */
/* Algoritmo                                                           */
/* ------------------------------------------------------------------ */

void init_water(float *W, const float *Z, long r0, long nr, const Params *p)
{
    for (long i = 0; i < nr; i++) {
        long gi = r0 + i;
        int edge_row = (gi == 0 || gi == p->rows - 1);
        for (long j = 0; j < p->cols; j++) {
            float z = Z[i * p->cols + j];
            int fixed = edge_row || j == 0 || j == p->cols - 1 || z <= NODATA_Z;
            /* borda (e sem dado): a água escoa livremente, nível = max(Z, nível imposto) */
            W[i * p->cols + j] = fixed ? fmaxf(z, p->level) : INFINITY;
        }
    }
}

static inline float min2(float a, float b) { return a < b ? a : b; }

long sweep(float *restrict W, const float *restrict Z, long nr, long cols,
           const float *wtop, const float *wbot,
           long i_first, long i_last, unsigned char *dirty,
           int dir, float eps, long *rows_visited)
{
    long changes = 0;
    int down  = (dir == 0 || dir == 2);   /* linhas de cima para baixo? */
    int right = (dir == 0 || dir == 3);   /* colunas da esquerda para a direita? */
    long nrows = i_last - i_first + 1;
    int prev_changed = 0;                 /* a linha anterior (no sentido da varredura) mudou? */

    for (long k = 0; k < nrows; k++) {
        long i = down ? i_first + k : i_last - k;
        /* só processa linhas marcadas ou logo depois de uma linha que acabou de mudar */
        if (!dirty[i] && !prev_changed) continue;
        dirty[i] = 0;
        (*rows_visited)++;

        const float *up = (i == 0)      ? wtop : W + (i - 1) * cols;
        const float *dn = (i == nr - 1) ? wbot : W + (i + 1) * cols;
        float *w = W + i * cols;
        const float *z = Z + i * cols;
        long rc = 0;

        for (long t = 1; t <= cols - 2; t++) {
            long j = right ? t : cols - 1 - t;
            float wc = w[j], zc = z[j];
            if (wc <= zc) continue;             /* já está no chão: não desce mais */
            float m = min2(min2(min2(up[j - 1], up[j]), min2(up[j + 1], w[j - 1])),
                           min2(min2(w[j + 1], dn[j - 1]), min2(dn[j], dn[j + 1])));
            float cand = m + eps;
            if (cand < zc) cand = zc;           /* max(Z, menor vizinho + eps) */
            if (cand < wc) { w[j] = cand; rc++; }
        }

        prev_changed = (rc > 0);
        if (rc > 0) {
            changes += rc;
            /* a própria linha e a linha "de trás" (já visitada nesta varredura)
             * podem ter ficado instáveis: revisitar na próxima varredura */
            dirty[i] = 1;
            long behind = down ? i - 1 : i + 1;
            if (behind >= i_first && behind <= i_last) dirty[behind] = 1;
        }
    }
    return changes;
}

long relax_block(float *W, const float *Z, long nr, long cols,
                 const float *wtop, const float *wbot,
                 long i_first, long i_last, unsigned char *dirty,
                 int *dir, float eps, int max_sweeps, long *rows_visited)
{
    long total = 0, c;
    int n = 0;
    if (i_first > i_last) return 0;
    do {
        c = sweep(W, Z, nr, cols, wtop, wbot, i_first, i_last, dirty, *dir, eps, rows_visited);
        *dir = (*dir + 1) & 3;
        total += c;
        n++;
    } while (c > 0 && (max_sweeps <= 0 || n < max_sweeps));
    /* parou por c == 0: nenhuma linha marcada => bloco estável.
     * parou pelo limite: ainda há linhas marcadas, que ficam para a próxima rodada. */
    return total;
}

void stats_rows(const float *W, const float *Z, long r0, long nr, long cols, Stats *s)
{
    s->checksum = 0; s->flooded = 0; s->volume = 0.0; s->max_depth = 0.0;
    for (long i = 0; i < nr; i++)
        for (long j = 0; j < cols; j++) {
            long k = i * cols + j;
            uint32_t bits; memcpy(&bits, &W[k], 4);
            uint64_t gidx = (uint64_t)(r0 + i) * (uint64_t)cols + (uint64_t)j;
            s->checksum += mix64(gidx * 0x9E3779B97F4A7C15ULL + bits);
            if (Z[k] <= NODATA_Z) continue;
            double d = (double)W[k] - (double)Z[k];
            if (d > 0) {
                s->flooded++;
                s->volume += d;
                if (d > s->max_depth) s->max_depth = d;
            }
        }
}

void print_header(void)
{
    printf("versao,np,linhas,colunas,rodadas,trabalho,"
           "t_total,t_local,t_verif,t_troca,checksum,inundadas,volume,prof_max\n");
}

void print_result(const char *version, int np, const Params *p,
                  long rounds, long rows_visited,
                  double t_total, double t_local, double t_check, double t_exch,
                  const Stats *s)
{
    printf("%s,%d,%ld,%ld,%ld,%.2f,%.4f,%.4f,%.4f,%.4f,%016llx,%ld,%.6e,%.3f\n",
           version, np, p->rows, p->cols, rounds, (double)rows_visited / p->rows,
           t_total, t_local, t_check, t_exch,
           (unsigned long long)s->checksum, s->flooded, s->volume, s->max_depth);
    fflush(stdout);
}
