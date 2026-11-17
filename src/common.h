/*
 * common.h — partes compartilhadas usadas pelo programa MPI (mpi.c).
 *
 * Problema: no mapa de relevo do Rio Grande do Sul (ANADEM), descobrir onde
 * a água da chuva consegue escoar até o oceano / para fora do estado e onde ela
 * fica presa, e com que profundidade. Algoritmo de Planchon & Darboux (2001).
 *
 * Para cada célula c do mapa, calcula-se o nível d'água W(c): a altura que a água
 * atinge em c antes de conseguir escapar. A profundidade empoçada é W(c) - Z(c).
 *
 * Regra de atualização (aplicada até nada mudar):
 *     W(c) <- min( W(c), max( Z(c), min_{n vizinho de c} W(n) + eps ) )
 *
 * Z(c)   = altitude do terreno na célula c (lida do mapa, nunca muda)
 * W(c)   = nível d'água em c (começa +inf no interior e max(Z, nível) na borda/oceano)
 * eps    = desnível mínimo imposto para a água escoar (0 = lagos planos)
 * vizinhos = as 8 células ao redor (vizinhança D8)
 */
#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stdio.h>

#define NODATA_Z (-1.0e30f)   /* células sem dado (ex.: oceano) funcionam como ralo */

typedef struct {
    long rows, cols;          /* dimensões do mapa (lidas do arquivo) */
    const char *dem_path;     /* mapa do RS no formato .dem (ver tools/prepare_dem.py) */
    const char *out_path;     /* saída opcional: profundidades (float32) */
    const char *anim_path;    /* saída opcional: quadros da animação, um por rodada */
    int max_sweeps;           /* máximo de varreduras por rodada na Fase 1 (0 = até estabilizar) */
    float level;              /* nível d'água imposto na borda do mapa (-inf = nenhum) */
    float eps;                /* desnível mínimo de escoamento */
} Params;

void   parse_args(int argc, char **argv, Params *p);

/* Arquivo .dem: "FLOODDEM" + int64 linhas + int64 colunas + float32 linha a linha */
int  dem_read_header(const char *path, long *rows, long *cols);
void dem_read_all(const char *path, float *Z, long rows, long cols);
void write_depth(const char *path, const float *D, long rows, long cols);

/* Inicializa W para as linhas [r0, r0+nr) da grade global. */
void init_water(float *W, const float *Z, long r0, long nr, const Params *p);

/*
 * Uma varredura Gauss-Seidel sobre um bloco de nr linhas. Só as linhas locais
 * i_first..i_last podem ser atualizadas (a primeira/última linha global e as
 * colunas 0 e cols-1 são borda fixa).
 *   wtop/wbot  linha fantasma acima/abaixo do bloco (valores de W do vizinho)
 *   dirty[i]   1 se a linha i pode estar instável (precisa ser visitada)
 *   dir (0..3) sentido da varredura; alternar propaga rápido em todas as direções
 * Retorna quantas células mudaram; *rows_visited conta as linhas processadas.
 */
long sweep(float *W, const float *Z, long nr, long cols,
           const float *wtop, const float *wbot,
           long i_first, long i_last, unsigned char *dirty,
           int dir, float eps, long *rows_visited);

/* Varre o bloco até ele ficar estável (nenhuma linha marcada), ou no máximo
 * max_sweeps varreduras se max_sweeps > 0.
 * *dir guarda o sentido da próxima varredura. Retorna o total de mudanças. */
long relax_block(float *W, const float *Z, long nr, long cols,
                 const float *wtop, const float *wbot,
                 long i_first, long i_last, unsigned char *dirty,
                 int *dir, float eps, int max_sweeps, long *rows_visited);

/* Estatísticas e checksum de linhas [r0, r0+nr). O checksum é uma soma módulo 2^64
 * (independe da ordem), então a soma dos checksums parciais é igual ao checksum
 * de 1 processo: serve para provar que qualquer número de processos dá o resultado idêntico. */
typedef struct {
    uint64_t checksum;
    long     flooded;     /* células com profundidade > 0 */
    double   volume;      /* soma das profundidades (m * célula) */
    double   max_depth;   /* maior profundidade (m) */
} Stats;

void stats_rows(const float *W, const float *Z, long r0, long nr, long cols, Stats *s);

/* Imprime a linha de resultado em CSV. "trabalho" = linhas processadas / linhas da
 * grade, ou seja, o trabalho total medido em varreduras completas equivalentes. */
void print_header(void);
void print_result(const char *version, int np, const Params *p,
                  long rounds, long rows_visited,
                  double t_total, double t_local, double t_check, double t_exch,
                  const Stats *s);

#endif
