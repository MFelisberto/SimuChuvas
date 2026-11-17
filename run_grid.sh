#!/bin/bash
# run_grid.sh — bateria de medições na atlantica (4 nós, 8 núcleos físicos / 16 threads por nó).
#
#   sbatch run_grid.sh                 # grava resultados/resultados_<jobid>.csv
#   REPS=5 sbatch run_grid.sh          # 5 repetições em vez de 3
#   EXPS=fraca sbatch run_grid.sh      # só alguns experimentos (escala divisao k fraca carga)
#
#SBATCH --job-name=simuchuvas
#SBATCH --partition=LAD_geral
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=simuchuvas_%j.out

set -u
cd "$SLURM_SUBMIT_DIR"
make -s mpi
make -s dados                               # remonta o dados/rs.dem a partir das partes do repositório
[ -f dados/rs.dem ] || { echo "falta dados/rs.dem"; exit 1; }

REPS=${REPS:-3}                             # repetições de cada configuração
EXPS=${EXPS:-"escala divisao k fraca carga"}   # experimentos a rodar
tem() { [[ " $EXPS " == *" $1 "* ]]; }
MAPA=dados/rs.dem
mkdir -p resultados
OUT=resultados/resultados_${SLURM_JOB_ID:-local}.csv
CARGA=resultados/carga_${SLURM_JOB_ID:-local}.txt

# Como espalhar N processos pelos 4 nós: ranks consecutivos no mesmo nó (as faixas vizinhas
# trocam bordas, então assim só 3 fronteiras passam pela rede).
#   até 32: N/4 por nó, um por núcleo físico;   64: 16 por nó, dois por núcleo físico (HT).
mapa_de() {
  if [ "$1" -le 32 ]; then echo "--map-by ppr:$(( $1 / 4 )):node --bind-to core"
  else echo "--map-by ppr:16:node --bind-to hwthread --use-hwthread-cpus"; fi
}

echo "experimento,rep,versao,np,k,divisao,linhas,colunas,rodadas,trabalho,t_total,t_local,t_verif,t_troca,desbal,checksum,inundadas,volume,prof_max" > "$OUT"

run() {   # run <experimento> <rep> <np> <opções do programa...>
  local exp="$1" rep="$2" np="$3"; shift 3
  if [ "$np" -eq 1 ]; then mpirun -np 1 ./mpi "$@" $MAPA
  else mpirun -np "$np" $(mapa_de "$np") ./mpi "$@" $MAPA; fi | tail -1 | sed "s|^|$exp,$rep,|" >> "$OUT"
}

for rep in $(seq 1 "$REPS"); do
  # 1) Escalabilidade forte: base com 1 processo e a melhor configuração (divisão balanceada, k = 2)
  if tem escala; then
    run escala "$rep" 1
    for np in 4 8 16 32 64; do run escala "$rep" "$np" -b -k 2; done
  fi

  # 2) Divisão das faixas: mesmo número de linhas (igual) contra mesmo número de células com dado
  if tem divisao; then for np in 32 64; do run divisao "$rep" "$np" -k 2; done; fi

  # 3) Varreduras por rodada (k), 32 processos, divisão balanceada (0 = até a faixa estabilizar)
  if tem k; then for k in 0 1 4 8 16; do run k "$rep" 32 -b -k "$k"; done; fi

  # 4) Escalabilidade fraca (só núcleos físicos): carga proporcional a p. Com p processos, o RS
  #    inteiro com células f = sqrt(32/p) vezes maiores (-r f), ou seja, p/32 das células do mapa.
  #    Para cada carga: tempo sequencial (1 processo) e tempo paralelo (p processos).
  if tem fraca; then
    for np in 1 4 8 16 32; do
      f=$(awk -v p="$np" 'BEGIN { printf "%.6f", sqrt(32 / p) }')
      run fraca "$rep" 1 -r "$f"
      [ "$np" -gt 1 ] && run fraca "$rep" "$np" -b -k 2 -r "$f"
    done
  fi
done

# Carga de cada processo (uma execução de cada divisão, 32 processos), para o gráfico por rank
if tem carga; then
  for div in "" "-b"; do
    echo "## divisao=$([ -n "$div" ] && echo balanceada || echo igual)" >> "$CARGA"
    mpirun -np 32 $(mapa_de 32) ./mpi -k 2 $div -v $MAPA | grep "^#" >> "$CARGA"
  done
fi

echo "pronto: $OUT$(tem carga && echo " e $CARGA")"
