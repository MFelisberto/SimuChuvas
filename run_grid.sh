#!/bin/bash
# run_grid.sh — bateria de medições no cluster (4 nós).
#
#   sbatch run_grid.sh                 # roda tudo e grava resultados/resultados_<jobid>.csv
#
# AJUSTEM antes de rodar: partição/fila e a forma de chamar o MPI usadas no LAD
# (o mesmo jeito que funcionou no T2). Os parâmetros abaixo assumem OpenMPI e nós
# com 8 núcleos físicos / 16 threads de hardware (HT) — confiram com `lscpu` num nó.
#
#SBATCH --job-name=alagamento
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=alagamento_%j.out
##SBATCH --partition=grid          # <- descomentem/ajustem o nome da partição

set -u
cd "$SLURM_SUBMIT_DIR"
make -s mpi

REPS=${REPS:-3}                             # repetições de cada configuração
mkdir -p resultados
OUT=resultados/resultados_${SLURM_JOB_ID:-local}.csv

# Entradas medidas: o mapa do RS (gerem antes com tools/prepare_dem.py).
#  - dados/rs.dem   : RS inteiro, células de ~60 m (~170 M células)
#  - dados/rs30.dem : RS inteiro, células de ~30 m (~690 M células), se existir:
#                     python3 tools/prepare_dem.py --res 30 -o dados/rs30.dem
make -s dados                               # remonta o dados/rs.dem a partir das partes do repositório
INPUTS=("dados/rs.dem")
[ -f dados/rs30.dem ] && INPUTS+=("dados/rs30.dem")
[ -f dados/rs.dem ] || { echo "falta dados/rs.dem"; exit 1; }

echo "entrada,rep,versao,np,linhas,colunas,rodadas,trabalho,t_total,t_local,t_verif,t_troca,checksum,inundadas,volume,prof_max" > "$OUT"

run() {   # run <entrada> <rep> <comando...>
  local inp="$1" rep="$2"; shift 2
  "$@" $inp | tail -1 | sed "s|^|\"$inp\",$rep,|" >> "$OUT"
}

for inp in "${INPUTS[@]}"; do
  for rep in $(seq 1 "$REPS"); do
    # 1 processo (base do speedup): o mapa inteiro num bloco só, sem paralelismo
    run "$inp" "$rep" mpirun -np 1 ./mpi

    # MPI com 4, 8, 16 e 32 processos, espalhados igualmente pelos 4 nós
    # (1, 2, 4 e 8 por nó; 8 por nó = todos os núcleos físicos)
    for np in 4 8 16 32; do
      run "$inp" "$rep" mpirun -np $np --map-by ppr:$((np / 4)):node --bind-to core ./mpi
    done
    # MPI com 64 processos = 16 por nó (usa o HT: 2 processos por núcleo físico)
    run "$inp" "$rep" mpirun -np 64 --map-by ppr:16:node --bind-to hwthread --use-hwthread-cpus ./mpi
  done
done

echo "pronto: $OUT"
