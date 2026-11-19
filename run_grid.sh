#!/bin/bash
# run_grid.sh — bateria de medições na atlantica (4 nós, 8 núcleos físicos / 16 threads por nó).
#
#   sbatch run_grid.sh                 # grava resultados/resultados_<jobid>.csv
#   REPS=5 sbatch run_grid.sh          # 5 repetições em vez de 3
#   EXPS=fraca sbatch run_grid.sh      # só alguns experimentos
#   experimentos: escala divisao k fraca carga (MPI puro) e diag hibrido hfraca hsocket hlivre (MPI + OpenMP)
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
make -s mpi hibrido
make -s dados                               # remonta o dados/rs.dem a partir das partes do repositório
[ -f dados/rs.dem ] || { echo "falta dados/rs.dem"; exit 1; }

REPS=${REPS:-3}                             # repetições de cada configuração
EXPS=${EXPS:-"escala divisao k fraca carga diag hibrido hfraca hsocket hlivre"}   # experimentos a rodar
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

echo "experimento,rep,versao,np,threads,k,divisao,linhas,colunas,rodadas,trabalho,t_total,t_local,t_verif,t_troca,desbal,checksum,inundadas,volume,prof_max" > "$OUT"

run() {   # run <experimento> <rep> <np> <opções do programa...>
  local exp="$1" rep="$2" np="$3"; shift 3
  if [ "$np" -eq 1 ]; then mpirun -np 1 ./mpi "$@" $MAPA
  else mpirun -np "$np" $(mapa_de "$np") ./mpi "$@" $MAPA; fi | tail -1 | sed "s|^|$exp,$rep,|" >> "$OUT"
}

# Versão híbrida: np processos MPI x T threads OpenMP (OMP_NUM_THREADS). Modos de distribuição:
#   no      (padrão) 1 processo por nó, que RESERVA o nó inteiro: 8 núcleos físicos (pe=8), ou as
#           16 threads de hardware com HT (pe=16); as threads se espalham por elas (spread)
#   soquete 1 processo por soquete (2 por nó), cada um com os 4 núcleos (ou 8 threads, HT) do
#           seu soquete
#   livre   1 processo por nó sem reserva (--bind-to none), como no job 45913, para comparar
runh() {  # runh <experimento> <rep> <modo> <np> <threads> <opções do programa...>
  local exp="$1" rep="$2" modo="$3" np="$4" T="$5"; shift 5
  local pe map places=cores bind=spread
  case "$modo" in
    no)      pe=8; [ "$T" -gt 8 ] && { pe=16; places=threads; }
             map="--map-by ppr:1:node:PE=$pe --bind-to core" ;;
    soquete) pe=4; bind=close; [ "$T" -gt 4 ] && { pe=8; places=threads; }
             map="--map-by ppr:1:socket:PE=$pe --bind-to core" ;;
    livre)   [ "$T" -gt 8 ] && places=threads
             map="--map-by ppr:1:node --bind-to none" ;;
  esac
  [ "$places" = threads ] && map="$map --use-hwthread-cpus"
  [ "$places" = threads ] && map="${map/--bind-to core/--bind-to hwthread}"
  local extra=""; [ -n "${OMP_DISPLAY_AFFINITY:-}" ] && extra="-x OMP_DISPLAY_AFFINITY -x OMP_AFFINITY_FORMAT"
  OMP_NUM_THREADS=$T OMP_PLACES=$places OMP_PROC_BIND=$bind \
    mpirun -np "$np" $map -x OMP_NUM_THREADS -x OMP_PLACES -x OMP_PROC_BIND $extra ./hibrido "$@" $MAPA \
    | tail -1 | sed "s|^|$exp,$rep,|" >> "$OUT"
}

# Diagnóstico (vai para o simuchuvas_<job>.out): topologia do nó e onde cada thread rodou
if tem diag; then
  echo "===== topologia do nó ====="
  lscpu | grep -iE "^(cpu\(s\)|thread|core|socket|numa)"
  command -v numactl >/dev/null && numactl -H | head -6
  for cfg in "no 4 8" "no 4 16" "soquete 8 4" "livre 4 8"; do
    set -- $cfg
    echo "===== afinidade: modo $1, $2 processos x $3 threads ====="
    OMP_DISPLAY_AFFINITY=true OMP_AFFINITY_FORMAT="host %H thread %n cpus %A" \
      runh diag 0 "$1" "$2" "$3" -b -k 2 -r 4 2>&1 | grep -E "^host|rror|not" | sort -V | head -70
    tail -1 "$OUT"
  done
fi

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

  # 5) Híbrida, escalabilidade forte: 4 processos (1 por nó) x 1, 2, 4, 8 e 16 threads
  #    (= 4 a 64 trabalhadores; 16 threads = HT), mapa inteiro, balanceada, k = 2
  if tem hibrido; then for T in 1 2 4 8 16; do runh hibrido "$rep" no 4 "$T" -b -k 2; done; fi

  # 6) Híbrida, escalabilidade fraca: p = 4 x T trabalhadores com a carga p/32 (mesmas cargas
  #    do experimento fraca, que dá o tempo sequencial de cada carga)
  if tem hfraca; then
    for T in 1 2 4 8; do
      f=$(awk -v p="$(( 4 * T ))" 'BEGIN { printf "%.6f", sqrt(32 / p) }')
      runh hfraca "$rep" no 4 "$T" -b -k 2 -r "$f"
    done
  fi

  # 7) Híbrida com 1 processo por soquete: 8 processos x 4 threads (32) e x 8 threads (64, HT)
  if tem hsocket; then for T in 4 8; do runh hsocket "$rep" soquete 8 "$T" -b -k 2; done; fi

  # 8) Híbrida sem reserva de núcleos (como no job 45913), para comparar com o modo "no"
  if tem hlivre; then for T in 4 8 16; do runh hlivre "$rep" livre 4 "$T" -b -k 2; done; fi
done

# Carga de cada processo (uma execução de cada divisão, 32 processos), para o gráfico por rank
if tem carga; then
  for div in "" "-b"; do
    echo "## divisao=$([ -n "$div" ] && echo balanceada || echo igual)" >> "$CARGA"
    mpirun -np 32 $(mapa_de 32) ./mpi -k 2 $div -v $MAPA | grep "^#" >> "$CARGA"
  done
fi

echo "pronto: $OUT$(tem carga && echo " e $CARGA")"
