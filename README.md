# Onde a água fica presa no Rio Grande do Sul — fases paralelas com MPI

O programa pega o **mapa de relevo do Rio Grande do Sul** e responde: se chover sobre o estado
inteiro, onde a água consegue escorrer até o oceano (ou para fora do estado) e onde ela fica
presa? E com que profundidade?

O cálculo é feito em paralelo com **MPI**, no modelo de **fases paralelas**. O objetivo do
trabalho é medir quanto ele acelera com 1, 4, 8, 16, 32 e 64 processos no cluster.

---

## Estrutura da pasta

```
alagamento/
├── README.md        este manual
├── Makefile         compila o programa (make)
├── run_grid.sh      job Slurm com as medições no cluster
├── src/             o cálculo, em C
│   ├── common.h     declarações e explicação da regra de atualização
│   ├── common.c     leitura do mapa, varredura, checksum, saída CSV
│   ├── mpi.c        versão MPI puro: distribuição do mapa e laço de fases com MPI
│   └── hibrido.c    versão híbrida MPI + OpenMP: faixas entre processos, sub-faixas entre threads
├── tools/           ferramentas em Python
│   ├── prepare_dem.py      baixa o relevo e gera os mapas em dados/
│   ├── render.py           desenha as figuras (RS inteiro ou só Porto Alegre)
│   ├── animar.py           monta o GIF das rodadas a partir dos quadros gravados com -a
│   └── limites_rs.geojson  contornos do RS e de Porto Alegre (IBGE)
├── dados/           mapas: rs.dem (RS, ~60 m), poa.dem (Porto Alegre, ~30 m) e seus .json
├── figuras/         figuras geradas pelo render.py
└── resultados/      CSVs das medições (criada pelo run_grid.sh)
```

Arquivos que aparecem com o uso e podem ser apagados quando quiser:
- **`mpi`** e **`hibrido`:** os executáveis, criados pelo `make`;
- **`dados/*.dep`:** profundidades gravadas com `-o`, refeitas em segundos;
- **`dados/*.anim`:** quadros da animação gravados com `-a`;
- **`dados/anadem_rs.tif`:** o arquivo original da ANA, baixado pelo `prepare_dem.py`. Sem ele, gerar um mapa novo baixa tudo de novo.

---

## Como usar

### 1. Instalar (no seu computador, dentro do ambiente conda)

```bash
pip install rasterio numpy pillow matplotlib
conda install -c conda-forge openmpi
export OMPI_CC=clang     # no Mac: faz o mpicc do conda usar o compilador do próprio Mac
```

- **rasterio, numpy:** ler o relevo e gerar os mapas;
- **pillow, matplotlib:** desenhar as figuras;
- **openmpi:** compilar e rodar o programa no seu computador. No cluster o MPI já existe.
- **`OMPI_CC=clang`:** sem isso, o `mpicc` do conda procura um compilador do conda que não vem
  instalado, e o `make` falha. Vale até fechar o terminal; num terminal novo, rode de novo.

### 2. Gerar os mapas (precisa de internet, uma vez só)

```bash
python tools/prepare_dem.py                  # RS inteiro, células de ~60 m  -> dados/rs.dem
python tools/prepare_dem.py --regiao poa     # Porto Alegre, células de ~30 m -> dados/poa.dem
python tools/prepare_dem.py --regiao enchente   # Taquari + Grande POA, ~30 m -> dados/enchente.dem
python tools/prepare_dem.py --res 30 -o dados/rs30.dem   # RS em ~30 m (opcional, ~2,8 GB)
```

O primeiro comando baixa o ANADEM das bacias do RS, publicado pela ANA, e guarda como
`dados/anadem_rs.tif` (o script mostra o progresso). Os outros comandos reaproveitam esse arquivo.
O arquivo da ANA vem em outra projeção (Albers, em metros); o script reprojeta para
latitude/longitude e recorta.

### 3. Compilar e rodar

```bash
make                                              # gera ./mpi e ./hibrido (make mpi: só o MPI puro)
mpirun -np 1 ./mpi                                # 1 processo (usa dados/rs.dem)
mpirun -np 8 ./mpi dados/rs.dem                   # 8 processos
mpirun -np 4 ./mpi dados/rs.dem -o dados/rs.dep   # também grava as profundidades
OMP_NUM_THREADS=4 mpirun -np 2 -x OMP_NUM_THREADS ./hibrido -b -k 2 dados/rs.dem   # híbrida: 2 processos x 4 threads
```

O `./hibrido` aceita as mesmas opções, menos `-a` e `-o` (figuras e animação são feitas com o `./mpi`).
No macOS, o `clang` da Apple não tem OpenMP: para compilar o `hibrido` no Mac, instale a `libomp`
(`conda install -c conda-forge llvm-openmp`) ou rode só no cluster.

| Opção | Significado |
|---|---|
| `mapa.dem` | o mapa a usar (padrão: `dados/rs.dem`) |
| `-L h` | impõe nível d'água de `h` metros no mar e nas bordas baixas (Guaíba e Lagoa "cheios", ver "Cenário da enchente") |
| `-e eps` | desnível mínimo para a água andar (padrão 0) |
| `-o arq.dep` | grava a profundidade da água em cada célula, para o `render.py` |
| `-a arq.anim` | grava um quadro por rodada, para a animação (`tools/animar.py`) |
| `-k n` | no máximo `n` varreduras por rodada na Fase 1 (`0` = até a faixa estabilizar, o padrão; com `-a`, o padrão vira 1). Na atlantica, `-k 2` foi o melhor |
| `-b` | divisão **balanceada**: faixas com o mesmo número de células com dado, em vez do mesmo número de linhas |
| `-v` | imprime a faixa, o tempo da Fase 1 e o trabalho de cada processo (para o gráfico de carga) |
| `-r f` | reduz o mapa: cada bloco de `f` x `f` células vira uma (a menor altitude do bloco), ou seja, o RS inteiro com `1/f²` das células. `f` pode ser fracionário. Usado na escalabilidade fraca |

### 4. Gerar as figuras

```bash
# RS inteiro
python tools/render.py dados/rs.dem dados/rs.dep -o figuras/rs_mapa.png

# Porto Alegre com o Guaíba a 5,35 m (use o mesmo número no -L e no --nivel)
mpirun -np 4 ./mpi -L 5.35 dados/poa.dem -o dados/poa_cheia.dep
python tools/render.py dados/poa.dem dados/poa_cheia.dep --recorte poa --nivel 5.35 \
    -o figuras/porto_alegre_cheia.png --subtitulo "Guaíba a 5,35 m (pico de 5/5/2024)"
```

- **`--recorte`:** `rs` (padrão, estado inteiro), `poa` (município de Porto Alegre), `poa_centro`
  (Centro Histórico e zona norte de Porto Alegre, com os bairros; use com o `dados/poa.dem`) ou
  `enchente` (Vale do Taquari + Grande Porto Alegre, a região mais atingida em 2024; use com o
  `dados/enchente.dem`).
- **`--nivel`:** pinta só o que foi alagado a partir do Guaíba, da Lagoa e do mar. Sem ele, o mapa
  mostra todas as poças do relevo; em Porto Alegre isso inclui muitas poças falsas entre prédios
  e árvores.

Para a região das enchentes (figuras e GIF):

```bash
mpirun -np 4 ./mpi -L 5.35 dados/enchente.dem -o dados/enchente_cheia.dep
python tools/render.py dados/enchente.dem dados/enchente_cheia.dep --recorte enchente --nivel 5.35 \
    -o figuras/enchente_cheia.png --subtitulo "Guaíba a 5,35 m (pico de 5/5/2024)"
mpirun -np 4 ./mpi dados/enchente.dem -a dados/enchente.anim
python tools/animar.py dados/enchente.anim dados/enchente.dem -o figuras/enchente_rodadas.gif \
    --titulo "Fases paralelas: Vale do Taquari e Grande Porto Alegre"
```

No Vale do Taquari a cheia de 2024 foi o rio subindo dezenas de metros, o que o modelo "banheira"
(`-L`, água entrando a partir do Guaíba) não representa: ali o mapa só mostra as baixadas onde a
água fica presa. O `-L` funciona bem na Grande Porto Alegre (delta do Jacuí, Caí, Sinos, Gravataí).

### 5. Animar as rodadas (opcional)

```bash
mpirun -np 4 ./mpi dados/rs.dem -a dados/rs.anim
python tools/animar.py dados/rs.anim dados/rs.dem -o figuras/rodadas.gif
```

O GIF mostra as fases paralelas acontecendo, rodada por rodada:
- **azul-escuro:** células que ainda não chegaram no nível final. O programa começa supondo
  tudo submerso (nível infinito), e essa "cortina" vai recuando conforme ele descobre por onde
  a água escapa: do oceano e das bordas para dentro do estado;
- **linhas vermelhas tracejadas:** as fronteiras entre as faixas dos processos (P0, P1...). Dá
  para ver a cortina parando em cada fronteira e só passar para a faixa vizinha na rodada seguinte;
- **relevo e azul:** o resultado final aparecendo.

Para o zoom no Centro Histórico e no Humaitá, com o Guaíba a 5,35 m e os nomes dos bairros:

```bash
python tools/prepare_dem.py --regiao centro
mpirun -np 4 ./mpi -L 5.35 dados/centro.dem -a dados/centro.anim
python tools/animar.py dados/centro.anim dados/centro.dem -o figuras/centro_rodadas.gif \
    --titulo "Fases paralelas: Centro Histórico e Humaitá (Guaíba a 5,35 m)" \
    --lugares poa_centro --nivel 5.35 --min-mudanca 0 --ms 350
```

A imagem estática do mesmo recorte:

```bash
mpirun -np 4 ./mpi -L 5.35 dados/centro.dem -o dados/centro_cheia.dep
python tools/render.py dados/centro.dem dados/centro_cheia.dep --recorte poa_centro --nivel 5.35 \
    -o figuras/centro_cheia.png --subtitulo "Guaíba a 5,35 m (pico de 5/5/2024)"
```

- **`--lugares poa_centro`:** escreve os nomes dos bairros no GIF;
- **`--nivel 5.35`:** como no `render.py`, pinta só o que foi alagado a partir do Guaíba (use o mesmo valor do `-L`);
- **`--min-mudanca 0 --ms 350`:** mostra todos os quadros, mais devagar (o recorte é pequeno e converge em poucas rodadas).

Com `-a`, cada processo faz **uma varredura por rodada** (`-k 1`), para a animação ter um quadro
por passo. Isso não muda o resultado (o checksum é o mesmo), mas muda o número de rodadas e o
tempo, então **não usem `-a` nas medições**. O arquivo `.anim` pode ter algumas centenas de MB;
`--passo 2` no `animar.py` usa 1 a cada 2 quadros e deixa o GIF menor.

### 6. Medir no cluster

**Pelo repositório (mais fácil):** o código e o mapa do RS estão no repositório. O `rs.dem` vai
comprimido e dividido em partes de 80 MB, porque o GitHub não aceita arquivos acima de 100 MB.

```bash
git clone https://github.com/MFelisberto/SimuChuvas.git
cd SimuChuvas
make dados        # remonta o dados/rs.dem (o run_grid.sh também faz isso sozinho)
sbatch run_grid.sh
```

Se o mapa mudar, rode `make partes` no seu computador para regerar as partes e faça commit.

**Sem o repositório:**

1. Copiar a pasta, sem o cache e sem os `.dep`:
   ```bash
   rsync -av --exclude anadem_rs.tif --exclude '*.dep' --exclude '*.anim' alagamento USUARIO@ENDERECO:~/
   ```
2. O `run_grid.sh` já está configurado para a atlantica (partição `LAD_geral`, Open MPI 4.1.1,
   4 nós de 8 núcleos físicos / 16 threads).
3. Submeter: `sbatch run_grid.sh`. Os resultados vão para `resultados/resultados_<job>.csv`, e a
   carga de cada processo para `resultados/carga_<job>.txt`.

O script compila, remonta o mapa, roda cada configuração 3 vezes (mude com `REPS=5 sbatch run_grid.sh`)
e grava a coluna `experimento`. Para rodar só alguns experimentos: `EXPS="fraca" sbatch run_grid.sh`.

| Experimento | O que mede |
|---|---|
| `escala` | 1 processo (base do speedup) e 4, 8, 16, 32 e 64 processos com divisão balanceada e `k = 2` |
| `divisao` | 32 e 64 processos com divisão **igual** (mesmo número de linhas), para comparar com a balanceada |
| `k` | 32 processos, divisão balanceada, com `k` = 0 (sem limite), 1, 4, 8 e 16 |
| `fraca` | escalabilidade fraca (só núcleos físicos): com `p` = 1, 4, 8, 16 e 32, o mapa com `p/32` das células (`-r` com `f = √(32/p)`); para cada carga, o tempo com 1 processo e com `p` processos |
| `carga` | uma execução de cada divisão com 32 processos e `-v` (vai para `carga_<job>.txt`) |
| `hibrido` | híbrida, escalabilidade forte: 4 processos (1 por nó, reservando os 8 núcleos do nó: `--map-by ppr:1:node:PE=8`) x 1, 2, 4, 8 e 16 threads (16 = HT, reservando as 16 threads de hardware) |
| `hfraca` | híbrida, escalabilidade fraca: 4 processos x 1, 2, 4 e 8 threads, com as mesmas cargas do `fraca` |
| `hsocket` | híbrida com 1 processo por soquete: 8 processos x 4 threads (32) e x 8 threads (64, HT) |
| `hlivre` | híbrida com 1 processo por nó sem reservar núcleos (`--bind-to none`), 4 x 4, 8 e 16 threads, para comparar com o modo padrão |
| `diag` | só imprime no `simuchuvas_<job>.out` a topologia do nó e o núcleo em que cada thread rodou |

Distribuição nos 4 nós: com até 32 processos, N/4 por nó, um por núcleo físico (`--map-by ppr:N:node
--bind-to core`); com 64, 16 por nó, dois por núcleo físico (HT). Ranks consecutivos ficam no mesmo
nó, para que só 3 fronteiras entre faixas passem pela rede.

---

## O problema

### O mapa

- **Fonte:** ANADEM, modelo digital de **terreno** de 30 m feito pelo IPH/UFRGS e pela ANA
  (Laipelt et al., 2024). Ele parte do Copernicus GLO-30 (relevo por satélite) e remove o viés
  da vegetação com aprendizado de máquina treinado com medições de laser de satélite (GEDI):
  o erro médio cai de 9,6 m para 1,5 m. A ANA publica um arquivo só com as bacias do RS.
- **Fora das bacias do RS** (oceano, Uruguai, Argentina) não há dado. Essas células funcionam
  como ralo: a água que chega nelas sai do mapa.
- **Recorte:** o retângulo que contém o RS, de Barra do Guarita (norte) ao Chuí (sul) e de Barra
  do Quaraí (oeste) a Torres (leste), uns 740 km × 770 km.
- **Grade:** cada célula é um quadradinho de terreno e guarda a altitude dele, de ~0 m no
  litoral e na Lagoa dos Patos até mais de 1.300 m nos Campos de Cima da Serra.
  - `rs.dem` (células de ~60 m): ~12,4 mil linhas × ~14,8 mil colunas ≈ **184 milhões de células**;
  - `rs30.dem` (células de ~30 m): ~24,8 mil × ~29,6 mil ≈ **735 milhões de células**.
- **Reamostragem:** para fazer células de 60 m a partir das de 30 m, o script fica com a mais
  baixa de cada grupo. Assim, os canais estreitos dos rios continuam abertos e não viram
  barragens falsas.

Taquari, Jacuí, Caí, Sinos e Gravataí desembocam no **Guaíba**. O Guaíba deságua na **Lagoa dos
Patos**, e a lagoa só chega ao oceano pela **barra de Rio Grande**. O norte e o oeste do estado
escoam pelo **Rio Uruguai**.

### O que é calculado

Para cada célula, o **nível d'água**: a altura que a água atinge ali antes de conseguir escapar.

**profundidade = nível d'água − altitude do terreno**

- **Profundidade 0:** existe um caminho sempre descendo até o oceano ou a borda do mapa, e a água escorre.
- **Profundidade maior que 0:** a célula está num buraco cercado de terreno mais alto; a água sobe
  até a beirada mais baixa e fica ali (lagos, banhados, várzeas).

### O algoritmo (Planchon & Darboux, 2001)

**Início:**
- na borda do mapa e no oceano, nível = altitude (a água escapa ali);
- no resto do estado, nível = infinito ("ainda não sei por onde a água sai daqui").

**Regra**, repetida em cada célula até nenhuma mudar:

**novo nível = o menor entre (nível atual) e (o maior entre a altitude e o menor nível dos 8 vizinhos + ε)**

- **menor nível dos 8 vizinhos:** o melhor caminho por onde a água pode sair;
- **o maior entre isso e a altitude:** a água nunca fica abaixo do chão;
- **o menor entre isso e o nível atual:** o nível só baixa, conforme aparecem caminhos de fuga melhores;
- **ε:** desnível mínimo para a água andar (opção `-e`, padrão 0, que gera lagos planos).

A informação "por aqui escoa" nasce no oceano e na borda e vai entrando no estado: passa pela
barra de Rio Grande, sobe a Lagoa dos Patos, chega ao Guaíba e se espalha pelos vales.

Detalhes da implementação (`sweep()` e `relax_block()` em `src/common.c`):
- **Varredura no lugar:** cada célula já usa o valor que o vizinho acabou de atualizar (estilo
  Gauss-Seidel), e o sentido alterna entre ↘ ↖ ↙ ↗ para a informação andar rápido em todas as direções.
- **Vetor `dirty`:** marca as linhas que ainda podem mudar; as varreduras seguintes só revisitam essas.
- **Resultado único:** como o nível só diminui e a resposta final é única, qualquer ordem de cálculo
  chega exatamente no mesmo mapa. Por isso o resultado com N processos pode ser comparado bit a bit
  com o de 1 processo (coluna `checksum`).

---

## As fases paralelas

O processo 0 lê o mapa e entrega a cada processo uma **faixa horizontal** do RS
(`MPI_Scatterv`): o processo 0 fica com o norte e o último com o extremo sul. Cada processo
guarda também duas **linhas fantasma**: a última linha da faixa de cima e a primeira da faixa de baixo.

Depois, todos repetem até terminar:

| Fase | O que acontece | Código |
|---|---|---|
| 1. Processamento local | cada processo aplica a regra na sua faixa até ela parar de mudar | `relax_block()` |
| 2. Verificação distribuída | `MPI_Allreduce` (OU lógico) de "minha faixa mudou?"; se ninguém mudou, terminou | `src/mpi.c` |
| 3. Troca com os vizinhos | `MPI_Sendrecv` das linhas de borda com a faixa de cima e a de baixo | `src/mpi.c` |

**Versão híbrida (`src/hibrido.c`).** O mesmo modelo em dois níveis: entre processos (um por nó)
continua igual; dentro do processo, a faixa é dividida de novo em **sub-faixas, uma por thread
OpenMP**, com o mesmo critério (`-b` ou não). Cada thread guarda a sua cópia das linhas vizinhas,
então na Fase 1 nenhuma thread lê o que outra está escrevendo. Fase 1: `omp parallel` com redução
do número de mudanças; Fase 2: `MPI_Allreduce` (só a thread mestre chama MPI, `MPI_THREAD_FUNNELED`);
Fase 3: a mestre troca as bordas com os processos vizinhos e depois cada thread atualiza as suas
cópias das linhas vizinhas na memória compartilhada. Antes do laço, cada thread copia as suas
linhas para memória nova, para elas ficarem no soquete da thread (cada nó tem 2 soquetes).

### O que esperar das medições

Mais processos dividem o cálculo em mais partes, mas também criam mais fronteiras. A informação
de escoamento atravessa **uma fronteira entre processos por rodada**, e cada rodada custa um
`MPI_Allreduce` e uma troca de mensagens. Por exemplo, a água do Guaíba só "descobre" que escapa
pela barra de Rio Grande depois que essa informação sobe, faixa por faixa, até Porto Alegre.
Estimativa pela posição das cidades no mapa:

| Processos | Altura de cada faixa | Fronteiras entre a barra de Rio Grande e Porto Alegre |
|---|---|---|
| 4 | ~185 km | ~1 |
| 8 | ~93 km | ~3 |
| 16 | ~46 km | ~5 |
| 32 | ~23 km | ~10 |
| 64 | ~12 km | ~20 |

Então o speedup não deve crescer sem limite: a partir de algum número de processos, o custo das
rodadas extras come o ganho de dividir o cálculo. O número real de rodadas sai na coluna `rodadas`.

### Como calcular speedup e eficiência

**speedup = t_total com 1 processo / t_total com N processos**

- quanto mais rápido ficou: speedup 20 = 20 vezes mais rápido que com 1 processo.

**eficiência = speedup / N**

- **N:** número de processos;
- quanto cada processo rendeu: 1,0 (100%) seria o ideal. Fica abaixo por causa de comunicação,
  espera e trabalho extra.

Usem a média das repetições de cada configuração.

---

## Saída (uma linha CSV por execução)

```
versao,np,threads,k,divisao,linhas,colunas,rodadas,trabalho,t_total,t_local,t_verif,t_troca,desbal,checksum,inundadas,volume,prof_max
```

(Os CSVs dos jobs 45565, 45572 e 45575 são de antes da coluna `threads`.)

| Coluna | Significado |
|---|---|
| `versao` | `mpi` ou `hibrido` |
| `np` | número de processos MPI |
| `threads` | threads OpenMP por processo (`1` no `mpi`) |
| `k` | máximo de varreduras por rodada (`0` = sem limite) |
| `divisao` | `igual` (mesmo número de linhas por faixa) ou `balanceada` (mesmo número de células com dado) |
| `linhas`, `colunas` | tamanho do mapa |
| `rodadas` | quantas vezes o laço de fases (1 → 2 → 3) rodou |
| `trabalho` | linhas processadas ÷ linhas do mapa = trabalho total em "varreduras completas do mapa". Comparado com o de 1 processo, mostra o trabalho extra causado pela divisão em faixas |
| `t_total` | tempo do laço de fases em segundos, sem a leitura do mapa (o maior entre os processos) |
| `t_local` | tempo na fase 1 |
| `t_verif` | tempo na fase 2; **inclui a espera pelas faixas mais lentas** (desbalanceamento: faixas com mais várzea e lagoa dão mais trabalho) |
| `t_troca` | tempo na fase 3 |
| `desbal` | desbalanceamento, a mesma métrica do T1/T2: quanto o maior `t_local` excede a média, em % (0% = perfeito). Na híbrida, calculado entre todas as threads de todos os processos |
| `checksum` | impressão digital do resultado; **tem que ser igual para qualquer número de processos** |
| `inundadas` | células com água presa |
| `volume` | soma das profundidades (m × célula) |
| `prof_max` | maior profundidade (m) |

---

## Cenário da enchente (`-L`)

Com `-L 5.35`, o mar e as bordas baixas do mapa ficam com água a 5,35 m. Essa água entra pela
barra de Rio Grande, sobe a Lagoa dos Patos e chega ao Guaíba; tudo abaixo de 5,35 m que está
ligado a eles fica embaixo d'água. É o modelo de inundação mais simples que existe (modelo
"banheira"). Para comparação, o Guaíba chegou a **5,35 m em 5 de maio de 2024**; a marca
anterior era 4,76 m, de 1941.

O `-L` muda o resultado (e o checksum), mas quase não muda o tempo. Para as medições, usem o mapa sem `-L`.

---

## Limitações (para o relatório)

- **Prédios e pontes continuam no mapa:** o ANADEM remove a vegetação, mas não corrige áreas
  urbanas (o próprio artigo diz que faltaram dados de treino para isso). Em Porto Alegre ainda
  podem aparecer poças falsas entre prédios, e uma ponte sobre um vale fundo pode virar uma
  "barragem". Com `--nivel` no `render.py`, essas poças não aparecem no mapa.
- **Resolução:** a 30–60 m, diques, comportas e o muro da Mauá praticamente não aparecem, e
  gargantas muito estreitas podem sumir do mapa. O ponto mais fundo do resultado (210 m) é o
  **cânion do Itaimbezinho** "cheio": no próprio ANADEM a saída do cânion aparece bloqueada,
  provavelmente porque a garganta do rio é mais estreita que a célula do mapa.
- **Cavas de mineração** aparecem como lagos fundos de verdade, como a mina de carvão de
  Candiota (~80 m).
- **Guaíba e Lagoa dos Patos sem dado:** o ANADEM não cobre essas lagoas. No cálculo, elas
  funcionam como ralo (a água que chega nelas escapa), o que equivale a supor que estão no
  nível normal; com `-L`, elas ficam no nível imposto. Nas figuras, aparecem em azul-claro.
- **`-L` aproximado:** o zero da régua do Cais Mauá não é o mesmo zero (nível do mar) do mapa,
  então `-L 5.35` é um cenário parecido com 2024, não a enchente exata.
- **Sem tempo nem vazão:** o modelo não tem chuva, correnteza nem duração; só a geometria de onde a água fica.

---

## Validação

- **Mesmo resultado com qualquer número de processos:** testado com 1, 2, 3, 4, 5, 8, 9 e 16
  processos, com e sem `-L` e `-e`; o `checksum` foi sempre idêntico.
- **Corretude contra um algoritmo independente:** as profundidades coincidem exatamente
  (erro 0,0) com as de um *priority-flood* (Barnes et al., 2014) implementado à parte em Python.
- **Valor de referência:** com o `dados/rs.dem` (ANADEM, células de ~60 m), o `checksum` é
  **`bb09d286e0336ace`** com 1 e com 4 processos. Todas as medições no cluster têm que dar
  exatamente esse valor (desde que usem o mesmo `rs.dem`).

Os tempos medidos no seu computador servem só para testar; para o relatório, valem as medições do cluster.
