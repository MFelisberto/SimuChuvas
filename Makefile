MPICC   ?= mpicc
CFLAGS  ?= -O3 -std=c11 -Wall -Wextra
LDLIBS  = -lm

all: mpi hibrido

mpi: src/mpi.c src/common.c src/common.h
	$(MPICC) $(CFLAGS) -o $@ src/mpi.c src/common.c $(LDLIBS)

# versão híbrida MPI + OpenMP (no macOS com clang, precisa da libomp)
hibrido: src/hibrido.c src/common.c src/common.h
	$(MPICC) $(CFLAGS) -fopenmp -o $@ src/hibrido.c src/common.c $(LDLIBS)

# O rs.dem (~700 MB) vai no repositório comprimido e dividido em partes de 80 MB
# (o GitHub não aceita arquivos acima de 100 MB). "make dados" remonta o arquivo.
dados: dados/rs.dem

dados/rs.dem:
	cat dados/rs.dem.gz.part_* | gunzip > $@

# Gera as partes a partir do dados/rs.dem (rodar de novo se o mapa mudar).
partes:
	rm -f dados/rs.dem.gz.part_*
	gzip -6 -c dados/rs.dem | split -b 80m - dados/rs.dem.gz.part_

clean:
	rm -f mpi hibrido

.PHONY: all dados partes clean
