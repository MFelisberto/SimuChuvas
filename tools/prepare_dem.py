#!/usr/bin/env python3
"""
prepare_dem.py — baixa o ANADEM do Rio Grande do Sul e gera os mapas lidos pelo programa.

ANADEM: modelo digital de TERRENO de 30 m feito pelo IPH/UFRGS e pela ANA a partir do
Copernicus GLO-30, com o viés da vegetação removido (Laipelt et al., 2024, Remote Sensing).
A ANA publica um arquivo só com as bacias hidrográficas do RS:
    https://metadados.snirh.gov.br/geonetwork/srv/api/records/06eb254d-e45a-4c2b-9226-8e30d6afaa0d

O arquivo original vem em projeção Albers (metros). Este script o reprojeta para uma grade
em latitude/longitude, recorta e grava no formato .dem:
    "FLOODDEM" (8 bytes) + int64 linhas + int64 colunas + float32 linha a linha (linha 0 = norte)
Fora das bacias do RS não há dado: essas células viram -1e30 e funcionam como ralo.

Uso (precisa de internet só no primeiro comando, que baixa o arquivo da ANA uma vez):
    pip install rasterio numpy
    python tools/prepare_dem.py                  # RS inteiro, células de ~60 m -> dados/rs.dem
    python tools/prepare_dem.py --regiao poa     # Porto Alegre, células de ~30 m -> dados/poa.dem
    python tools/prepare_dem.py --regiao enchente   # Taquari + Grande POA, ~30 m -> dados/enchente.dem
    python tools/prepare_dem.py --regiao centro     # Centro Histórico e Humaitá, ~30 m -> dados/centro.dem
    python tools/prepare_dem.py --res 30 -o dados/rs30.dem    # RS inteiro em ~30 m
"""
import argparse
import json
import os
import struct
import sys
import time
import urllib.request

import numpy as np

URL_ANADEM_RS = ("https://metadados.snirh.gov.br/files/"
                 "06eb254d-e45a-4c2b-9226-8e30d6afaa0d/dem.tif")
ORIGEM = os.path.join("dados", "anadem_rs.tif")    # cópia local do arquivo da ANA
NODATA_Z = -1.0e30
M_POR_GRAU = 111320.0

# Recortes prontos (graus decimais: lat_min, lat_max, lon_min, lon_max) e tamanho da célula (m).
# rs:  retângulo que contém o estado, de Barra do Guarita ao Chuí e de Barra do Quaraí a Torres.
# poa: o município de Porto Alegre com uma margem.
REGIOES = {
    "rs":  {"bbox": [-33.76, -27.07, -57.66, -49.68], "res": 60},
    "poa": {"bbox": [-30.30, -29.90, -51.34, -50.98], "res": 30},
    # enchente: a região mais atingida em 2024, do Vale do Taquari (Muçum, Encantado, Lajeado)
    # até Porto Alegre, Canoas, Eldorado do Sul e o Vale do Sinos.
    "enchente": {"bbox": [-30.25, -29.05, -52.15, -50.95], "res": 30},
    # centro: zoom no Centro Histórico, Navegantes, Humaitá e Arquipélago (Porto Alegre)
    "centro": {"bbox": [-30.055, -29.955, -51.285, -51.145], "res": 30},
}


def baixar(url, destino):
    if os.path.exists(destino):
        return
    os.makedirs(os.path.dirname(os.path.abspath(destino)), exist_ok=True)
    print(f"baixando o ANADEM do RS (uma vez só):\n  {url}", file=sys.stderr)
    t0 = [time.time()]

    def progresso(blocos, tam_bloco, total):
        feito = blocos * tam_bloco
        if time.time() - t0[0] > 2 or feito >= total > 0:
            t0[0] = time.time()
            if total > 0:
                print(f"  {feito / 2**20:,.0f} de {total / 2**20:,.0f} MiB "
                      f"({100 * min(feito, total) / total:.0f}%)", file=sys.stderr)
            else:
                print(f"  {feito / 2**20:,.0f} MiB", file=sys.stderr)

    try:
        urllib.request.urlretrieve(url, destino + ".part", progresso)
    except Exception as e:
        sys.exit(f"não consegui baixar o arquivo da ANA: {e}")
    os.replace(destino + ".part", destino)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regiao", choices=sorted(REGIOES), default="rs",
                    help="rs = estado inteiro (padrão); poa = município de Porto Alegre; "
                         "enchente = Vale do Taquari + Grande Porto Alegre; "
                         "centro = Centro Histórico e Humaitá (Porto Alegre)")
    ap.add_argument("--bbox", nargs=4, type=float, metavar=("LAT_MIN", "LAT_MAX", "LON_MIN", "LON_MAX"),
                    help="recorte próprio em graus decimais (substitui o da --regiao)")
    ap.add_argument("--res", type=float, help="tamanho da célula em metros (padrão: 60 no rs, 30 no poa)")
    ap.add_argument("-o", "--out", help="arquivo de saída (padrão: dados/<regiao>.dem)")
    ap.add_argument("--origem", default=ORIGEM, help="cópia local do ANADEM (baixada se não existir)")
    a = ap.parse_args()

    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_origin
    from rasterio.vrt import WarpedVRT

    reg = REGIOES[a.regiao]
    bbox = a.bbox or reg["bbox"]
    res_m = a.res or reg["res"]
    out = a.out or os.path.join("dados", (a.regiao if a.bbox is None else "recorte") + ".dem")
    lat_min, lat_max, lon_min, lon_max = bbox
    if not (lat_min < lat_max and lon_min < lon_max):
        sys.exit("bbox inválida: use LAT_MIN LAT_MAX LON_MIN LON_MAX")

    baixar(URL_ANADEM_RS, a.origem)

    # grade de saída em lat/lon: células de res_m metros na direção norte-sul
    res_deg = res_m / M_POR_GRAU
    rows = int(round((lat_max - lat_min) / res_deg))
    cols = int(round((lon_max - lon_min) / res_deg))
    transform = from_origin(lon_min, lat_max, res_deg, res_deg)

    with rasterio.open(a.origem) as src:
        nodata = src.nodata
        # Reamostragem pelo MÍNIMO: ao juntar várias células de 30 m numa maior, fica a mais
        # baixa. Isso mantém abertos os canais estreitos dos rios (com média ou vizinho mais
        # próximo eles podem "sumir" e virar barragens falsas).
        resamp = Resampling.min if res_m > 31 else Resampling.bilinear
        with WarpedVRT(src, crs="EPSG:4326", transform=transform, width=cols, height=rows,
                       resampling=resamp, src_nodata=nodata, nodata=np.nan,
                       dtype="float32") as vrt:
            z = vrt.read(1).astype(np.float32)

    sem_dado = ~np.isfinite(z)
    if nodata is not None and np.isfinite(nodata):
        sem_dado |= (z == nodata)
    z[sem_dado] = NODATA_Z

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "wb") as f:
        f.write(b"FLOODDEM")
        f.write(struct.pack("<qq", rows, cols))
        f.write(np.ascontiguousarray(z, dtype="<f4").tobytes())

    meta = {"bbox": bbox, "rows": rows, "cols": cols, "res_m": res_m,
            "transform": list(transform)[:6], "fonte": "ANADEM (ANA / IPH-UFRGS)"}
    with open(os.path.splitext(out)[0] + ".json", "w") as f:
        json.dump(meta, f, indent=2)

    valid = z[~sem_dado]
    if valid.size == 0:
        sys.exit("o recorte não tem nenhum dado: confira a --bbox")
    print(f"{out}: {rows} x {cols} = {rows * cols / 1e6:.1f} M células "
          f"({100 * valid.size / z.size:.0f}% com dado), "
          f"altitude {valid.min():.1f} a {valid.max():.1f} m, {rows * cols * 4 / 2**20:.0f} MiB")


if __name__ == "__main__":
    main()
