#!/usr/bin/env python3
"""
animar.py — transforma os quadros gravados pelo programa (opção -a) num GIF que mostra
as fases paralelas rodada por rodada.

    mpirun -np 4 ./mpi dados/rs.dem -a dados/rs.anim     # 1 quadro por rodada (-k 1 automático)
    python tools/animar.py dados/rs.anim dados/rs.dem -o figuras/rodadas.gif

O que aparece em cada quadro:
  azul-escuro  células que ainda não chegaram no nível final: o programa começa supondo tudo
               submerso (nível infinito) e o nível vai baixando conforme ele descobre por onde
               a água escapa; essa "cortina" vai recuando até sobrar só o resultado
  azul         água que fica presa no relevo (mais escuro = mais fundo)
  relevo       onde a água já escoa
  linhas       fronteiras entre as faixas dos processos: a informação atravessa uma por rodada

Requer: numpy, matplotlib, pillow.
"""
import argparse
import json
import os
import struct

import numpy as np
from matplotlib.colors import LightSource, LogNorm
from PIL import Image, ImageDraw, ImageFont

from render import AGUA, AQUI, COR_FUNDO, COR_MAR, M_PER_DEG, NODATA_Z, RECORTES, TERRENO, \
    mascara_feicao, open_grid

COR_DESCONHECIDO = np.array([0.12, 0.20, 0.33])


def fonte(tam):
    """Fonte com acentos: a DejaVu Sans que vem junto com o matplotlib."""
    from matplotlib import font_manager
    try:
        return ImageFont.truetype(font_manager.findfont("DejaVu Sans"), tam)
    except OSError:
        return ImageFont.load_default()


def ler_anim(path):
    with open(path, "rb") as f:
        if f.read(8) != b"FLOODANI":
            raise SystemExit(f"{path} não é um arquivo de animação (gere com: ./mpi ... -a {path})")
        rs, cs, F, R, C, np_ = struct.unpack("<6q", f.read(48))
        inicios = struct.unpack(f"<{np_}q", f.read(8 * np_))
        cab = f.tell()
    tam_quadro = 8 + 4 * rs * cs
    n = (os.path.getsize(path) - cab) // tam_quadro
    quadros = np.memmap(path, dtype=np.dtype([("rodada", "<i8"), ("d", "<f4", (rs, cs))]),
                        mode="r", offset=cab, shape=(n,))
    return dict(rs=rs, cs=cs, F=F, R=R, C=C, np=np_, inicios=inicios, quadros=quadros)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("anim", help="quadros gravados com -a (.anim)")
    ap.add_argument("dem", help="o mesmo mapa usado na execução (.dem)")
    ap.add_argument("-o", "--out", default="figuras/rodadas.gif")
    ap.add_argument("--passo", type=int, default=1, help="usa 1 a cada N quadros (GIF menor)")
    ap.add_argument("--ms", type=int, default=150, help="duração de cada quadro (ms)")
    ap.add_argument("--min-mudanca", type=float, default=0.05,
                    help="pula quadros em que menos de X%% do mapa mudou desde o último quadro "
                         "mostrado (o fim da convergência quase não muda nada); 0 = mostra todos")
    ap.add_argument("--titulo", default="Fases paralelas no mapa do RS")
    ap.add_argument("--nivel", type=float, help="o mesmo valor do -L: pinta só o que foi alagado "
                    "a partir do Guaíba/Lagoa, sem as poças locais (como no render.py)")
    ap.add_argument("--lugares", choices=sorted(RECORTES), help="escreve no mapa os nomes dos "
                    "lugares desse recorte do render.py (ex.: poa_centro = bairros de Porto Alegre)")
    a = ap.parse_args()

    A = ler_anim(a.anim)
    F, rs, cs = A["F"], A["rs"], A["cs"]
    Zm = open_grid(a.dem, b"FLOODDEM")
    if Zm.shape != (A["R"], A["C"]):
        raise SystemExit("o mapa não é o mesmo usado para gravar a animação")
    z = np.array(Zm[::F, ::F][:rs, :cs], dtype=np.float64)

    # relevo de fundo (igual ao render.py), calculado uma vez só
    meta_path = os.path.splitext(a.dem)[0] + ".json"
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else None
    res = meta["transform"][0] if meta else 1 / 1800
    lat_c = (meta["bbox"][0] + meta["bbox"][1]) / 2 if meta else -30.0
    dy = res * F * M_PER_DEG
    dx = dy * np.cos(np.radians(lat_c))
    nodata = z < NODATA_Z / 2
    z[nodata] = 0.0
    cores = TERRENO(np.clip(z / 1300.0, 0, 1))[..., :3]
    relevo = LightSource(azdeg=315, altdeg=40).shade_rgb(cores, z, vert_exag=2.0, dx=dx, dy=dy,
                                                          blend_mode="soft")
    relevo[nodata] = COR_FUNDO
    gy, gx = np.gradient(z)
    mar = ~nodata & (np.abs(gx) + np.abs(gy) == 0) & (z < 6.0)      # rios largos e lagoas
    relevo[mar] = COR_MAR
    if meta:
        limites = json.load(open(os.path.join(AQUI, "limites_rs.geojson")))
        px = res * F
        lagoas = nodata & mascara_feicao(limites, "Rio Grande do Sul com lagoas", z.shape,
                                         meta["transform"][2], meta["transform"][5], px)
        relevo[lagoas] = COR_MAR

    # proporção certa: em latitude ~30°, um grau de longitude é ~13% menor que um de latitude
    larg = int(round(cs * np.cos(np.radians(lat_c))))
    escala = max(1, int(np.ceil(700 / max(larg, rs))))
    larg, alt = larg * escala, rs * escala
    topo, base = 64, 44
    f_tit, f_txt = fonte(20), fonte(14)

    norm = LogNorm(vmin=0.1, vmax=10)
    quadros = A["quadros"]
    idx = list(range(0, len(quadros), a.passo))
    if idx[-1] != len(quadros) - 1:
        idx.append(len(quadros) - 1)
    total_rodadas = int(quadros[-1]["rodada"])
    validos = (~nodata).sum()

    # o último quadro é o resultado final; em cada quadro, uma célula está "resolvida" quando
    # o nível dela já chegou no valor final (o nível só baixa, nunca sobe)
    final = np.nan_to_num(np.array(quadros[-1]["d"], dtype=np.float64), nan=0.0)
    presa = (final > 0.1) & ~nodata & ~mar
    if a.nivel is not None:      # só o que foi alagado a partir do Guaíba/Lagoa (como no render.py)
        presa &= np.abs(z + final - a.nivel) < 0.01
    agua = AGUA(norm(np.clip(final, 0.1, 10)))[..., :3]
    resultado = relevo.copy()
    resultado[presa] = 0.1 * relevo[presa] + 0.9 * agua[presa]

    # escolhe os quadros: pula os que quase não mudam nada em relação ao último mostrado
    escolhidos, ultimo = [], None
    for k in idx:
        d = np.nan_to_num(np.array(quadros[k]["d"], dtype=np.float64), nan=0.0)
        desc = (d > final + 1e-3) & ~nodata
        if ultimo is not None and k != idx[-1] and \
                100.0 * (desc != ultimo).sum() / max(1, validos) < a.min_mudanca:
            continue
        escolhidos.append(k); ultimo = desc
    idx = escolhidos

    imagens = []
    for k in idx:
        d = np.nan_to_num(np.array(quadros[k]["d"], dtype=np.float64), nan=0.0)
        desconhecido = (d > final + 1e-3) & ~nodata
        rgb = resultado.copy()
        rgb[desconhecido] = 0.2 * relevo[desconhecido] + 0.8 * COR_DESCONHECIDO

        mapa = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8)).resize((larg, alt), Image.NEAREST)
        img = Image.new("RGB", (larg, alt + topo + base), "white")
        img.paste(mapa, (0, topo))
        dr = ImageDraw.Draw(img)

        # fronteiras entre processos
        for p, ini in enumerate(A["inicios"]):
            y = topo + ini / F * escala
            if p > 0:
                for x0 in range(0, larg, 12):                       # linha tracejada
                    dr.line([(x0, y), (min(x0 + 7, larg), y)], fill=(230, 40, 40), width=2)
            if A["np"] <= 16:
                dr.text((6, y + 3), f"P{p}", font=f_txt, fill=(230, 40, 40),
                        stroke_width=2, stroke_fill="white")

        # nomes dos lugares (posição em graus -> pixel do quadro)
        if a.lugares and meta:
            lon0, lat0 = meta["transform"][2], meta["transform"][5]
            for nome, (lat, lon) in RECORTES[a.lugares]["lugares"].items():
                x = (lon - lon0) / (res * F) * (larg / cs)
                y = topo + (lat0 - lat) / (res * F) * escala
                if 0 <= x < larg and topo <= y < topo + alt:
                    w = dr.textlength(nome, font=f_txt)
                    dr.text((x - w / 2, y - 8), nome, font=f_txt, fill=(20, 20, 20),
                            stroke_width=2, stroke_fill="white")

        rodada = int(quadros[k]["rodada"])
        pct = 100.0 * desconhecido.sum() / max(1, validos)
        dr.text((10, 8), a.titulo, font=f_tit, fill="black")
        estado = "estado inicial" if rodada == 0 else f"rodada {rodada} de {total_rodadas}"
        dr.text((10, 36), f"{estado}  ·  {A['np']} processos  ·  ainda sem o nível final: "
                          f"{pct:.0f}% do mapa", font=f_txt, fill=(60, 60, 60))

        # legenda
        yb = topo + alt + 12
        x = 10
        for cor, nome in [(COR_DESCONHECIDO, "ainda sem o nível final"),
                          (AGUA(0.6)[:3], "alagado a partir do Guaíba" if a.nivel is not None else "água presa"), ((0.9, 0.3, 0.3), "fronteira entre processos")]:
            dr.rectangle([x, yb, x + 18, yb + 18], fill=tuple(int(255 * v) for v in cor))
            dr.text((x + 24, yb + 1), nome, font=f_txt, fill="black")
            x += 24 + int(dr.textlength(nome, font=f_txt)) + 22
        imagens.append(img)

    # uma paleta só para todos os quadros (feita a partir do primeiro e do último), para as
    # cores não "piscarem" de um quadro para o outro
    amostra = Image.new("RGB", (imagens[0].width, imagens[0].height * 2))
    amostra.paste(imagens[0], (0, 0)); amostra.paste(imagens[-1], (0, imagens[0].height))
    paleta = amostra.quantize(colors=192, method=Image.MEDIANCUT)
    imagens = [im.quantize(palette=paleta, dither=Image.Dither.NONE) for im in imagens]

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    # os primeiros quadros (onde quase tudo acontece) ficam mais tempo na tela; o último fica 2,5 s
    duracoes = [max(a.ms, 500) if i < 6 else a.ms for i in range(len(imagens) - 1)] + [2500]
    imagens[0].save(a.out, save_all=True, append_images=imagens[1:], duration=duracoes, loop=0)
    print(f"{a.out}: {len(imagens)} quadros, {larg}x{alt + topo + base} px, "
          f"{os.path.getsize(a.out) / 2**20:.1f} MiB")


if __name__ == "__main__":
    main()
