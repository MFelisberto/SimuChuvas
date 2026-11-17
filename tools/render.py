#!/usr/bin/env python3
"""
render.py — desenha o resultado como um mapa: relevo colorido e sombreado, água presa
em tons de azul, contorno do estado (ou da cidade), cidades, legenda e escala.

    # Rio Grande do Sul inteiro
    mpirun -np 4 ./mpi dados/rs.dem -o dados/rs.dep
    python3 tools/render.py dados/rs.dem dados/rs.dep -o figuras/rs_mapa.png

    # só a cidade de Porto Alegre (gere antes: python3 tools/prepare_dem.py --regiao poa)
    mpirun -np 4 ./mpi dados/poa.dem -o dados/poa.dep
    python3 tools/render.py dados/poa.dem dados/poa.dep --recorte poa -o figuras/porto_alegre.png

    # Porto Alegre com o Guaíba cheio: use o mesmo valor no -L e no --nivel
    mpirun -np 4 ./mpi -L 5.35 dados/poa.dem -o dados/poa_cheia.dep
    python3 tools/render.py dados/poa.dem dados/poa_cheia.dep --recorte poa --nivel 5.35 \
        -o figuras/porto_alegre_cheia.png --subtitulo "Guaíba a 5,35 m (pico de 5/5/2024)"

Requer: numpy, matplotlib, pillow. Usa o .json do mapa (gerado pelo prepare_dem.py)
para as coordenadas e tools/limites_rs.geojson para os contornos.
"""
import argparse
import json
import os
import struct

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LightSource, LinearSegmentedColormap, LogNorm
from matplotlib.patches import Patch
from PIL import Image, ImageDraw

NODATA_Z = -1.0e30
M_PER_DEG = 111320.0
AQUI = os.path.dirname(os.path.abspath(__file__))

# relevo: verde nas baixadas -> bege -> marrom na serra -> quase branco no alto
TERRENO = LinearSegmentedColormap.from_list("terreno", [
    (0.00, "#6e9f5c"), (0.04, "#8fb36f"), (0.12, "#bccb8a"), (0.30, "#e3d6a0"),
    (0.55, "#c9a67a"), (0.80, "#9c7b5e"), (1.00, "#efe9e1")])
# água presa: azul-claro (raso) -> azul-marinho (fundo)
AGUA = LinearSegmentedColormap.from_list("agua", ["#9ecae1", "#4292c6", "#08519c", "#08306b"])
COR_MAR = np.array([0.73, 0.83, 0.91])       # lagoas, Guaíba e rios largos (água que já existe)
COR_FUNDO = np.array([0.95, 0.95, 0.95])     # sem dado (fora das bacias do RS: oceano, países vizinhos)

CIDADES_RS = {
    "Porto Alegre": (-30.03, -51.23), "Caxias do Sul": (-29.17, -51.18),
    "Pelotas": (-31.77, -52.34), "Rio Grande": (-32.03, -52.10),
    "Santa Maria": (-29.68, -53.81), "Passo Fundo": (-28.26, -52.41),
    "Uruguaiana": (-29.75, -57.09), "Sant'Ana do Livramento": (-30.89, -55.53),
    "Bagé": (-31.33, -54.11), "Lajeado": (-29.47, -51.96), "Torres": (-29.34, -49.73),
    "Chuí": (-33.69, -53.46),
}
LUGARES_POA = {
    "Centro Histórico": (-30.031, -51.230), "Aeroporto": (-29.994, -51.171),
    "Arena do Grêmio": (-29.974, -51.196), "Beira-Rio": (-30.066, -51.236),
    "Arquipélago": (-29.990, -51.262), "Sarandi": (-29.985, -51.125),
    "Restinga": (-30.150, -51.140), "Lago Guaíba": (-30.120, -51.275),
}
BAIRROS_CENTRO_NORTE = {
    "Centro Histórico": (-30.030, -51.229), "Cidade Baixa": (-30.042, -51.221),
    "Praia de Belas": (-30.049, -51.233), "Menino Deus": (-30.057, -51.224),
    "Floresta": (-30.017, -51.208), "São Geraldo": (-30.006, -51.204),
    "Navegantes": (-29.992, -51.206), "Humaitá": (-29.981, -51.193),
    "Anchieta": (-29.977, -51.168), "Sarandi": (-29.985, -51.125),
    "Aeroporto": (-29.994, -51.171), "Moinhos de Vento": (-30.025, -51.202),
    "Petrópolis": (-30.040, -51.185), "Cristal": (-30.095, -51.245),
    "Arquipélago": (-29.990, -51.262), "Lago Guaíba": (-30.075, -51.268),
}
CIDADES_ENCHENTE = {
    "Muçum": (-29.165, -51.869), "Encantado": (-29.236, -51.870), "Roca Sales": (-29.284, -51.868),
    "Lajeado": (-29.467, -51.961), "Estrela": (-29.502, -51.950), "Taquari": (-29.800, -51.865),
    "Montenegro": (-29.689, -51.461), "Novo Hamburgo": (-29.678, -51.131),
    "São Leopoldo": (-29.760, -51.147), "Canoas": (-29.918, -51.184),
    "Porto Alegre": (-30.033, -51.230), "Eldorado do Sul": (-30.085, -51.619),
    "Guaíba": (-30.114, -51.325),
}
# Cada recorte: título, contorno desenhado, lugares marcados e alguns ajustes visuais.
#   bbox: caixa fixa (lat_min, lat_max, lon_min, lon_max); sem ela, usa a caixa do contorno.
#   pontos: True = cidades com ponto e nome; False = só o nome do lugar.
RECORTES = {
    "rs":  {"titulo": "Onde a água fica presa no Rio Grande do Sul", "contorno": "Rio Grande do Sul",
            "lugares": CIDADES_RS, "max_px": 2600, "margem": 0.0, "pontos": True,
            "km": 100, "casas": 1, "zmax": 1300.0, "min_depth": 0.05, "legenda": "lower right"},
    "poa": {"titulo": "Onde a água fica presa em Porto Alegre", "contorno": "Porto Alegre",
            "lugares": LUGARES_POA, "max_px": 2000, "margem": 0.015, "pontos": False,
            "km": 5, "casas": 2, "zmax": None, "min_depth": 0.3, "legenda": "upper right"},
    "poa_centro": {"titulo": "Porto Alegre: Centro Histórico e zona norte", "contorno": "Porto Alegre",
                   "bbox": [-30.11, -29.93, -51.30, -51.08],
                   "lugares": BAIRROS_CENTRO_NORTE, "max_px": 2000, "margem": 0.0, "pontos": False,
                   "km": 2, "casas": 2, "zmax": None, "min_depth": 0.3, "legenda": "upper right"},
    "enchente": {"titulo": "Enchentes de 2024: Vale do Taquari e Grande Porto Alegre",
                 "contorno": "Rio Grande do Sul", "bbox": [-30.25, -29.05, -52.15, -50.95],
                 "lugares": CIDADES_ENCHENTE, "max_px": 2200, "margem": 0.0, "pontos": True,
                 "km": 20, "casas": 1, "zmax": 700.0, "min_depth": 0.3, "legenda": "upper right"},
}


def open_grid(path, magic):
    with open(path, "rb") as f:
        m = f.read(8)
        if m != magic:
            raise SystemExit(f"{path}: cabeçalho {m!r}, esperado {magic!r}")
        rows, cols = struct.unpack("<qq", f.read(16))
    return np.memmap(path, dtype="<f4", mode="r", offset=24, shape=(rows, cols))


def polygons(geom):
    """Lista de anéis externos (lon, lat) de um Polygon/MultiPolygon GeoJSON."""
    if geom["type"] == "Polygon":
        return [np.array(geom["coordinates"][0])]
    return [np.array(p[0]) for p in geom["coordinates"]]


def mask_from_polygons(rings, shape, lon_left, lat_top, px_deg):
    """Rasteriza o contorno na grade da imagem (True = dentro)."""
    img = Image.new("L", (shape[1], shape[0]), 0)
    draw = ImageDraw.Draw(img)
    for ring in rings:
        xy = [((lon - lon_left) / px_deg, (lat_top - lat) / px_deg) for lon, lat in ring]
        draw.polygon(xy, fill=1)
    return np.array(img, dtype=bool)


def mascara_feicao(limites, nome, shape, left, top, px):
    """Máscara (True = dentro) de uma feição do limites_rs.geojson na grade da imagem."""
    feat = next(f for f in limites["features"] if f["properties"]["nome"] == nome)
    return mask_from_polygons(polygons(feat["geometry"]), shape, left, top, px)


def scale_bar(ax, lon0, lat0, km, deg_per_km_lon):
    w = km * deg_per_km_lon
    ax.plot([lon0, lon0 + w], [lat0, lat0], color="black", lw=3, solid_capstyle="butt")
    ax.plot([lon0, lon0 + w / 2], [lat0, lat0], color="white", lw=1.6, solid_capstyle="butt")
    ax.text(lon0 + w / 2, lat0, f"{km:g} km\n", ha="center", va="bottom", fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dem", help="mapa (.dem)")
    ap.add_argument("dep", help="profundidades gravadas pelo programa com -o (.dep)")
    ap.add_argument("-o", "--out", default="figuras/rs_mapa.png")
    ap.add_argument("--recorte", choices=sorted(RECORTES), default="rs",
                    help="rs = estado inteiro (padrão); poa = município de Porto Alegre; "
                         "poa_centro = Centro Histórico e zona norte de Porto Alegre; "
                         "enchente = Vale do Taquari + Grande Porto Alegre")
    ap.add_argument("--titulo", help="título (padrão depende do recorte)")
    ap.add_argument("--subtitulo", default="", help='ex.: "Guaíba a 5,35 m"')
    ap.add_argument("--nivel", type=float,
                    help="o mesmo valor usado no -L do programa: desenha só o que foi alagado a partir "
                         "do Guaíba / da Lagoa / do mar (sem as poças locais)")
    ap.add_argument("--min-depth", type=float, help="profundidade mínima desenhada, em m "
                    "(padrão: 0,05 no rs; 0,3 no poa, que esconde as poças entre prédios e árvores)")
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()
    cfg = RECORTES[a.recorte]

    Zm, Dm = open_grid(a.dem, b"FLOODDEM"), open_grid(a.dep, b"FLOODDEP")
    if Zm.shape != Dm.shape:
        raise SystemExit("o .dem e o .dep têm tamanhos diferentes")
    meta_path = os.path.splitext(a.dem)[0] + ".json"
    if not os.path.exists(meta_path):
        raise SystemExit(f"falta {meta_path} (é gerado pelo prepare_dem.py)")
    res, _, lon_left, _, nres, lat_top = json.load(open(meta_path))["transform"]
    R, C = Zm.shape

    limites = json.load(open(os.path.join(AQUI, "limites_rs.geojson")))
    feat = next(f for f in limites["features"] if f["properties"]["nome"] == cfg["contorno"])
    rings = polygons(feat["geometry"])

    # recorte: caixa do contorno (+ margem), limitada ao mapa disponível
    if "bbox" in cfg:
        lat_min, lat_max, lon_min, lon_max = cfg["bbox"]
    else:
        allpts = np.vstack(rings)
        m = cfg["margem"]
        lon_min, lat_min = allpts.min(axis=0) - m
        lon_max, lat_max = allpts.max(axis=0) + m
    r0 = max(0, int((lat_top - lat_max) / -nres)); r1 = min(R, int(np.ceil((lat_top - lat_min) / -nres)))
    c0 = max(0, int((lon_min - lon_left) / res));  c1 = min(C, int(np.ceil((lon_max - lon_left) / res)))
    step = max(1, int(np.ceil(max(r1 - r0, c1 - c0) / cfg["max_px"])))

    z = np.array(Zm[r0:r1:step, c0:c1:step], dtype=np.float64)
    d = np.array(Dm[r0:r1:step, c0:c1:step], dtype=np.float64)
    px = res * step
    left, top = lon_left + c0 * res, lat_top + r0 * nres
    right, bottom = left + z.shape[1] * px, top - z.shape[0] * px
    lat_c = (top + bottom) / 2
    coslat = np.cos(np.radians(lat_c))
    dy_m, dx_m = px * M_PER_DEG, px * M_PER_DEG * coslat

    inside = mask_from_polygons(rings, z.shape, left, top, px)
    nodata = z < NODATA_Z / 2
    z[nodata] = 0.0
    gy, gx = np.gradient(z)
    sea = ~nodata & (np.abs(gx) + np.abs(gy) == 0) & (z < 6.0)   # superfície de água plana e baixa

    # relevo colorido por altitude + sombreado
    zmax = cfg["zmax"] or max(60.0, float(np.percentile(z[inside & ~nodata], 99.5)))
    colors = TERRENO(np.clip(z / zmax, 0, 1))[..., :3]
    ls = LightSource(azdeg=315, altdeg=40)
    rgb = ls.shade_rgb(colors, z, vert_exag=2.0 if cfg["pontos"] else 1.5,
                       dx=dx_m, dy=dy_m, blend_mode="soft")

    # fora do contorno: cinza claro e apagado
    gray = rgb.mean(axis=-1, keepdims=True)
    rgb = np.where(inside[..., None], rgb, 0.80 + 0.17 * gray)
    rgb[sea] = COR_MAR
    rgb[nodata] = COR_FUNDO
    # o ANADEM não tem dado no Guaíba nem na Lagoa dos Patos: dentro do contorno "com lagoas"
    # do RS, célula sem dado é lagoa
    lagoas = nodata & mascara_feicao(limites, "Rio Grande do Sul com lagoas", z.shape, left, top, px)
    rgb[lagoas] = COR_MAR

    # água presa (só dentro do contorno)
    min_depth = a.min_depth if a.min_depth is not None else cfg["min_depth"]
    wet = (d > min_depth) & ~sea & ~nodata & inside
    if a.nivel is not None:
        # alagado "de fora": o nível final da água é exatamente o nível imposto com -L
        wet &= np.abs(z + d - a.nivel) < 0.01
    norm = LogNorm(vmin=0.1, vmax=10)
    wcol = AGUA(norm(np.clip(d, 0.1, 10)))[..., :3]
    rgb[wet] = 0.08 * rgb[wet] + 0.92 * wcol[wet]

    # área com água dentro do contorno (cada pixel = dx_m x dy_m)
    area_km2 = wet.sum() * dx_m * dy_m / 1e6

    # ---------------- figura ----------------
    h, w = rgb.shape[:2]
    aspect = 1 / coslat
    fig_w = 9.0
    fig_h = fig_w * (h * aspect) / w + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(np.clip(rgb, 0, 1), extent=(left, right, bottom, top), aspect=aspect,
              interpolation="lanczos" if step > 1 else "bilinear")
    for ring in rings:
        ax.plot(ring[:, 0], ring[:, 1], color="#222222", lw=0.9 if a.recorte == "rs" else 1.2)

    halo = [pe.withStroke(linewidth=2.5, foreground="white")]
    for nome, (lat, lon) in cfg["lugares"].items():
        if not (left < lon < right and bottom < lat < top):
            continue
        if cfg["pontos"]:
            dx_txt = 0.012 * (right - left)
            ax.plot(lon, lat, "o", ms=3.5, color="#222222", mec="white", mew=0.8)
            ax.text(lon + dx_txt, lat + 0.6 * dx_txt, nome, fontsize=7.5, path_effects=halo)
        else:
            italic = nome.startswith("Lago")
            ax.text(lon, lat, nome, fontsize=8, ha="center", va="center",
                    style="italic" if italic else "normal",
                    color="#08306b" if italic else "#111111", path_effects=halo)

    ax.set_xlim(left, right); ax.set_ylim(bottom, top)
    ax.tick_params(labelsize=7)
    casas = cfg["casas"]
    ax.xaxis.set_major_formatter(lambda v, _: f"{abs(v):.{casas}f}°O")
    ax.yaxis.set_major_formatter(lambda v, _: f"{abs(v):.{casas}f}°S")

    km = cfg["km"]
    deg_per_km_lon = 1 / (M_PER_DEG / 1000 * coslat)
    scale_bar(ax, left + 0.05 * (right - left), bottom + 0.05 * (top - bottom), km, deg_per_km_lon)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=AGUA)
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.035, pad=0.06, aspect=40)
    cb.set_label("profundidade da água (m)", fontsize=8)
    cb.set_ticks([0.1, 0.3, 1, 3, 10]); cb.set_ticklabels(["0,1", "0,3", "1", "3", "10 ou mais"])
    cb.ax.tick_params(labelsize=7)
    rotulo = "alagado a partir do Guaíba/Lagoa" if a.nivel is not None else "água presa no relevo"
    ax.legend(handles=[Patch(color=COR_MAR, label="lagoas, Guaíba e rios largos"),
                       Patch(color=AGUA(0.6), label=rotulo)],
              loc=cfg["legenda"], fontsize=7.5, framealpha=0.9)

    titulo = a.titulo or cfg["titulo"]
    sub = (a.subtitulo + " · " if a.subtitulo else "") + \
          f"{'área alagada' if a.nivel is not None else 'área com água presa'} ≈ {area_km2:,.0f} km²".replace(",", ".")
    ax.set_title(titulo, fontsize=13, fontweight="bold", loc="left", pad=22)
    ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=9, color="#444444")
    cb.ax.text(1.0, -2.6, "Relevo: ANADEM (ANA / IPH-UFRGS) · Limites: IBGE · "
               "Modelo: preenchimento de depressões (Planchon & Darboux)",
               transform=cb.ax.transAxes, ha="right", va="top", fontsize=6, color="#666666")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=a.dpi, bbox_inches="tight")
    print(f"{a.out}: {'área alagada' if a.nivel is not None else 'área com água presa'} "
          f"≈ {area_km2:.0f} km² (1 px = {step} células)")


if __name__ == "__main__":
    main()
