import csv, collections
import numpy as np, matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
# Latin Modern (a mesma fonte do texto do relatório, como no T1/T2), tamanho de projeto 10
import glob, os, subprocess
def acha_lm():   # pasta com lmroman10-regular.otf (TeX Live / MacTeX)
    try:
        p = subprocess.run(["kpsewhich", "lmroman10-regular.otf"], capture_output=True, text=True).stdout.strip()
        if p: return os.path.dirname(p)
    except OSError: pass
    for d in glob.glob("/usr/share/tex*/fonts/opentype/public/lm") + glob.glob("/usr/local/texlive/*/texmf-dist/fonts/opentype/public/lm"):
        return d
LM = acha_lm()
FONTE = "serif"
if LM:
    fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if "Latin Modern" not in f.name]
    for nome in ("lmroman10-regular", "lmroman10-bold", "lmroman10-italic", "lmroman10-bolditalic"):
        fm.fontManager.addfont(os.path.join(LM, nome + ".otf"))
    FONTE = "Latin Modern Roman"
plt.rcParams.update({"axes.formatter.use_locale":False,"font.family":FONTE,"mathtext.fontset":"cm","axes.unicode_minus":False,"font.size":8,"axes.spines.top":False,"axes.spines.right":False,
  "axes.edgecolor":"#888","axes.linewidth":0.6,"xtick.color":"#444","ytick.color":"#444","axes.grid":True,
  "grid.color":"#e3e3e3","grid.linewidth":0.5,"legend.frameon":False,"pdf.fonttype":3})
B,O,A,G="#2a78d6","#eb6834","#1baf7a","#9a9a9a"
rows=list(csv.DictReader(open("../resultados/resultados_45565.csv")))
g=collections.defaultdict(list)
for r in rows: g[(r["experimento"],int(r["np"]),int(r["k"]),r["divisao"])].append(r)
best={k:min(v,key=lambda r:float(r["t_total"])) for k,v in g.items()}
T1=float(best[("escala",1,0,"igual")]["t_total"])
nps=[1,4,8,16,32,64]
def b(n): return best[("escala",n,0,"igual")] if n==1 else best[("escala",n,2,"balanceada")]
S=[T1/float(b(n)["t_total"]) for n in nps]
W=3.3
# 1) convenção das planilhas da aula: colunas = eficiência (eixo da direita), linhas = speedup (esquerda)
import glob
fr=[r for f in sorted(glob.glob("../resultados/resultados_*.csv")) if not f.endswith("45565.csv") for r in csv.DictReader(open(f)) if r["experimento"]=="fraca"]
def combo(ax,ps,S,E,titulo,extra=None):
    x=np.arange(len(ps)); ax2=ax.twinx()
    ax2.bar(x,E,0.55,color="#b7d3f6",zorder=1,label="eficiência (%)")
    pass  # valores da eficiência ficam nas tabelas
    ax2.set_ylim(0,115); ax2.set_ylabel("eficiência (%)",color="#1c5cab"); ax2.grid(False)
    ax2.spines["right"].set_visible(True); ax2.tick_params(axis="y",colors="#1c5cab")
    ax.set_zorder(ax2.get_zorder()+1); ax.patch.set_visible(False)
    ax.plot(x,ps,"--",color=G,lw=1,label="speedup ideal")
    ax.plot(x,S,"-o",color=O,lw=2,ms=4,label="speedup")
    for i,v in enumerate(S): ax.annotate(f"{v:.1f}".replace(".",","),(i,v),textcoords="offset points",xytext=(0,5),ha="center",fontsize=6.5)
    if extra: extra(ax,x)
    ax.set_xticks(x,[str(p) for p in ps]); ax.set_xlabel("processos"); ax.set_ylabel("speedup")
    ax.set_ylim(0,max(ps)*1.05); ax.set_title(titulo,fontsize=8.5)
    h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
    ax.legend(h1+h2,l1+l2,loc="upper center",bbox_to_anchor=(0.5,-0.3),ncol=3,fontsize=6.5)
fig,ax=plt.subplots(1,2,figsize=(W*2,2.05))
E=[100*s/n for s,n in zip(S,nps)]
def ht(a,x): a.axvline(4.5,color=G,lw=0.8,ls=":"); a.text(4.6,55,"HT",color="#666",fontsize=7)
combo(ax[0],nps,S,E,"Escalabilidade forte (mapa inteiro, $T_1$ = 18,56 s)",ht)
if fr:
    gg=collections.defaultdict(list)
    for r in fr: gg[(int(r["linhas"]),int(r["np"]))].append(float(r["t_total"]))
    cargas=sorted({k[0] for k in gg})
    pw=[1,4,8,16,32]; Sw=[];Ew=[]
    for p,L in zip(pw,cargas):
        ts=min(gg[(L,1)]); tp=min(gg[(L,p)]) if p>1 else ts
        Sw.append(ts/tp); Ew.append(100*ts/tp/p)
    combo(ax[1],pw,Sw,Ew,"Escalabilidade fraca (carga = $p$/32 do mapa)")
else:
    ax[1].text(0.5,0.5,"aguardando execução",ha="center",transform=ax[1].transAxes)
fig.tight_layout(); fig.savefig("fig/speedup.pdf")
# 2) decomposicao do tempo
fig,ax=plt.subplots(figsize=(W,1.65))
lab=[str(n) for n in nps[1:]]
med=[];esp=[];tro=[]
for n in nps[1:]:
    r=b(n); tl=float(r["t_local"]); d=float(r["desbal"])/100; mean=tl/(1+d)
    med.append(mean); tro.append(float(r["t_troca"])); esp.append(float(r["t_total"])-mean-float(r["t_troca"]))
x=np.arange(len(lab))
ax.bar(x,med,0.6,color=B,label="cálculo (média)",edgecolor="white",lw=1)
ax.bar(x,esp,0.6,bottom=med,color=O,label="espera (sincronização)",edgecolor="white",lw=1)
ax.bar(x,tro,0.6,bottom=np.add(med,esp),color=A,label="troca de bordas",edgecolor="white",lw=1)
ideal=[T1/n for n in nps[1:]]
ax.plot(x,ideal,"_",color="black",ms=14,mew=1.5,label="$T_1/p$")
ax.set_xticks(x,lab); ax.set_xlabel("processos"); ax.set_ylabel("s"); ax.set_ylim(0,8.2); ax.legend(fontsize=6.5,loc="upper right")
ax.set_title("Onde vai o tempo (balanceada, $k=2$)",fontsize=8.5)
fig.tight_layout(); fig.savefig("fig/tempo.pdf")
# 3) carga por rank
def lecarga():
    d={};cur=None
    for l in open("../resultados/carga_45565.txt"):
        if l.startswith("## "): cur=l.split("=")[1].strip(); d[cur]=[]
        elif l[0].isdigit(): d[cur].append([float(v) for v in l.split(",")])
    return {k:np.array(v) for k,v in d.items()}
c=lecarga()
fig,ax=plt.subplots(figsize=(W*2,1.65))
r=np.arange(32); w=0.4
for off,(k,col,nm) in zip((-w/2,w/2),[("igual",O,"igual (388 linhas cada)"),("balanceada",B,"balanceada (mesmo nº de células com dado)")]):
    t=c[k][:,1]; ax.bar(r+off,t,w,color=col,label=f"{nm}: máx {t.max():.2f} s, média {t.mean():.2f} s".replace(".",","))
    ax.axhline(t.mean(),color=col,lw=0.9,ls="--")
ax.set_xticks(r,[str(i) for i in r],fontsize=6.5); ax.set_xlabel("processo (0 = norte do RS, 31 = Chuí)"); ax.set_ylabel(r"$t_{\mathrm{local}}$ (s)")
ax.set_xlim(-0.7,31.7); ax.set_ylim(0,1.75); ax.legend(fontsize=7,loc="upper center",ncol=2,bbox_to_anchor=(0.5,1.2)); ax.set_title("Tempo de cálculo de cada processo, 32 processos, $k=2$",fontsize=8.5,pad=20)
fig.tight_layout(); fig.savefig("fig/carga.pdf")
# 4) efeito do k
ks=[0,1,2,4,8,16]
def kb(k): return best[("escala",32,2,"balanceada")] if k==2 else best[("k",32,k,"balanceada")]
fig,ax=plt.subplots(1,2,figsize=(W,1.65))
lbl=[r"$\infty$" if k==0 else str(k) for k in ks]
t=[float(kb(k)["t_total"]) for k in ks]; rd=[int(kb(k)["rodadas"]) for k in ks]
ax[0].bar(range(6),t,0.6,color=[O if k==2 else B for k in ks]); ax[0].set_xticks(range(6),lbl)
ax[0].set_title("tempo (s)",fontsize=8); ax[0].set_xlabel("$k$")
for i,v in enumerate(t): ax[0].text(i,v+0.08,f"{v:.1f}".replace(".",","),ha="center",fontsize=6.5)
ax[1].bar(range(6),rd,0.6,color=[O if k==2 else B for k in ks]); ax[1].set_xticks(range(6),lbl)
ax[1].set_title("rodadas",fontsize=8); ax[1].set_xlabel("$k$")
for i,v in enumerate(rd): ax[1].text(i,v+1.5,str(v),ha="center",fontsize=6.5)
for a in ax: a.grid(axis="x",visible=False)
fig.suptitle("Varreduras por rodada ($k$), 32 processos",fontsize=8.5,y=0.98)
fig.tight_layout(); fig.savefig("fig/k.pdf")
print("ok")
