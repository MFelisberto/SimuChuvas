import csv, collections
import numpy as np, matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.rcParams.update({"axes.formatter.use_locale":False,"font.family":"serif","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
  "axes.edgecolor":"#888","axes.linewidth":0.6,"xtick.color":"#444","ytick.color":"#444","axes.grid":True,
  "grid.color":"#e3e3e3","grid.linewidth":0.5,"legend.frameon":False,"pdf.fonttype":42})
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
# 1) speedup + eficiencia
fig,ax=plt.subplots(1,2,figsize=(W*2,1.9))
ax[0].plot(nps,nps,"--",color=G,lw=1,label="ideal")
ax[0].plot(nps,S,"-o",color=B,lw=2,ms=4,label="balanceada, $k=2$")
ax[0].plot([32,64],[T1/float(best[("divisao",n,2,"igual")]["t_total"]) for n in (32,64)],"s",color=O,ms=4,label="igual, $k=2$")
for n,s in zip(nps,S):
    if n>1: ax[0].annotate(f"{s:.1f}".replace(".",","),(n,s),textcoords="offset points",xytext=(0,6),ha="center",fontsize=7)
ax[0].set_xscale("log",base=2); ax[0].set_yscale("log",base=2)
ax[0].set_xticks(nps,[str(n) for n in nps]); ax[0].set_yticks([1,2,4,8,16,32,64],["1","2","4","8","16","32","64"])
ax[0].set_ylim(0.8,80); ax[0].set_xlabel("processos"); ax[0].set_title("Speedup ($T_1$ = 18,56 s)",fontsize=8.5)
ax[0].legend(loc="upper left",fontsize=7)
E=[100*s/n for s,n in zip(S,nps)]
trab=[100*27.71/float(b(n)["trabalho"]) for n in nps]
ax[1].plot(nps,E,"-o",color=B,lw=2,ms=4,label="eficiência $E=S/p$")
ax[1].plot(nps,trab,"-^",color=A,lw=2,ms=4,label="eficiência de trabalho")
ax[1].set_xscale("log",base=2); ax[1].set_xticks(nps,[str(n) for n in nps]); ax[1].set_ylim(0,105)
ax[1].set_ylabel("%"); ax[1].set_xlabel("processos"); ax[1].legend(loc="lower left",fontsize=7)
ax[1].set_title("Eficiência",fontsize=8.5)
for a in ax:
    a.axvline(32,color=G,lw=0.8,ls=":"); a.text(45,a.get_ylim()[1]*(0.92 if a is ax[1] else 0.55),"HT",color="#666",fontsize=7,ha="center")
fig.tight_layout(); fig.savefig("fig/speedup.pdf")
# 2) decomposicao do tempo
fig,ax=plt.subplots(figsize=(W,1.9))
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
ax.set_xticks(x,lab); ax.set_xlabel("processos"); ax.set_ylabel("s"); ax.legend(fontsize=6.5,loc="upper right")
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
fig,ax=plt.subplots(figsize=(W*2,1.85))
r=np.arange(32); w=0.4
for off,(k,col,nm) in zip((-w/2,w/2),[("igual",O,"igual (388 linhas cada)"),("balanceada",B,"balanceada (mesmo nº de células com dado)")]):
    t=c[k][:,1]; ax.bar(r+off,t,w,color=col,label=f"{nm}: máx {t.max():.2f} s, média {t.mean():.2f} s".replace(".",","))
    ax.axhline(t.mean(),color=col,lw=0.9,ls="--")
ax.set_xticks(r,[str(i) for i in r],fontsize=6.5); ax.set_xlabel("processo (0 = norte do RS, 31 = Chuí)"); ax.set_ylabel("$t_{local}$ (s)")
ax.set_xlim(-0.7,31.7); ax.set_ylim(0,1.75); ax.legend(fontsize=7,loc="upper center",ncol=2,bbox_to_anchor=(0.5,1.2)); ax.set_title("Tempo de cálculo de cada processo, 32 processos, $k=2$",fontsize=8.5,pad=20)
fig.tight_layout(); fig.savefig("fig/carga.pdf")
# 4) efeito do k
ks=[0,1,2,4,8,16]
def kb(k): return best[("escala",32,2,"balanceada")] if k==2 else best[("k",32,k,"balanceada")]
fig,ax=plt.subplots(1,2,figsize=(W,1.7))
lbl=["∞" if k==0 else str(k) for k in ks]
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
