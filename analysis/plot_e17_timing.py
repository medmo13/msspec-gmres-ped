import csv,statistics as st,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
r=list(csv.DictReader(open('runs/e17_timing/rows.csv')))
r+=[dict(x,case='B01') for x in csv.DictReader(open('runs/e17_big/rows.csv')) if x['rep']!='w']   # LiF: one run per method
names={'H42':'Ru hcp\n0001, 26 atoms','H17':'PbB rocksalt\n011, 48 atoms','H28':'TiO$_2$ rutile\n85 atoms','H48':'Hf hcp\nsemipolar\n106 atoms','B01':'LiF rocksalt\nm201, 277 atoms\n(1 run)'}
M=[('lu','dense LU','#555'),('gmres','tree-GMRES','#1f77b4'),('series3','series order 3','#e69f00')]
f,ax=plt.subplots(1,2,figsize=(10.5,4))
for k,key,lab in [(0,'wall_s','Wall time of the whole 19-energy scan (s)'),(1,'peak_rss_mib','Peak memory (MiB)')]:
    for j,(m,l,c) in enumerate(M):
        v=[st.median(float(x[key]) for x in r if x['case']==cs and x['method']==m and x['rep']!='w') for cs in names]
        ax[k].bar(np.arange(5)+(j-1)*0.27,v,0.27,color=c,label=l)
    ax[k].set_xticks(range(5)); ax[k].set_xticklabels(names.values(),fontsize=6.2); ax[k].set_yscale('log'); ax[k].set_ylabel(lab,fontsize=8)
    ax[k].spines[['top','right']].set_visible(False)
ax[0].legend(fontsize=7,frameon=False,loc='upper left')
with open('analysis/out/e17_timing.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['case','natoms','method','n_runs','wall_median_s','wall_min_s','wall_max_s','peak_rss_mib'])
    for cs in names:
        for m,_,_ in M:
            t=[x for x in r if x['case']==cs and x['method']==m and x['rep']!='w']; wl=[float(x['wall_s']) for x in t]
            w.writerow([cs,t[0]['natoms'],m,len(t),st.median(wl),min(wl),max(wl),st.median(float(x['peak_rss_mib']) for x in t)])
plt.tight_layout(); plt.savefig('figures/e17_timing.png',dpi=200); plt.savefig('figures/e17_timing.pdf')
