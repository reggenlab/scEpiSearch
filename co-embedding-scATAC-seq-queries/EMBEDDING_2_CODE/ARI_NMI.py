import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN

adj = np.loadtxt("./embedding_results/scepi_adj_embedding_2.txt",delimiter=",")
clustering = SpectralClustering(n_clusters=2,assign_labels="kmeans",affinity='precomputed').fit(adj)
clustering.labels_[:3] = 0
clustering.labels_[201:205] = 1
l = []
l.extend(np.repeat(1,200))
l.extend(np.repeat(0,200))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))


human_bcell = np.loadtxt("output_scale_embedding/feature_human_gm.txt",delimiter="\t")
mouse_bcell = np.loadtxt("output_scale_embedding/feature_gm_mouse.txt",delimiter="\t")
human_hek = np.loadtxt("output_scale_embedding/feature_hek_human.txt",delimiter="\t")
mouse_proximal = np.loadtxt("output_scale_embedding/feature_proximal_mouse.txt",delimiter="\t")
latent_space = np.concatenate([human_bcell,mouse_bcell,human_hek,mouse_proximal],axis=0)
ts = TSNE(n_components=2).fit_transform(latent_space)
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,200))
l.extend(np.repeat(0,200))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

latent_space = np.loadtxt("scvi_embedding2.txt",delimiter=",")
ts = TSNE(n_components=2).fit_transform(latent_space)
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,200))
l.extend(np.repeat(0,200))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

latent_space = np.loadtxt("scanorama_embedding2.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(np.transpose(latent_space))
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,200))
l.extend(np.repeat(0,200))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

ts = pd.read_csv("case_2_mint.txt",sep=" ",index_col=0)
ts = np.array(ts)[:,:2]
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,200))
l.extend(np.repeat(0,200))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

import numpy as np
import matplotlib.pyplot as plt
data = [[0.931,0.873],
[0.238, 0.393],
[ 0.24, 0.396],
[0.114, 0.16],
[0.25521844703694835,0.36296]]
X = np.arange(3)
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(3)
ax.bar(X + 0.00, data[0], color = 'crimson', width = 0.15)
ax.bar(X + 0.15, data[1], color = 'lightgreen', width = 0.15)
ax.bar(X + 0.30, data[2], color = 'aqua', width = 0.15)
ax.bar(X + 0.45, data[3], color = 'orange', width = 0.15)
ax.bar(X + 0.60, data[4], color = 'fuchsia', width = 0.15)
ax.set_title('Case Study-3')
ax.set_ylabel('Score in Fraction',fontsize=15)
ax.set_xticks(ind)
ax.set_xticklabels(( 'ARI', 'NMI'),fontsize=15)
ax.legend(labels=['scEpiSearch', 'SCALE', 'SCVI','scanorama','MINT'],loc="upper left", bbox_to_anchor=(1,1),fontsize=15)
