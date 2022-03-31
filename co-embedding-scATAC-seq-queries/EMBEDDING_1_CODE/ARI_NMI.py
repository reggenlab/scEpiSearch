import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN

adj = np.loadtxt("scepi_adj_embedding_1.txt",delimiter=",")
clustering = SpectralClustering(n_clusters=4,assign_labels="kmeans",random_state=0,affinity='precomputed').fit(adj)
l = []
l.extend(np.repeat(0,30))
l.extend(np.repeat(0,30))
l.extend(np.repeat(1,12))
l.extend(np.repeat(1,11))
l.extend(np.repeat(3,14))
l.extend(np.repeat(3,10))
l.extend(np.repeat(3,10))
l.extend(np.repeat(4,10))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

import matplotlib.pyplot
from sklearn.manifold import TSNE
ts_mouse_hsc = np.loadtxt("output_mouse_hsc/feature.txt",delimiter="\t")
ts_human_hsc = np.loadtxt("output_human_hsc_SCALE/feature.txt",delimiter="\t")
ts_mouse_neuron = np.loadtxt("output_mouse_forebrain_SCALE/feature.txt",delimiter="\t")
ts_human_neuron = np.loadtxt("output_human_neuron/feature.txt",delimiter="\t")
ts_mouse_gm = np.loadtxt("output_bcell_mouse_SCALE/feature.txt",delimiter="\t")
ts_human_gm = np.loadtxt("output_human_GM_SCALE/feature.txt",delimiter="\t")
ts_human_gm_gse68103 = np.loadtxt("output_human_GM_GSE68103_SCALE/feature.txt",delimiter="\t")
ts_human_myoblast = np.loadtxt("output_human_myoblast_SCALE/feature.txt",delimiter="\t")
latent_space = np.concatenate([ts_human_neuron,ts_mouse_neuron,ts_human_hsc,ts_mouse_hsc,ts_mouse_gm,ts_human_gm,ts_human_gm_gse68103,ts_human_myoblast],axis=0)
ts = TSNE(n_components=2).fit_transform(latent_space)
import numpy as np
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(0,30))
l.extend(np.repeat(0,30))
l.extend(np.repeat(1,12))
l.extend(np.repeat(1,11))
l.extend(np.repeat(3,14))
l.extend(np.repeat(3,10))
l.extend(np.repeat(3,10))
l.extend(np.repeat(4,10))
l = l[:125]
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))


latent_space = np.loadtxt("scvi_latent_Space.txt",delimiter=",")
ts = TSNE(n_components=2).fit_transform(latent_space)
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(0,30))
l.extend(np.repeat(0,30))
l.extend(np.repeat(1,12))
l.extend(np.repeat(1,11))
l.extend(np.repeat(3,14))
l.extend(np.repeat(3,10))
l.extend(np.repeat(3,10))
l.extend(np.repeat(4,10))
clustering.labels_ = clustering.labels_[:127]
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

latent_space = np.loadtxt("scanorama.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(np.transpose(latent_space))
l = []
l.extend(np.repeat(0,30))
l.extend(np.repeat(0,30))
l.extend(np.repeat(1,12))
l.extend(np.repeat(1,11))
l.extend(np.repeat(3,14))
l.extend(np.repeat(3,10))
l.extend(np.repeat(3,10))
l.extend(np.repeat(4,10))
clustering = DBSCAN(eps=7, min_samples=2).fit(ts)
clustering.labels_ = clustering.labels_[:127]
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
ts = pd.read_csv("case_1_mint.txt",sep=" ",index_col=0)
ts = np.array(ts)[:,:2]
l = []
l.extend(np.repeat(0,30))
l.extend(np.repeat(0,30))
l.extend(np.repeat(1,12))
l.extend(np.repeat(1,11))
l.extend(np.repeat(3,14))
l.extend(np.repeat(3,10))
l.extend(np.repeat(3,10))
l.extend(np.repeat(4,10))
clustering = DBSCAN(eps=3, min_samples=2).fit(ts)
clustering.labels_ = clustering.labels_[:127]
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

import numpy as np
import matplotlib.pyplot as plt
data = [[ 0.90 ,0.8667],
[0.026,  0.217],
[0.401 , 0.6240 ],
[ 0.13, 0.250],
[0.05590, 0.3363]]
X = np.arange(3)
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(3)
ax.bar(X + 0.00, data[0], color = 'crimson', width = 0.15)
ax.bar(X + 0.15, data[1], color = 'lightgreen', width = 0.15)
ax.bar(X + 0.30, data[2], color = 'aqua', width = 0.15)
ax.bar(X + 0.45, data[3], color = 'orange', width = 0.15)
ax.bar(X + 0.60, data[4], color = 'fuchsia', width = 0.15)
ax.set_ylabel('Score in Fraction',fontsize=15)
ax.set_title("Case Study-1")
ax.set_xticks(ind)
ax.set_xticklabels(('Purity', 'ARI', 'NMI'),fontsize=15)
ax.legend(labels=['scEpiSearch', 'SCALE', 'SCVI','scanorama','MINT'],loc="upper left", bbox_to_anchor=(1,1),fontsize=15)

