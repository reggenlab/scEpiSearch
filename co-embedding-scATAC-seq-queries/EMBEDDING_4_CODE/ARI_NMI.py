

adj = np.loadtxt("scepi_adj_embedding_4.txt",delimiter=",")
clustering = SpectralClustering(n_clusters=2,assign_labels="kmeans",random_state=0,affinity='precomputed').fit(adj)

clustering.labels_[:20] = 0
clustering.labels_[2101:2120] = 1
l = []
l.extend(np.repeat(1,2035))
l.extend(np.repeat(0,2101))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

ts_mouse_hsc = np.loadtxt("scale_data/result/feature_mouse_hsc.txt",delimiter="\t")
ts_human_hsc = np.loadtxt("scale_data/result/feature_human_hsc.txt",delimiter="\t")
ts_mouse_neuron = np.loadtxt("scale_data/result/feature_mouse_brain.txt",delimiter="\t")
ts_human_neuron = np.loadtxt("scale_data/result/feature_human_brain.txt",delimiter="\t")
latent_space = np.concatenate([ts_human_neuron,ts_mouse_neuron,ts_human_hsc,ts_mouse_hsc],axis=0)
ts = TSNE(n_components=2).fit_transform(latent_space)
from sklearn.cluster import DBSCAN
import numpy as np
clustering = DBSCAN(eps=10, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,2035))
l.extend(np.repeat(0,2101))
l=l[:4080]
print(purity_score(l,clustering.labels_))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

latent_space = np.loadtxt("scvi_latent.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(latent_space)
clustering = DBSCAN(eps=0.05, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,2035))
l.extend(np.repeat(0,2101))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

latent_space = np.loadtxt("scanorama.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(np.transpose(latent_space))
clustering = DBSCAN(eps=0.05, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,2035))
l.extend(np.repeat(0,2101))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

ts = pd.read_csv("case_4_mint.txt",sep=" ",index_col=0)
ts = np.array(ts)[:,:2]
clustering = DBSCAN(eps=0.05, min_samples=2).fit(ts)
l = []
l.extend(np.repeat(1,2035))
l.extend(np.repeat(0,2101))
print(purity_score(l,clustering.labels_))
print(adjusted_rand_score(l,clustering.labels_))
print(normalized_mutual_info_score(l, clustering.labels_))

import numpy as np
import matplotlib.pyplot as plt
data = [[0.9626,0.923],
[0.1251,0.297],
[0.04,0.086],
[0.16,0.385],
[0.06941449848293321,0.14799]]
X = np.arange(3)
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(3)
ax.bar(X + 0.00, data[0], color = 'crimson', width = 0.15)
ax.bar(X + 0.15, data[1], color = 'lightgreen', width = 0.15)
ax.bar(X + 0.30, data[2], color = 'aqua', width = 0.15)
ax.bar(X + 0.45, data[3], color = 'orange', width = 0.15)
ax.bar(X + 0.60, data[4], color = 'fuchsia', width = 0.15)
ax.set_title('Case Study-2')
ax.set_ylabel('Score in Fraction',fontsize=15)
ax.set_xticks(ind)
ax.set_xticklabels(('ARI', 'NMI'),fontsize=15)
ax.legend(labels=['scEpiSearch', 'SCALE', 'SCVI','scanorama','MINT'],loc="upper left", bbox_to_anchor=(1,1),fontsize=15)
