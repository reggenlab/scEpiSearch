import itertools, json, subprocess, sys, os, csv, random, glob, re, gzip, time, matplotlib, multiprocessing
import numpy as np
import pandas as pd
import statistics as sc
import numpy_indexed as npi
from scipy import stats
from scipy.stats import hypergeom
from scipy.stats import spearmanr
import seaborn as sns
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from collections import OrderedDict 
from sklearn.manifold import TSNE
import sklearn.preprocessing
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
from numpy.random import seed
import sklearn as sk
import errno , subprocess
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import SimpSOM.hexagons as hx
import SimpSOM.densityPeak as dp
import SimpSOM.qualityThreshold as qt
from sklearn import cluster
from collections import OrderedDict

def reg_gen_dist(chr_start , chr_end , ref_start , ref_end, strand):
		rs = chr_start + 1
		re = chr_end
		gs = ref_start + 1
		ge = ref_end
		if strand == '+':
			l = [gs - rs , gs - re]
			l.sort()
			return l[0]
		else:
			l = [ge - rs , ge - re]
			l.sort()
			return l[0]

def is_intronic(chr_start , chr_end , ref_start , ref_end):
		rs = chr_start + 1
		re = chr_end
		gs = ref_start + 1
		ge = ref_end
		if (rs > ge or gs>re):
			return 1
		else:
			return 0

def foreground_calc(d_chr , d ,nearest_gene):
		larg_dist = 1e10
		for i in d_chr.keys():
			for val in d_chr[i]:
				cdL = - larg_dist   # left
				cdR =  larg_dist    # right
				cdI = larg_dist     # intronic
				cgL, cgR, cgI = "","",""
				for val2 in d[i]:
					if(is_intronic(val['Start'] , val['End'] , val2['txStart'] , val2['txEnd'])):
						dist = abs(reg_gen_dist(val['Start'] , val['End'] , val2['txStart'] , val2['txEnd'], val2['strand']))
						if (dist < cdI):
							cdI = dist
							cgI = val2['name2']
					else:
						dist = reg_gen_dist(val['Start'] , val['End'] , val2['txStart'] , val2['txEnd'], val2['strand'])
						if (dist <= 0):
							if (dist > cdL):
								cdL = dist
								cgL = val2['name2']
						elif (dist < cdR):
								cdR = dist
								cgR = val2['name2']
				cdist = []
				cdist.extend((abs(cdL) , cdR , cdI))
				cgene = []
				cgene.extend((cgL , cgR , cgI))
				cd_idx = np.argsort(cdist)
				cd = np.array(cdist)[cd_idx[:2]]
				cg = np.array(cgene)[cd_idx[:2]]
				if( cg[1] != ''):
					tmp_g = cgene[cd_idx[2]]
					cg[1] == cg[0] and (tmp_g != '') and cg[1] == tmp_g and cd[1] == cdist[cd_idx[2]]
				else:
					cg[1] = cg[0]
					cd[1] = cd[0]
				nearest_gene.iloc[val['index'],0] = cg[0]
				nearest_gene.iloc[val['index'],1] = cd[0]
				nearest_gene.iloc[val['index'],2] = cg[1]
				nearest_gene.iloc[val['index'],3] = cd[1]
		return nearest_gene

def nearest_gene_accurate(query_type, chr_file, acc_fast, query_file):
    epi = np.loadtxt(query_file,delimiter=",")
#     epi = np.reshape(epi, (epi.shape[0], -1))
#     epi=epi[:,:20]
    if (query_type == 1):
        ref = pd.read_csv('scepisearch_integration/human/refseq-hg19.txt' , sep = '\t')
        ref.loc[:,'chrom'] = (ref['chrom'].str.split("_", expand=True)).iloc[: , 0]
        chr = pd.read_csv(chr_file, sep ='\t', header = None)

        d = OrderedDict()
        for i in ref['chrom'].unique():
            d[i] = [{'name' : ref['name'][j] , 'strand' : ref['strand'][j] , 'txStart' : ref['txStart'][j] ,'txEnd' : ref['txEnd'][j] ,'exonCount' : ref['exonCount'][j] ,'name2': ref['name2'][j]} for j in ref[ref['chrom']==i].index]

        cmd = ['Rscript scepisearch_integration/human/accessibility_score_faster/global_score.R '+chr_file+' scepisearch_integration/acc_score.csv scepisearch_integration/foreground.csv']
        print(cmd)
        process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        ref = pd.read_csv('scepisearch_integration/mouse/refGene.txt' , sep = '\t')
        ref.loc[:,'chrom'] = (ref['chrom'].str.split("_", expand=True)).iloc[: , 0]
        chr = pd.read_csv(chr_file, sep ='\t', header = None)

        d = OrderedDict()
        for i in ref['chrom'].unique():
            d[i] = [{'name' : ref['name'][j] , 'strand' : ref['strand'][j] , 'txStart' : ref['txStart'][j] ,'txEnd' : ref['txEnd'][j] ,'exonCount' : ref['exonCount'][j] ,'name2': ref['name2'][j]} for j in ref[ref['chrom']==i].index]

        cmd = ['Rscript scepisearch_integration/mouse/accessibility_score_faster/global_score.R '+chr_file+' scepisearch_integration/acc_score.csv scepisearch_integration/foreground.csv']
        process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if(process == 0):
        nearest_gene_epi_path = 'scepisearch_integration/foreground.csv'
        nearest_gene = pd.read_csv(nearest_gene_epi_path, sep="\t", header=None)
        if (acc_fast==2):
            acc_score = np.loadtxt("scepisearch_integration/acc_score.csv")
            epi = epi[(nearest_gene!=0).all(axis=1)]
            acc_score = acc_score[(nearest_gene!=0).all(axis=1)]
            nearest_gene = nearest_gene[(nearest_gene!=0).all(axis=1)]
            np.savetxt("scepisearch_integration/acc_score.csv",acc_score)
            return nearest_gene, epi
        else:
            ind_bool = (nearest_gene!=0).all(axis=1)
            ind_zero = np.where(ind_bool)[0]

            d_chr = OrderedDict()
            for i in chr.iloc[:,0].unique():
                d_chr[i] = [{'Start' : chr.iloc[j,1] ,'End' : chr.iloc[j,2] ,'index' : j} for j in chr[chr.iloc[:,0]==i].index if j not in ind_zero]
            nearest_gene = foreground_calc(d_chr , d , nearest_gene)
            nearest_gene.to_csv("scepisearch_integration/foreground.csv",header=False,sep="\t")
            return nearest_gene, epi
    else:
        nearest_gene = pd.DataFrame({})
        return nearest_gene,epi


import multiprocessing, itertools
from scipy.stats import hypergeom
import numpy as np
import pandas as pd
import numpy_indexed as npi

def apply_par(args):
		d, nearest_gene, gene_list = args
		m_gene = d.apply(lambda x : gene_enrich(x , nearest_gene, gene_list), axis = 0)
		return m_gene

def hypergeometric_test(foreground , background):
		fgc1 = list(foreground.iloc[:,0])
		fgd1 = list(foreground.iloc[:,1])
		fgc2 = list(foreground.iloc[:,2])
		fgd2 = list(foreground.iloc[:,3])

		for i in range(len(fgd1)):
			if fgd1[i] <= fgd2[i]:
				fgd2[i] = 2000000
			else:
				fgd1[i] = 2000000

		afc = np.append(fgc1, fgc2)
		afd = np.append(fgd1,fgd2)
		pos = np.where(afd < 1000000)
		afc = afc[pos]

		#---------------------------------
		bgc1 = list(background.iloc[:,0])
		bgd1 = list(background.iloc[:,1])
		bgc2 = list(background.iloc[:,2])
		bgd2 = list(background.iloc[:,3])


		for i in range(len(bgd1)):
			if bgd1[i] <= bgd2[i]:
				bgd2[i] = 2000000
			else:
				bgd1[i] = 2000000

		abc = np.append(bgc1, bgc2)
		abd = np.append(bgd1,bgd2)
		pos = np.where(abd < 1000000)
		abc = abc[pos]

		#---------------------------------
		N = int(len(abc)) + 1
		n = int(len(afc)) + 1

		fgGenes =  set(afc)
		fgGenes = list(fgGenes)
		apvals = list()

		for i in range(len(fgGenes)):
			pos = np.where(afc == fgGenes[i])[0]
			pos1 = np.where(abc == fgGenes[i])[0]
			k = len(pos)
			K = len(pos1)
			mn = min(K,n) +1
			pval = 0

			for j in range(k,mn):
				pval = pval + hypergeom.pmf(j, N, n, K)
			apvals.append(pval)

		summary = np.hstack((np.expand_dims(apvals,axis=1),np.expand_dims(fgGenes, axis=1)))
		summaryDF = pd.DataFrame(summary, columns=["P-value", "Genes"])
		return summaryDF

def gene_enrich(x, nearest_gene, gene_list):
		col_ind = int(x.name)
		fore = np.take(nearest_gene , x , axis = 0)
		gene_enrich = hypergeometric_test(fore, nearest_gene)
		gene_pval = list(gene_enrich['Genes'])
		pval_total = gene_enrich['P-value'].astype(float)
		pval_total = list(pval_total)
		pval_total.append(1)
		pval_total = -np.log2(pval_total)
		index = pval_total.argsort()[-50:][::-1]
		mark_gene = np.asarray(gene_pval)[index]
		ind = npi.indices(gene_pval ,gene_list, missing = -1)
		pval_corr = np.asarray(pval_total)[ind]
		return pd.Series([mark_gene, pval_corr])

def gene_enrichment_calc(epi,gene_list,nearest_gene):
		gene_enriched = np.zeros([epi.shape[1], len(gene_list)])
		marker_gene = np.empty([epi.shape[1], len(range(50))] ,dtype = 'object')

		Z =10000
		sorted_col_idx = np.argsort(epi, axis=0)[epi.shape[0]-Z::,:]

		workers = 10
		pool = multiprocessing.Pool(processes=workers)
		split_dfs = np.array_split(pd.DataFrame(sorted_col_idx), workers, axis=1)

		res = pd.concat(pool.map(apply_par, [(d ,nearest_gene, gene_list) for d in split_dfs]), axis = 1)
		pool.close()
		res_0 = list(itertools.chain.from_iterable(res.iloc[0]))
		res_1 = list(itertools.chain.from_iterable(res.iloc[1]))
		marker_gene = pd.DataFrame(np.array(res_0).reshape(sorted_col_idx.shape[1],50))

		gene_enriched = pd.DataFrame(np.array(res_1).reshape(sorted_col_idx.shape[1], len(gene_list)))
		gene_enriched = np.transpose(np.array(gene_enriched))

		np.savetxt('scepisearch_integration/enrichment_scores.txt', gene_enriched , delimiter=" ", fmt='%f')
		np.savetxt('scepisearch_integration/marker_gene.txt', marker_gene , delimiter=" ", fmt = '%s')
		return gene_enriched
	
	
#!/usr/bin/python
import numpy as np , pandas as pd
from multiprocessing import Process
import numpy_indexed as npi

def get_digit(x):
    return(int(re.search(r'\d+', x).group()))

def process_query(chr_file, epi_path, top_study, query_type, acc_fast,active_poised,imputation):
    if(query_type==1):
        sc_gene_path = 'scepisearch_integration/human/genes_21159.txt'
    else:
        sc_gene_path = 'scepisearch_integration/mouse/gene_mouse.csv'

    gene_list = list()
    with open(sc_gene_path, 'r') as fl:
        for l in fl.readlines():
            gene_list.append(l.rstrip('\n'))

    nearest_gene,epi = nearest_gene_accurate(query_type,chr_file,acc_fast,epi_path)
    # nearest_gene=pd.read_csv("./foreground.csv",header=None,sep=",")
    # epi=np.loadtxt(epi_path,delimiter=",")
    if not nearest_gene.empty:
        acc_score_path = 'scepisearch_integration/acc_score.csv'
        acc_score = np.loadtxt(acc_score_path)
        # epi = np.loadtxt(epi_path, delimiter=",")
        acc_score[acc_score == 0] = 1

        epi = np.array(epi)/(acc_score[:, np.newaxis])

        start_nearest = list(nearest_gene.iloc[:,1])
        end_nearest = list(nearest_gene.iloc[:,3])
        genename_start = list(nearest_gene.iloc[:,0])
        genename_end = list(nearest_gene.iloc[:,2])

        #gene names from nearest genes from either of the two gene locations
        epi_gene = [0] * len(start_nearest)
        for k in range(len(start_nearest)):
            if start_nearest[k] <= end_nearest[k]:
                epi_gene[k] = genename_start[k]
            else:
                epi_gene[k] = genename_end[k]

        gene_enriched = gene_enrichment_calc(epi,gene_list,nearest_gene)
    
    return gene_enriched,epi_gene,epi

def median_calc_null(x , corr_mat, exp_ref):
#     ind_top = [x for x in x if x != -1]
    mat_top = np.take(exp_ref , x , axis = 0)
    med_top = np.median(mat_top , axis=0)
#     print(med_top.shape)
    corr_mat[x.name,:] = med_top
	
def cluster_exp_par(args):
    val,key,epi,final_res,lock,epi_gene,mouse_gene_lower = args
    f0 = gzip.GzipFile('scepisearch_integration/MCA_reference_scepisearch/clust_'+str(key)+'.npy.gz', "r")
    epi_ref = np.load(f0)
#     f3 = gzip.GzipFile('.//clust_'+str(key)+'.npy.gz', "r")
#     corr_mat = np.load(f3)
    
#     for i in range(epi_ref.shape[1]):
#         epi_ref[:,i]=epi_ref[:,i]/np.mean(epi_ref[:,i])
    
    sorted_col = np.loadtxt("scepisearch_integration/gene_name_exp_loc_mouse.txt",delimiter=",")
    sorted_col=sorted_col.astype(int)
    
    clust_infor = pd.read_csv("scepisearch_integration/MCA_reference_scepisearch/clusters_final.txt",sep=" ",dtype='str')
#     metadata = pd.read_csv("./searchProject/meta_mouse/metadata_exp.csv", header=None,sep="@")
#     metadata = metadata.ix[:,2]
#     nan_rows = metadata[metadata.isnull()]
#     unknown = metadata[metadata == 'Unknown']
#     white_cells = metadata[metadata == 'White blood cells']
#     unanno =  np.hstack((pd.Series(nan_rows.index) , pd.Series(unknown.index)))
#     unanno = np.hstack((pd.Series(unanno),pd.Series(white_cells.index)))
#     print(len(unanno))
    
    ind_ref = np.array(clust_infor.loc[clust_infor['cluster_no'] == str(key), 'id'])
    
#     epi_ref = epi_ref[:,~np.isin(ind_ref.astype(int),unanno.astype(int))]
#     ind_ref = ind_ref[~np.isin(ind_ref.astype(int),unanno.astype(int))]

    query = np.take(epi, list(val), axis=1)
    print(query.shape)
    
    #pick top 1000 sites from epi data
    N = 500
    sorted_col_idx_epi = np.argsort(query, axis=0)[query.shape[0]-N::,:]
    top_epi_gene = pd.DataFrame(sorted_col_idx_epi).apply(lambda x: np.take(epi_gene, indices = x) , axis = 0)
    gene_name_exp_loc = top_epi_gene.apply(lambda x : npi.indices(mouse_gene_lower ,x, missing = -1)  , axis = 0)
    
    res = gene_name_exp_loc.apply(lambda x : median_calc(x , epi_ref) ,axis = 0)
    res = np.array(res)
    
    corr_mat = np.zeros([1000,epi_ref.shape[1]])
    pd.DataFrame(sorted_col).apply(lambda x : median_calc_null(x , corr_mat ,epi_ref) ,axis = 0)
    corr_mat = np.array(corr_mat)
      
    pval = []
    sorted_10_idx = np.argsort(res, axis=0)[res.shape[0]-50::,:]
    sorted_raveled = sorted_10_idx.ravel()
    col_idx = np.arange(res.shape[1])
    val_top = res[sorted_10_idx, col_idx]
    val_top_raveled = val_top.ravel()
    for a, b in zip(sorted_raveled, val_top_raveled):
        pval.append(sum(corr_mat[:,a] >= b))
    pval = np.reshape(pval , (sorted_10_idx.shape[0] , sorted_10_idx.shape[1]))
    
    pval = pval/float(1000)
    pval = pval + 0.001

    p_adjust_epi = stats.p_adjust(FloatVector(pval.ravel()), method = 'BH')
    p_adjust_epi = np.array(p_adjust_epi).reshape(pval.shape[0], pval.shape[1])

    result = np.array(ind_ref)[sorted_10_idx.ravel()]
    result = result.reshape(sorted_10_idx.shape[0],sorted_10_idx.shape[1])
    
    with lock:
        for i,j in enumerate(val):
            final_res[j,'corr'] = np.append(final_res[j,'corr'], val_top[:,i])
            final_res[j,'index'] = np.append(final_res[j,'index'], result[:,i])
            final_res[j,'pval'] = np.append(final_res[j,'pval'], pval[:,i])
            final_res[j,'adj_pval'] = np.append(final_res[j,'adj_pval'],p_adjust_epi[:,i])

##########################Select query peak file and count file of query cells #################################
chr_file='queries_scepisearch/neuron_GSM2579603_peaks_hg19_final.bed'
query_file='queries_scepisearch/neuron_GSE97942_SCEPI.csv'
top_study=5
##############################select query_type=1 if human else 2
query_type=1
acc_fast=2
active_poised=2
imputation=1

gene_enriched,epi_gene,epi = process_query(chr_file,query_file,top_study,query_type,acc_fast,active_poised,imputation)
print(epi.shape)

gene_mouse = list()
with open("scepisearch_integration/mouse/gene_mouse.csv", 'r') as fl:
    for l in fl.readlines():
        gene_mouse.append(l.rstrip('\n'))
	
sc_gene_path = 'scepisearch_integration/human/genes_21159.txt'
gene_list = list()
with open(sc_gene_path, 'r') as fl:
    for l in fl.readlines():
        gene_list.append(l.rstrip('\n'))
	
human_gene_lower = [x.lower() for x in gene_list]
mouse_gene_lower = [x.lower() for x in gene_mouse]
intersect_gene = list(set(human_gene_lower) & set(mouse_gene_lower))

gene_enriched_mouse = gene_enriched

data = np.loadtxt(query_file,delimiter=",")
cols = epi.shape[1]
gene_enriched_mouse = np.zeros((23595,cols)) 
for i in range(len(human_gene_lower)):
    if human_gene_lower[i] in mouse_gene_lower:
        ind = mouse_gene_lower.index(human_gene_lower[i])
        gene_enriched_mouse[ind,:] = gene_enriched[i,:]
    else:
        continue
	
mean_array = np.load("scepisearch_integration/MCA_reference_scepisearch/mean_array.npy")
gene_enriched_mouse = np.array(gene_enriched_mouse)

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
        
nr, nc = gene_enriched_mouse.shape
xvec = ro.FloatVector(gene_enriched_mouse.transpose().reshape((gene_enriched_mouse.size)))
xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)
    
#     print(np.array(xr).shape)

nry, ncy = mean_array.shape
xvecy = ro.FloatVector(mean_array.transpose().reshape((mean_array.size)))
yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)
    
#     print(np.array(yr).shape)

res = ro.r.cor(xr,yr, method="spearman")
res = np.array(res)
        
res = np.transpose(res)

index_value = np.genfromtxt("scepisearch_integration/MCA_reference_scepisearch/mean_array_subclusterindex.txt",dtype='str')

N = 10
sorted_clust = np.argsort(res, axis=0)[res.shape[0]-N::,:]
sorted_clust = index_value[sorted_clust.ravel()].reshape((sorted_clust.shape[0],sorted_clust.shape[1]))

gene_enriched_mouse = pd.DataFrame(gene_enriched_mouse)
gene_enriched_mouse.index = mouse_gene_lower
epi_gene = [x.lower() for x in epi_gene]

def median_calc(x , epi_ref):
    ind_top = [x for x in x if x != -1]
    mat_top = np.take(epi_ref , ind_top , axis = 0)
    med_top = np.median(mat_top , axis=0)
    return med_top

un = np.unique(sorted_clust)
id_clust = OrderedDict()
for i,j in enumerate(un):
    id_clust[str(j)]=[]
    id_clust[str(j)][:] = np.where(sorted_clust == j)[1]
	
# final_res = OrderedDict()
manager = multiprocessing.Manager()
final_res = manager.dict()
lock = manager.Lock()
for i in range(cols):
#     final_res[i] = {}
    final_res[i,'corr'] = np.array([])
    final_res[i,'index'] = np.array([])
    final_res[i,'pval'] = np.array([])
    final_res[i,'adj_pval'] = np.array([])

workers= 30
p = multiprocessing.Pool(processes=workers)                     #number of processes = number of CPUs
keys, values= zip(*id_clust.items())

p.map(cluster_exp_par, [(val,key,epi,final_res,lock,epi_gene,mouse_gene_lower) for val,key in zip(values,keys)])
p.close()


# correlation_matrix = np.zeros([10,81173])
metadata = pd.read_csv("scepisearch_integration/MCA_reference_labels.csv", header=None, sep="@")
final_corr = np.zeros([cols,5], dtype='int')
pval_epi = np.zeros([cols,5])
final_fdr = np.zeros([cols,5])

for i in range(268):
    ind = np.argsort(final_res[i,'corr'])[::-1][:5]
#     correlation_matrix[i,res[i,'index'].astype(int)] = res[i,'corr']
    final_corr[i,:] = final_res[i,'index'][ind].astype(int)
    pval_epi[i,:] = final_res[i,'pval'][ind]
    final_fdr[i,:] = final_res[i,'adj_pval'][ind]
	
# metadata = metadata.iloc[:,0]
final_corr = final_corr.astype(int)
cells = metadata[np.ravel(final_corr)]
cells_forebrain = pd.DataFrame(np.reshape(np.array(cells), (cols,5)))

