import numpy as np
import pandas as pd
from collections import OrderedDict
import subprocess,pickle

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
#     epi = pd.DataFrame(epi)
#     epi = epi.drop(epi.columns[[48,97,74,28,90,22,12,6]], axis=1)
#     epi = epi.iloc[:,69:]
#     print(epi.shape)
#     epi = np.array(epi)
    if (query_type == 1):
        ref = pd.read_csv('./searchProject/storage/scepisearch/human/refseq-hg19.txt' , sep = '\t')
        ref.loc[:,'chrom'] = (ref['chrom'].str.split("_", expand=True)).iloc[: , 0]
        chr = pd.read_csv(chr_file, sep ='\t', header = None)

        d = OrderedDict()
        for i in ref['chrom'].unique():
            d[i] = [{'name' : ref['name'][j] , 'strand' : ref['strand'][j] , 'txStart' : ref['txStart'][j] ,'txEnd' : ref['txEnd'][j] ,'exonCount' : ref['exonCount'][j] ,'name2': ref['name2'][j]} for j in ref[ref['chrom']==i].index]

        cmd = ['Rscript ./human/accessibility_score_faster/global_score.R '+chr_file+' ./acc_score.csv ./foreground.csv']
        process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        ref = pd.read_csv('./searchProject/storage/scepisearch/mouse/refGene.txt' , sep = '\t')
        ref.loc[:,'chrom'] = (ref['chrom'].str.split("_", expand=True)).iloc[: , 0]
        chr = pd.read_csv(chr_file, sep ='\t', header = None)

        d = OrderedDict()
        for i in ref['chrom'].unique():
            d[i] = [{'name' : ref['name'][j] , 'strand' : ref['strand'][j] , 'txStart' : ref['txStart'][j] ,'txEnd' : ref['txEnd'][j] ,'exonCount' : ref['exonCount'][j] ,'name2': ref['name2'][j]} for j in ref[ref['chrom']==i].index]

        cmd = ['Rscript ./mouse/accessibility_score_faster/global_score.R '+chr_file+' ./acc_score.csv ./foreground.csv']
        process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if(process == 0):
        nearest_gene_epi_path = './foreground.csv'
        nearest_gene = pd.read_csv(nearest_gene_epi_path, sep="\t", header=None)
        if (acc_fast==2):
            acc_score = np.loadtxt("./acc_score.csv")
            epi = epi[(nearest_gene!=0).all(axis=1)]
            acc_score = acc_score[(nearest_gene!=0).all(axis=1)]
            nearest_gene = nearest_gene[(nearest_gene!=0).all(axis=1)]
            np.savetxt("./acc_score.csv",acc_score)
            return nearest_gene, epi
        else:
            ind_bool = (nearest_gene!=0).all(axis=1)
            ind_zero = np.where(ind_bool)[0]
#             print(len(ind_zero))
#             print(nearest_gene.shape)

            d_chr = OrderedDict()
            for i in chr.iloc[:,0].unique():
                d_chr[i] = [{'Start' : chr.iloc[j,1] ,'End' : chr.iloc[j,2] ,'index' : j} for j in chr[chr.iloc[:,0]==i].index if j not in ind_zero]
            nearest_gene = foreground_calc(d_chr , d , nearest_gene)
            nearest_gene.to_csv("./foreground.csv",header=False,sep="\t")
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

def gene_enrichment_calc(epi,gene_list,nearest_gene,query_type):
        gene_enriched = np.zeros([epi.shape[1], len(gene_list)])
        marker_gene = np.empty([epi.shape[1], len(range(50))] ,dtype = 'object')

        Z =5000
        sorted_col_idx = np.argsort(epi, axis=0)[epi.shape[0]-Z::,:]

        workers = 30
        pool = multiprocessing.Pool(processes=workers)
        split_dfs = np.array_split(pd.DataFrame(sorted_col_idx), workers, axis=1)

        res = pd.concat(pool.map(apply_par, [(d ,nearest_gene, gene_list) for d in split_dfs]), axis = 1)
        pool.close()
        res_0 = list(itertools.chain.from_iterable(res.iloc[0]))
        res_1 = list(itertools.chain.from_iterable(res.iloc[1]))
        marker_gene = pd.DataFrame(np.array(res_0).reshape(sorted_col_idx.shape[1],50))
        
        if query_type == 1:
            a_file = open("./searchProject/meta_human/markers_human.pkl", "rb")
            dict_markers = pickle.load(a_file)
            a_file.close()
        else:
            a_file = open("./searchProject/meta_mouse/markers_mouse.pkl", "rb")
            dict_markers = pickle.load(a_file)
            a_file.close()
        celltype_markers = pd.DataFrame(index = range(marker_gene.shape[0]), columns = range(marker_gene.shape[1]))
        
        for i in range(marker_gene.shape[0]):
            marker_gene.iloc[i,:] = marker_gene.iloc[i,:].str.upper()
            res = marker_gene.iloc[i,:].map(dict_markers)
            celltype_markers.iloc[i,:] = res
        
        a_f = open("./MGI_database.pkl", "rb")
        dict_mgi = pickle.load(a_f)
        a_f.close()
        gene_mgi = pd.DataFrame(index = range(marker_gene.shape[0]), columns = range(marker_gene.shape[1]))
        
        for i in range(marker_gene.shape[0]):
            marker_gene.iloc[i,:] = marker_gene.iloc[i,:].str.upper()
            r = marker_gene.iloc[i,:].map(dict_mgi)
#             print(r)
            gene_mgi.iloc[i,:] = r
        
        
        for i in range(marker_gene.shape[0]):
            marker_gene.iloc[i,:] = marker_gene.iloc[i,:].str.upper()
            res = marker_gene.iloc[i,:].map(dict_markers)
            celltype_markers.iloc[i,:] = res
        
        
        gene_enriched = pd.DataFrame(np.array(res_1).reshape(sorted_col_idx.shape[1], len(gene_list)))
        gene_enriched = np.transpose(np.array(gene_enriched))
        
        celltype_markers.to_csv('./scepisearch_query_results/celltype_markers.txt',sep="\t",header=False,index=False)
        gene_mgi.to_csv('./scepisearch_query_results/gene_mgi.txt',sep="\t",header=False,index=False)
        np.savetxt('./scepisearch_query_results/enrichment_scores.txt', gene_enriched , delimiter=" ", fmt='%f')
        np.savetxt('./scepisearch_query_results/marker_gene.txt', marker_gene , delimiter=" ", fmt = '%s')
        return gene_enriched
		
		
		
import sklearn.preprocessing, multiprocessing, gzip, matplotlib
import matplotlib.pyplot as plt  , seaborn as sns , SimpSOM as sps, matplotlib.pyplot as plt
from collections import OrderedDict
from sklearn import cluster
# from keras.models import model_from_json
# from keras.models import Model
from tensorflow.keras import optimizers
sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import numpy as np
import pandas as pd
import numpy_indexed as npi
matplotlib.use('Agg')
plt.switch_backend('agg')
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

def rand_jitter(arr):
		stdev = .1*(max(arr)-min(arr))
		return arr + np.random.randn(len(arr)) * stdev

def epi_matching(epi, query_type, gene_list, nearest_gene, top_study, epi_gene, gene_enriched, gene_name_exp_loc):
        pinf = float('+inf')
        ninf = float('-inf')
        gene_enriched[gene_enriched == ninf] = 0.0
        gene_enriched[gene_enriched == pinf] = np.finfo(np.float16).max
        
        where_are_NaNs = np.isnan(gene_enriched)
        gene_enriched[where_are_NaNs] = 0
               
        if(query_type == 1):
            correlation_matrix = np.zeros([epi.shape[1],742297])
        else:
            correlation_matrix = np.zeros([epi.shape[1],81173])

#autoencoder output
#         reduced_gene_enriched = autoencoder_output(query_type, gene_enriched)
        res = clusters_corr_epi(query_type,gene_enriched)

# N = 10
        sorted_clust = np.argsort(res, axis=0)[res.shape[0]-50::,:]
    
        if query_type==1:
            index_value = np.genfromtxt("./searchProject/storage/scepisearch/human/epi_new/mean_array_subclusterindex.txt",dtype='str')
            #print(index_value)
            sorted_clust = index_value[sorted_clust.ravel()].reshape((sorted_clust.shape[0],sorted_clust.shape[1]))
            print(sorted_clust.shape)
        
#         sorted_clust = [item for sublist in sorted_clust for item in sublist]
#         print(sorted_clust)
#         arr = np.arange(0,800)
#         arr = arr.astype('str')
#         sorted_clust = [x for x in sorted_clust if str(x) not in arr]
#         print(sorted_clust)
#         sorted_clust = np.reshape(sorted_clust, (len(list(sorted_clust)), -1))
        
        un = np.unique(sorted_clust)
        id_clust = OrderedDict()
        for i,j in enumerate(un):
            id_clust[str(j)]=[]
            id_clust[str(j)][:] = np.where(sorted_clust == j)[1]
#         print(ind_clust)

        workers= 20
        p = multiprocessing.Pool(processes=workers)
        keys, values= zip(*id_clust.items())

        manager = multiprocessing.Manager()
        final_res = manager.dict()
        lock = manager.Lock()

        for i in range(epi.shape[1]):
            final_res[i,'corr'] = np.array([])
            final_res[i,'index'] = np.array([])
            final_res[i,'pval'] = np.array([])
            final_res[i,'adj_pval'] = np.array([])

        p.map(cluster_epi_par, [(val,key,final_res,lock,query_type,gene_enriched) for val,key in zip(values,keys)])
        p.close()

        final_corr = np.zeros([epi.shape[1],top_study], dtype='int')
        pval_epi = np.zeros([epi.shape[1],top_study])
        final_fdr = np.zeros([epi.shape[1],top_study])

        for i in range(epi.shape[1]):
            ind = np.argsort(final_res[i,'corr'])[::-1][:top_study]
            correlation_matrix[i,final_res[i,'index'].astype(int)] = final_res[i,'corr']
            final_corr[i,:] = final_res[i,'index'][ind].astype(int)
            pval_epi[i,:] = final_res[i,'pval'][ind]
            final_fdr[i,:] = final_res[i,'adj_pval'][ind]
        
        ind_corrmat_nonzero = ~np.all(correlation_matrix == 0, axis=0)
        ind_corrmat_nonzero = np.where(ind_corrmat_nonzero)[0]
        correlation_matrix = correlation_matrix[:,~np.all(correlation_matrix == 0, axis=0)]
#         print(correlation_matrix.shape)
        
        pinf = float('+inf')
        ninf = float('-inf')
        correlation_matrix[correlation_matrix == ninf] = 0.0
        correlation_matrix[correlation_matrix == pinf] = np.finfo(np.float16).max
        
        where_are_NaNs = np.isnan(correlation_matrix)
        correlation_matrix[where_are_NaNs] = 0
        
        np.savetxt('./scepisearch_query_results/fdr_epi.txt', final_fdr , delimiter=" ", fmt='%f')
        np.savetxt('./scepisearch_query_results/epi.txt', final_corr, fmt='%d', delimiter=" ")
        np.savetxt('./scepisearch_query_results/pval_epi.txt', pval_epi, delimiter=" ", fmt='%f')
#         np.savetxt('./scepisearch_query_results/autoencoder.txt', reduced_gene_enriched , delimiter=" ", fmt='%f')
        np.savetxt('./scepisearch_query_results/enrichment_scores.txt', gene_enriched , delimiter=" ", fmt='%f')
        

        net_size = int(np.ceil(np.sqrt(5*np.sqrt(epi.shape[1]))))
        net = sps.somNet(net_size, net_size, correlation_matrix, PBC=False)
        net.train(0.01, 1000)
        bmuList1 = net.project(correlation_matrix, show=False, printout=False)
        clust_infor = cluster_new(net, correlation_matrix, bmuList1, type= 'DBSCAN',cutoff=1, min_samples=1)
        clust_infor.iloc[:,3] = pd.Series(rand_jitter(np.array(clust_infor.iloc[:,3])))
        clust_infor.iloc[:,4] = pd.Series(rand_jitter(np.array(clust_infor.iloc[:,4])))
        clust_infor.to_csv('./scepisearch_query_results/tsne.txt',sep =" ", index = False , header = False)
        
        correlation_matrix = pd.DataFrame(correlation_matrix)
        #print(correlation_matrix.shape)
        if query_type == 1:
            meta_exp = pd.read_csv("./searchProject/meta_human/metadata_epi_new.csv", header=None,sep="@")
        else:
            meta_exp = pd.read_csv("./searchProject/meta_mouse/metadata_epi.csv", header=None,sep="@")
        #print(meta_exp.shape)
        meta_exp = meta_exp.iloc[:,0]
        meta_exp = np.array(meta_exp)
        meta_exp = meta_exp[ind_corrmat_nonzero]
        correlation_matrix.columns = meta_exp
        #print(correlation_matrix)
        ax = sns.clustermap(correlation_matrix)
        ax.savefig('./scepisearch_query_results/heatmap.png')

        correlation_matrix.to_csv('./scepisearch_query_results/correlation_matrix.txt', sep=" ", index=True,header=True)

        return gene_enriched, gene_name_exp_loc
    
def clusters_corr_null(null_model,exp_ref):
		nr, nc = exp_ref.shape
		xvec = ro.FloatVector(exp_ref.transpose().reshape((exp_ref.size)))
		xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

		nry, ncy = null_model.shape
		xvecy = ro.FloatVector(null_model.transpose().reshape((null_model.size)))
		yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

		res = ro.r.cor(xr,yr, method="spearman")
		res = np.array(res)
		res = np.transpose(res)
		return res


def cluster_epi_par(args):
        val,key,final_res,lock,query_type,gene_enriched = args
        if query_type==1:
            f0 = gzip.GzipFile('./searchProject/storage/scepisearch/human/epi_new/clusters_allgene_new_subclusters/clust_'+str(key)+'.npy.gz', "r")
            exp_ref = np.load(f0)
            null_model = np.load('./searchProject/storage/scepisearch/human/epi_new/null_model.npy')
            #corr_mat = np.load(f3)
            clust_infor = pd.read_csv("./searchProject/storage/scepisearch/human/epi_new/clusters_subclusters.txt",sep=" ")
        else:
            f0 = gzip.GzipFile('./searchProject/storage/scepisearch/mouse/epi/clusters_allgene/clust_'+str(key)+'.npy.gz', "r")
            exp_ref = np.load(f0)
            f3 = gzip.GzipFile('./searchProject/storage/scepisearch/mouse/epi/null_model/clust_'+str(key)+'.npy.gz', "r")
            corr_mat = np.load(f3)
            clust_infor = pd.read_csv("./searchProject/storage/scepisearch/mouse/epi/clusters.txt",sep=" ")

        pinf = float('+inf')
        ninf = float('-inf')
#         null_model[null_model == ninf] = 0.0
#         null_model[null_model == pinf] = np.finfo(np.float16).max
#         where_are_NaNs = np.isnan(null_model)
#         null_model[where_are_NaNs] = 0
        
        exp_ref[exp_ref == ninf] = 0.0
        exp_ref[exp_ref == pinf] = np.finfo(np.float16).max
        where_are_NaNs = np.isnan(exp_ref)
        exp_ref[where_are_NaNs] = 0
# print(np.argwhere(np.isnan(exp_ref)))
        
        if query_type == 1:
            corr_mat = clusters_corr_null(null_model,exp_ref)
            corr_mat[np.isnan(corr_mat)] = 0
        
#         print("key:",key)
#         print(clust_infor['cluster_no'] == str(key))
        if query_type == 1:
            ind_ref = list(clust_infor.loc[clust_infor['cluster_no'] == str(key), 'id'])
        else:
            ind_ref = list(clust_infor.loc[clust_infor['cluster_no'] == int(key), 'id'])
#         print(clust_infor.head())
#         print("ind_ref:",ind_ref)

        query = np.take(gene_enriched, list(val), axis=1)

        nr, nc = query.shape
        xvec = ro.FloatVector(query.transpose().reshape((query.size)))
        xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

        nry, ncy = exp_ref.shape
        xvecy = ro.FloatVector(exp_ref.transpose().reshape((exp_ref.size)))
        yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

        res = ro.r.cor(xr,yr, method="spearman")
        res = np.array(res)
        res = np.transpose(res)
        
#         res[np.isnan(res)] = 0

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
        pval = pval + 0.00001

        p_adjust_epi = stats.p_adjust(FloatVector(pval.ravel()), method = 'BH')
        p_adjust_epi = np.array(p_adjust_epi).reshape(pval.shape[0], pval.shape[1])
        
        print(sorted_10_idx.shape)
        print(len(ind_ref))
        result = np.array(ind_ref)[sorted_10_idx.ravel()]
        result = result.reshape(sorted_10_idx.shape[0],sorted_10_idx.shape[1])

        with lock:
            for i,j in enumerate(val):
                final_res[j,'corr'] = np.append(final_res[j,'corr'], val_top[:,i])
                final_res[j,'index'] = np.append(final_res[j,'index'], result[:,i])
                final_res[j,'pval'] = np.append(final_res[j,'pval'], pval[:,i])
                final_res[j,'adj_pval'] = np.append(final_res[j,'adj_pval'],p_adjust_epi[:,i])

def median_calc( x , exp_ref):
		ind_top = [x for x in x if x != -1]
		mat_top = np.take(exp_ref , ind_top , axis = 0)
		med_top = np.median(mat_top , axis=0)
		return med_top

def clusters_corr_epi(query_type,reduced_gene_enriched):
        if query_type==1:
            mean_array_reduced = np.load("./searchProject/storage/scepisearch/human/epi_new/mean_array.npy")
        else:
            mean_array_reduced = np.load("./searchProject/storage/scepisearch/mouse/epi/mean_array.npy")
            
        where_are_NaNs = np.isnan(mean_array_reduced)
        mean_array_reduced[where_are_NaNs] = 0
    
        nr, nc = reduced_gene_enriched.shape
        xvec = ro.FloatVector(reduced_gene_enriched.transpose().reshape((reduced_gene_enriched.size)))
        xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

        nry, ncy = mean_array_reduced.shape
        xvecy = ro.FloatVector(mean_array_reduced.transpose().reshape((mean_array_reduced.size)))
        yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

        res = ro.r.cor(xr,yr, method="spearman")
        res = np.array(res)
        res = np.transpose(res)
#         res[np.isnan(res)] = 0
        return res

def autoencoder_output(query_type,gene_enriched):
		gene_enriched = pd.DataFrame(np.transpose(gene_enriched))
		if(query_type==1):
			json_file = open('./searchProject/storage/scepisearch/human/epi/model_human.json', 'r')
			loaded_model_json = json_file.read()
			json_file.close()
			loaded_model = model_from_json(loaded_model_json)
			loaded_model.load_weights("./searchProject/storage/scepisearch/human/epi/model_human.h5")
		else:
			json_file = open('./searchProject/storage/scepisearch/mouse/epi/model_mouse.json', 'r')
			loaded_model_json = json_file.read()
			json_file.close()
			loaded_model = model_from_json(loaded_model_json)
			loaded_model.load_weights("./searchProject/storage/scepisearch/mouse/epi/model_mouse.h5")
		loaded_model.compile(loss = 'mean_squared_error', optimizer=sgd)
		encoder = Model(inputs=loaded_model.input, outputs=loaded_model.get_layer("dense_5").output)

		highest_non_inf = gene_enriched.max().loc[lambda v: v<np.Inf].max()
		gene_enriched_new = gene_enriched.replace(np.Inf, highest_non_inf + 100)
		gene_enriched_new = sklearn.preprocessing.scale(gene_enriched_new, axis = 0)
		reduced_gene_enriched = encoder.predict(gene_enriched_new)

		return np.transpose(np.array(reduced_gene_enriched))

def cluster_new(net, array, bmuList, type='DBSCAN', cutoff=1, min_samples=5):
		cluster_info = pd.DataFrame(index=range(array.shape[0]), columns=['id','color','cluster_no','xc','yc'])

		if type=='DBSCAN':
			cl = cluster.DBSCAN(eps=cutoff, min_samples=min_samples).fit(bmuList)

		j=0
		randCl = []
		for k in range(len(np.unique(cl.labels_))):
			randCl.append("#%06x" % np.random.randint(0, 0xFFFFFF))
		for i in range(len(cl.labels_)):
			cluster_info.iloc[j,0] = i
			cluster_info.iloc[j, 1] = randCl[cl.labels_[i]]
			cluster_info.iloc[j,2] = cl.labels_[i]
			cluster_info.iloc[j,3] = bmuList[i][0]
			cluster_info.iloc[j,4] = net.netHeight-bmuList[i][1]
			j = j+1

		return cluster_info
		
		
import pandas as pd
import numpy as np
import sklearn.preprocessing, multiprocessing, gzip
# from keras.models import model_from_json
# from keras.models import Model
from tensorflow.keras import optimizers
sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import numpy_indexed as npi
from collections import OrderedDict

def exp_matching(epi, gene_enriched, query_type, top_study, gene_name_exp_loc):
    #autoencoder output
#     reduced_gene_enriched = autoencoder_output(query_type, gene_enriched)

    #find relevant clusters
    res = clusters_corr_exp(query_type,gene_enriched)
#     print(res.shape)
    if query_type==1:
        index_value = np.genfromtxt("./human/clusters_celltypewise/subclusters/mean_array_subclusterindex.txt",dtype='str')
#         print(index_value)
    else:
        index_value = np.genfromtxt("./mouse_new_exp/clusters_celltype/mean_array_subclusterindex.txt",dtype='str')
    
    N = 10
    sorted_clust = np.argsort(res, axis=0)[res.shape[0]-N::,:]
#     print(sorted_clust)
    sorted_clust = index_value[sorted_clust.ravel()].reshape((sorted_clust.shape[0],sorted_clust.shape[1]))
    np.savetxt("./scepisearch_query_results/top_clusters.txt",sorted_clust,delimiter=" ",fmt='%s')

#     print(sorted_clust)
    un = np.unique(sorted_clust)
    id_clust = OrderedDict()
    for i,j in enumerate(un):
        id_clust[str(j)]=[]
        id_clust[str(j)][:] = np.where(sorted_clust == j)[1]
    
    workers= 10
    p = multiprocessing.Pool(processes=workers)
    keys, values= zip(*id_clust.items())

    manager = multiprocessing.Manager()
    final_res = manager.dict()
    lock = manager.Lock()

    for i in range(epi.shape[1]):
        final_res[i,'corr'] = np.array([])
        final_res[i,'index'] = np.array([])
        final_res[i,'pval'] = np.array([])
        final_res[i,'adj_pval'] = np.array([])

    p.map(cluster_exp_par, [(val,key,final_res,lock,query_type,gene_name_exp_loc) for val,key in zip(values,keys)])
    p.close()

    final_corr = np.zeros([epi.shape[1],top_study], dtype='int')
    pval_epi = np.zeros([epi.shape[1],top_study])
    final_fdr = np.zeros([epi.shape[1],top_study])
    
#     print(final_res)

    for i in range(epi.shape[1]):
        ind = np.argsort(final_res[i,'corr'])[::-1][:top_study]
        final_corr[i,:] = final_res[i,'index'][ind].astype(int)
#         final_med[i,:] = final_res[i,''][ind].astype(int)
        pval_epi[i,:] = final_res[i,'pval'][ind]
        final_fdr[i,:] = final_res[i,'adj_pval'][ind]

    np.savetxt('./scepisearch_query_results/fdr_exp.txt', final_fdr , delimiter=" ", fmt='%f')
    np.savetxt('./scepisearch_query_results/exp.txt', final_corr, fmt='%d', delimiter=" ")
    np.savetxt('./scepisearch_query_results/pval_exp.txt', pval_epi, delimiter=" ", fmt='%f')
#     np.savetxt('./scepisearch_query_results/median.txt',final_corr,delimiter=" ",fmt='%f')

def autoencoder_output(query_type, gene_enriched):
    gene_enriched = pd.DataFrame(np.transpose(gene_enriched))
    if(query_type==1):
        json_file = open('./searchProject/storage/scepisearch/human/exp/model_human.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights("./searchProject/storage/scepisearch/human/exp/model_human.h5")
    else:
        json_file = open('./searchProject/storage/scepisearch/mouse/exp/model_mouse.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights("./searchProject/storage/scepisearch/mouse/exp/model_mouse.h5")
    loaded_model.compile(loss = 'mean_squared_error', optimizer=sgd)
    encoder = Model(inputs=loaded_model.input, outputs=loaded_model.get_layer("dense_5").output)

    highest_non_inf = gene_enriched.max().loc[lambda v: v<np.Inf].max()
    gene_enriched_new = gene_enriched.replace(np.Inf, highest_non_inf + 100)
    gene_enriched_new = sklearn.preprocessing.scale(gene_enriched_new, axis = 0)
    reduced_gene_enriched = encoder.predict(gene_enriched_new)

    return np.transpose(np.array(reduced_gene_enriched))

def clusters_corr_exp(query_type,reduced_gene_enriched):
    if query_type==1:
        mean_array_reduced = np.load("./human/clusters_celltypewise/subclusters/mean_array.npy")
    else:
        mean_array_reduced = np.load("./mouse_new_exp/clusters_celltype/mean_array.npy")
        
#     where_are_NaNs = np.isnan(mean_array_reduced)
#     mean_array_reduced[where_are_NaNs] = 0
            
    for i in range(mean_array_reduced.shape[1]):
        mean_array_reduced[:,i]=mean_array_reduced[:,i]/np.mean(mean_array_reduced[:,i])
        
#     where_are_NaNs = np.isnan(mean_array_reduced)
#     mean_array_reduced[where_are_NaNs] = 0
    
    where_are_NaNs = np.isnan(reduced_gene_enriched)
    reduced_gene_enriched[where_are_NaNs] = 0

    nr, nc = reduced_gene_enriched.shape
    xvec = ro.FloatVector(reduced_gene_enriched.transpose().reshape((reduced_gene_enriched.size)))
    xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

    nry, ncy = mean_array_reduced.shape
    xvecy = ro.FloatVector(mean_array_reduced.transpose().reshape((mean_array_reduced.size)))
    yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

    res = ro.r.cor(xr,yr, method="spearman")
    res = np.array(res)
    res = np.transpose(res)
    return res

def median_calc(x , exp_ref):
		ind_top = [x for x in x if x != -1]
		mat_top = np.take(exp_ref , ind_top , axis = 0)
		med_top = np.median(mat_top , axis=0)
		return med_top
    
def median_calc_null(x , corr_mat, exp_ref):
#     ind_top = [x for x in x if x != -1]
    mat_top = np.take(exp_ref , x , axis = 0)
#     print(mat_top.shape)
    med_top = np.median(mat_top , axis=0)
#     print(med_top.shape)
    corr_mat[x.name,:] = med_top

def cluster_exp_par(args):
    val,key,final_res,lock,query_type,gene_name_exp_loc = args
    if(query_type == 1):
        f0 = gzip.GzipFile('./human/clusters_celltypewise/subclusters/clust_'+str(key)+'.npy.gz', "r")
        exp_ref = np.load(f0)
#         f3 = gzip.GzipFile('./searchProject/storage/scepisearch/human/exp/null_model/clust_'+str(key)+'.npy.gz', "r")
#         corr_mat = np.load(f3)
        sorted_col = np.load("./null_idx_enrichment.npy")
        clust_infor = pd.read_csv("./human/clusters_celltypewise/subclusters/clusters_final.txt",sep=" ",dtype='str')
#         print(sorted_col.shape)
        unanno = np.loadtxt("./searchProject/meta_human/exp_unknown_human.txt", delimiter=",")
    else:
        f0 = gzip.GzipFile('./mouse_new_exp/clusters_celltype/subclusters_twentymillion/clust_'+str(key)+'.npy.gz', "r")
        exp_ref = np.load(f0)
#         f3 = gzip.GzipFile('./searchProject/storage/scepisearch/mouse/exp/null_model/clust_'+str(key)+'.npy.gz', "r")
#         corr_mat = np.load(f3)
        sorted_col = np.loadtxt("./gene_name_exp_loc_mouse.txt",delimiter=",")
        sorted_col=sorted_col.astype(int)
        clust_infor = pd.read_csv("./mouse_new_exp/clusters_celltype/clusters_final.txt",sep=" ",dtype='str')
        metadata = pd.read_csv("./searchProject/meta_mouse/metadata_exp.csv", header=None,sep="@")
        metadata = metadata.iloc[:,2]
        nan_rows = metadata[metadata.isnull()]
        unanno = np.array(nan_rows.index)
    ind_ref = np.array(clust_infor.loc[clust_infor['cluster_no'] == str(key), 'id'])
    exp_ref = exp_ref[:,~np.isin(ind_ref.astype(int),unanno.astype(int))]
    ind_ref = ind_ref[~np.isin(ind_ref.astype(int),unanno.astype(int))]
#     print(len(ind_ref))
#     print(exp_ref.shape)
    
    if exp_ref.shape[1]==0:
        return
    
    for i in range(exp_ref.shape[1]):
        exp_ref[:,i]=exp_ref[:,i]/np.mean(exp_ref[:,i])

    where_are_NaNs = np.isnan(exp_ref)
    exp_ref[where_are_NaNs] = 0

    query = np.take(gene_name_exp_loc, list(val), axis=1)
    res = query.apply(lambda x : median_calc(x , exp_ref) ,axis = 0)
    res = np.array(res)
    
    corr_mat = np.zeros([1000,exp_ref.shape[1]])
    pd.DataFrame(sorted_col).apply(lambda x : median_calc_null(x , corr_mat ,exp_ref) ,axis = 0)
    corr_mat = np.array(corr_mat)
#     print(corr_mat)

    #DO RANKING BASED ON corr_mat
    score_ranking = np.argsort(corr_mat, axis=1)[:,corr_mat.shape[1]-50::]
    
    top_ref = np.unique(score_ranking)
    
    res = res[top_ref,:]
    
    pval = []
    sorted_10_idx = np.argsort(res, axis=0)[res.shape[0]-10::,:]
    sorted_raveled = sorted_10_idx.ravel()
    col_idx = np.arange(res.shape[1])
    val_top = res[sorted_10_idx, col_idx]
    val_top_raveled = val_top.ravel()
    for a, b in zip(sorted_raveled, val_top_raveled):
        pval.append(sum(corr_mat[:,a] > b))
    pval = np.reshape(pval , (sorted_10_idx.shape[0] , sorted_10_idx.shape[1]))

    pval = pval/float(1000)
    pval = pval + 0.00001
#     print(pval)

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
			
			
import linecache, os
import numpy_indexed as npi
from multiprocessing import Process
import itertools
def pdf_gen(query_file,query_type):
		if query_type == 1:
			metadata_epi = './searchProject/meta_human/metadata_epi.csv'
			metadata_exp = './searchProject/meta_human/metadata_exp_with_tissue.csv'
		elif query_type == 2:
			metadata_epi = './searchProject/meta_mouse/metadata_epi.csv'
			metadata_exp = './searchProject/meta_mouse/metadata_exp.csv'
		else:
			metadata_exp = './searchProject/meta_mouse/metadata_exp.csv'
		if query_type == 1 or query_type == 2:
			output_file_exp = open('./scepisearch_query_results/exp.txt')
			output_file_epi = open('./scepisearch_query_results/epi.txt')
			file_pval_exp = open('./scepisearch_query_results/pval_exp.txt')
			file_pval_epi = open('./scepisearch_query_results/pval_epi.txt')
			file_fdr_epi = open('./scepisearch_query_results/fdr_epi.txt')
			file_fdr_exp = open('./scepisearch_query_results/fdr_exp.txt')
			count = 0
			for line1,line2,line3,line4,line5,line6 in zip(output_file_exp,file_pval_exp,file_fdr_exp,output_file_epi,file_pval_epi,file_fdr_epi):
				matrix1 = list()
				matrix2 = list()
				line1 = line1.split(' ')
				line2 = line2.split(' ')
				line3 = line3.split(' ')
				line4 = line4.split(' ')
				line5 = line5.split(' ')
				line6 = line6.split(' ')
				line1[-1] = line1[-1].strip()
				line2[-1] = line2[-1].strip()
				line3[-1] = line3[-1].strip()
				line4[-1] = line4[-1].strip()
				line5[-1] = line5[-1].strip()
				line6[-1] = line6[-1].strip()
				temp_filename_exp = "./scepisearch_query_results/query_exp_"+str(count)+".txt"
				temp_filename_epi = "./scepisearch_query_results/query_epi_"+str(count)+".txt"
				with open(temp_filename_exp, "w") as f1 , open(temp_filename_epi, "w") as f2:
					for i,j,k,l,m,n in zip(line1, line2, line3, line4, line5, line6):
						line = linecache.getline(metadata_exp, int(i)+1)
						line = line.strip()
						line = line+"@"+j+"@"+k
						f1.write(line+'\n')
						line_epi = linecache.getline(metadata_epi, int(l)+1)
						line_epi = line_epi.strip()
						line_epi = line_epi+"@"+m+"@"+n
						f2.write(line_epi+'\n')
				count = count + 1
			filename = generate_pdf(count,query_type)
			for i in range(int(count)):
				os.remove('./scepisearch_query_results/query_exp_' + str(i) + '.txt')
				os.remove('./scepisearch_query_results/query_epi_' + str(i) + '.txt')
		else:
			output_file_exp = open('./scepisearch_query_results/exp.txt')
			file_pval_exp = open('./scepisearch_query_results/pval_exp.txt')
			file_fdr_exp = open('./scepisearch_query_results/fdr_exp.txt')
			count = 0
			for line1,line2,line3 in zip(output_file_exp,file_pval_exp,file_fdr_exp):
				matrix1 = list()
				matrix2 = list()
				line1 = line1.split(' ')
				line2 = line2.split(' ')
				line3 = line3.split(' ')
				line1[-1] = line1[-1].strip()
				line2[-1] = line2[-1].strip()
				line3[-1] = line3[-1].strip()
				temp_filename_exp = "./scepisearch_query_results/query_exp_"+str(count)+".txt"
				with open(temp_filename_exp, "w") as f1:
					for i,j,k in zip(line1, line2, line3):
						line = linecache.getline(metadata_exp, int(i)+1)
						line = line.strip()
						line = line+"@"+j+"@"+k
						f1.write(line+'\n')
				count = count + 1
			filename = generate_pdf(count,query_type)
			for i in range(int(count)):
				os.remove('./scepisearch_query_results/query_exp_' + str(i) + '.txt')
				
				
## pdf report_generator
# from reportlab.pdfgen import canvas
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import Table, TableStyle
from reportlab.lib.units import cm, inch
from reportlab.lib import colors
from reportlab.lib.units import inch
import matplotlib.pyplot as plt
import numpy as np
import csv
import time
from wordcloud import WordCloud, STOPWORDS 
import matplotlib.pyplot as plt 

def generate_pdf(query_cells,query_type):
	filename = './scepisearch_query_results/results.pdf'
	doc = SimpleDocTemplate(filename,
                        rightMargin=72,leftMargin=72,
                        topMargin=72,bottomMargin=18)

	Story=[]

	logo = './searchProject/logo.jpg'
	im = Image(logo, 6*inch, 1.5*inch)
	Story.append(im)

	full_name = "ScEpiSearch Report"
	address_parts = ["IIIT Delhi", "New Delhi"]
	formatted_time = time.ctime()

	styles=getSampleStyleSheet()
	styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

	ptext = '<font size=12>%s</font>' % formatted_time
	Story.append(Paragraph(ptext, styles["Normal"]))

	Story.append(Spacer(1, 12))

	# Create return address_parts
	ptext = '<font size=12>%s</font>' % full_name
	Story.append(Paragraph(ptext, styles["Normal"]))

	Story.append(Spacer(1, 12))

	for part in address_parts:
	    ptext = '<font size=12>%s</font>' % part.strip()
	    Story.append(Paragraph(ptext, styles["Normal"]))

	Story.append(Spacer(1, 12))
	Story.append(Spacer(1, 12))

	all_cells = [(0, 0), (-1, -1)]
	header = [(0, 0), (-1, 0)]
	column0 = [(0, 0), (0, -1)]
	column1 = [(1, 0), (1, -1)]
	column2 = [(2, 0), (2, -1)]
	column3 = [(3, 0), (3, -1)]
	column4 = [(4, 0), (4, -1)]
	table_style = TableStyle([('VALIGN', all_cells[0], all_cells[1], 'TOP'),
		('FONTSIZE', all_cells[0], all_cells[1], 6),
	    ('LINEBELOW', header[0], header[1], 1, colors.black),
	    ('ALIGN', column0[0], column0[1], 'LEFT'),
	    ('ALIGN', column1[0], column1[1], 'LEFT'),
	    ('ALIGN', column2[0], column2[1], 'LEFT'),
	     ('ALIGN', column3[0], column3[1], 'LEFT'),
    	('ALIGN', column4[0], column4[1], 'LEFT'),
	])

	# PDF Table - Column Widths
	# colWidths = [2 * cm,2* cm, 10 * cm,2.5 * cm,2.5 * cm]

	ptext = '<font size=8>Matching in RNA-Seq Data : </font>'
	Story.append(Paragraph(ptext, styles["Normal"]))
	cells_all = list()
	
	if query_type == 1 or query_type ==2:
		for i in range(int(query_cells)):
			ptext = '<font size=8>Query Cell Number : %s</font>' % str(i)
			Story.append(Paragraph(ptext, styles["Normal"]))
			data_exp = list()
			data_exp.append(['CellId', 'StudyId', 'Phenotype', 'P value', 'Adjusted P value'])
			temp_filename = './scepisearch_query_results/query_exp_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
			for i in range(len(reader)):
				if query_type == 1:
					reader[i] = reader[i][1:]
				cells_all.append(reader[i][2])	
			reader = data_exp + reader

			for index, row in enumerate(reader):
				for col, val in enumerate(row):
					# if col == 4 or col==0 or index == 0:
						# reader[index][col] = val.strip("'[]")
					# else:
					ptext = '<font size=6>%s</font>' % val
					reader[index][col] = Paragraph(ptext, styles['Normal'])

			t = Table(reader)
			t.setStyle(table_style)
			Story.append(t)
			Story.append(Spacer(1, 12))
			Story.append(Spacer(1, 12))

		Story.append(Spacer(1, 12))
		Story.append(Spacer(1, 12))

		ptext = '<font size=8>Matching in Single Cell Epigenome Data: </font>'
		Story.append(Paragraph(ptext, styles["Normal"]))

		for i in range(int(query_cells)):
			ptext = '<font size=8>Query Cell Number : %s</font>' % str(i)
			Story.append(Paragraph(ptext, styles["Normal"]))
			data_exp = list()
			data_exp.append(['CellId', 'StudyId', 'Phenotype', 'P value', 'Adjusted P value'])
			temp_filename = './scepisearch_query_results/query_epi_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
			for i in range(1,len(reader)):
				cells_all.append(reader[i][2])
			reader = data_exp + reader

			for index, row in enumerate(reader):
				for col, val in enumerate(row):
					# if col != 4 or index == 0:
						# reader[index][col] = val.strip("'[]()")
					# else:
					ptext = '<font size=6>%s</font>' % val
					reader[index][col] = Paragraph(ptext, styles['Normal'])

			t = Table(reader)
			t.setStyle(table_style)
			Story.append(t)
			Story.append(Spacer(1, 12))
			Story.append(Spacer(1, 12))
	else:
		for i in range(int(query_cells)):
			ptext = '<font size=8>Query Cell Number : %s</font>' % str(i)
			Story.append(Paragraph(ptext, styles["Normal"]))
			data_exp = list()
			data_exp.append(['CellId', 'StudyId', 'Phenotype', 'P value', 'Adjusted P value'])
			temp_filename = './scepisearch_query_results/query_exp_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
			for i in range(len(reader)):
				cells_all.append(reader[i][2])	
			reader = data_exp + reader

			for index, row in enumerate(reader):
				for col, val in enumerate(row):
					# if col == 4 or col==0 or index == 0:
						# reader[index][col] = val.strip("'[]")
					# else:
					ptext = '<font size=6>%s</font>' % val
					reader[index][col] = Paragraph(ptext, styles['Normal'])

			t = Table(reader)
			t.setStyle(table_style)
			Story.append(t)
			Story.append(Spacer(1, 12))
			Story.append(Spacer(1, 12))

		Story.append(Spacer(1, 12))
		Story.append(Spacer(1, 12))

		ptext = '<font size=8>Matching in Single Cell Epigenome Data: </font>'
		Story.append(Paragraph(ptext, styles["Normal"]))
		
	doc.build(Story)
	comment_words = '' 
	stopwords = set(STOPWORDS) 
	cells_all = np.array(cells_all)
	#print(cells_all)
	cells_all = cells_all.astype(str)
	comment_words += " ".join(cells_all)+" "
	wordcloud = WordCloud(width = 800, height = 800, background_color ='white', stopwords = stopwords, min_font_size = 10).generate(comment_words) 
  
	# plot the WordCloud image                        
	plt.figure(figsize = (8, 8), facecolor = None) 
	plt.imshow(wordcloud) 
	plt.axis("off") 
	plt.tight_layout(pad = 0) 
	#plt.show()
	plt.savefig('./scepisearch_query_results/wordcloud.png')
	 
	return filename


#!/usr/bin/python
import numpy as np , pandas as pd
from multiprocessing import Process
import numpy_indexed as npi

def get_digit(x):
    return(int(re.search(r'\d+', x).group()))

def process_query(chr_file, epi_path, top_study, query_type, acc_fast,active_poised,imputation):
    if(query_type==1):
        sc_gene_path = './searchProject/storage/scepisearch/human/genes_21159.txt'
    else:
        sc_gene_path = './searchProject/storage/scepisearch/mouse/gene_mouse.csv'

    gene_list = list()
    with open(sc_gene_path, 'r') as fl:
        for l in fl.readlines():
            gene_list.append(l.rstrip('\n'))

    nearest_gene,epi = nearest_gene_accurate(query_type,chr_file,acc_fast,epi_path)
#     print(nearesht_gene)
#     nearest_gene=pd.read_csv("./foreground.csv",header=None,sep=",")
    # epi=np.loadtxt(epi_path,delimiter=",")
    if not nearest_gene.empty:
        acc_score_path = './acc_score.csv'
        acc_score = np.loadtxt(acc_score_path)
        # epi = np.loadtxt(epi_path, delimiter=",")
        acc_score[acc_score == 0] = 1

        epi = np.array(epi)/(acc_score[:, np.newaxis])
#         np.savetxt("./scepisearch_query_results/k562_100cells/normalized_count.txt",epi,delimiter=",")

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
        
#         gene_enriched = np.loadtxt("./scepisearch_query_results/enrichment_scores.txt",delimiter=" ")
        gene_enriched = gene_enrichment_calc(epi,gene_list,nearest_gene,query_type)
        if int(active_poised)==1:
            N = 500
        else:
            N = 2000
        sorted_col_idx_exp = np.argsort(epi, axis=0)[epi.shape[0]-N::,:]
        top_epi_gene_exp = pd.DataFrame(sorted_col_idx_exp).apply(lambda x: np.take(epi_gene, indices = x) , axis = 0)
        gene_name_exp_loc_exp = top_epi_gene_exp.apply(lambda x : npi.indices(gene_list ,x, missing = -1)  , axis = 0)
        
        sorted_col_idx_epi = np.argsort(gene_enriched, axis=0)[gene_enriched.shape[0]-N::,:]
        top_epi_gene_epi = pd.DataFrame(sorted_col_idx_epi).apply(lambda x: np.take(gene_list, indices = x) , axis = 0)
        gene_name_exp_loc_epi = top_epi_gene_epi.apply(lambda x : npi.indices(gene_list ,x, missing = -1)  , axis = 0)


        p1 = Process(target=epi_matching, args=(epi,query_type,gene_list,nearest_gene,top_study,epi_gene,gene_enriched,gene_name_exp_loc_epi))
        p2 = Process(target=exp_matching, args=(epi,gene_enriched,query_type,top_study,gene_name_exp_loc_exp))
        p1.start()
        p2.start()
        p1.join()
        p2.join()

        print('cell type search done')
    else:
        print("error")
		
		
chr_file='./human/human_epigenome/5_queries/h1esc/h1esc.bed'
query_file='./human/human_epigenome/5_queries/h1esc/h1esc_human.txt'
top_study=5
query_type=1
acc_fast=2
active_poised=2
imputation=1

process_query(chr_file,query_file,top_study,query_type,acc_fast,active_poised,imputation)
pdf_gen(query_file,query_type=3)
