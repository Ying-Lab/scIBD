# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 14:49:41 2021

@author: User
"""
######final version

import sys
import numpy as np
import pandas as pd
import anndata
import scipy 
import math
from scipy import spatial
import os
import random
import time
import tracemalloc
import multiprocessing as mp
import psutil
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import scipy.io as scio
import umap
import skbio
from mpl_toolkits.mplot3d import axes3d
from scipy import stats
import seaborn as sns
import scanpy as sc
###plt setting
import warnings
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve, auc
from sklearn import preprocessing
from sklearn.metrics import calinski_harabasz_score
from sklearn.feature_extraction.text import TfidfTransformer
import matplotlib
warnings.filterwarnings("ignore")
random.seed(1111)
np.random.seed(1111)


# def TF_IDF_Pro(count_mat):
#     nfreqs = (1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))).A
#     tfidf_mat = nfreqs * (np.tile(np.log(1 + 1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1]))).A
#     if np.sum(np.isnan(tfidf_mat) == True) != 0:
#         #print("The results of TF-IDF contain nan, replacing it with zero!")
#         tfidf_mat[np.isnan(tfidf_mat)] = 0 
#     if np.sum(np.isinf(tfidf_mat) == True) != 0:
#         print("The results of TF-IDF contain inf, replacing it with mean value!")
#         tfidf_mat = fill_ndarray(tfidf_mat)
# #     print('tf-idf done')
#     return scipy.sparse.csr_matrix(tfidf_mat.T)

# Perform TF-IDF (count_mat: peak*cell)
def TF_IDF_Pro(count_mat): 
    
    if not scipy.sparse.issparse(count_mat):
        count_mat = scipy.sparse.coo_matrix(count_mat)
        
    nfreqs = count_mat.multiply(1.0 / count_mat.sum(axis=0))
    tfidf_mat = nfreqs.multiply(np.log(1 + 1.0 * count_mat.shape[1] / count_mat.sum(axis=1)).reshape(-1,1)).tocoo()
    
    return tfidf_mat

# def TF_IDF_Pro(count_mat): ##cell*peak
#     model = TfidfTransformer(smooth_idf=False, norm="l2")
#     model = model.fit(np.transpose(count_mat))
#     model.idf_ -= 1
#     tf_idf = np.transpose(model.transform(np.transpose(count_mat)))
#     return scipy.sparse.csr_matrix(tf_idf)

def CalJac(A,B=None):##scipy mat
    A = A.astype(bool).astype(int)
    if B is None:
        B = A
    else:
        B = B.astype(bool).astype(int)
    intersect = A.dot(B.T)
    a_sum = A.sum(axis=1).A1
    b_sum = B.sum(axis=1).A1
    xx,yy = np.meshgrid(a_sum, b_sum)
    union = ((xx+yy).T - intersect)
    jac_distance = (1-intersect/union).A
    return jac_distance

def F1(score, label):
    assert len(score) == len(label)
    P  = len(label[label=='DOUB'])
    PP = len(score[np.all([score!= 0,score!= 0.05],axis = 0)]) #TP+FP
    TP = len(label[np.all([score!= 0,score!= 0.05],axis = 0)][label[np.all([score!= 0,score!= 0.05],axis = 0)]=='DOUB'])
    precision = TP/PP
    recall = TP/P
    F1 = 2* precision*recall/(precision+recall)
    return precision,recall,F1

def random_color():
    colors1 = '0123456789ABCDEF'
    num = "#"
    for i in range(6):
        num += random.choice(colors1)
    return num

####testing umap then dbscan

## get clusters based on PCA+UMAP+DBSCAN 
def GetCluster(mat):
    PCA_umap = GetUMAP(mat)
    PCA_umap = StandardScaler().fit_transform(PCA_umap) 
    
    distance = metrics.pairwise_distances(PCA_umap, Y=None, metric='euclidean', n_jobs=64)
    Parameters ={'eps':0.1,'min_samples':int(mat.shape[0]*0.005)}
    db = DBSCAN(eps = Parameters['eps'], min_samples = Parameters['min_samples'],metric = 'precomputed').fit(distance) 
    cl = db.labels_
    return cl,distance

def random_weight(i,weight_data):
    random.seed(i)
    total = sum(weight_data.values())    # 权重求和
    ra = random.uniform(0, total)   # 在0与权重和之前获取一个随机数 
    curr_sum = 0
    ret = None
    keys = weight_data.keys()        # 使用Python3.x中的keys
    for k in keys:
        curr_sum += weight_data[k]             # 在遍历中，累加当前权重值
        if ra <= curr_sum:          # 当随机数<=当前权重和时，返回权重key
            ret = k
            break
    return ret

####variable peaks based on cluster 
def findpeaks(mat,cl,rate):
#     cl,cd = GetCluster(jac)
    cl_uni = list(set(cl))
    cl_score = np.array([]).reshape(0,mat.shape[1])
    for item in cl_uni:
        cl_score = np.vstack((cl_score,np.mean(mat[cl == item],axis = 0))).A
    var_peak = np.var(np.log2(cl_score+1),axis=0)
    idxbool = var_peak > np.quantile(var_peak,rate,interpolation='higher')
    return idxbool

##import scipy
def GetUMAP(mat):
    var_names = pd.DataFrame(np.array(range(mat.shape[1])), columns=['peak_names'])
    col_names = pd.DataFrame(np.array(range(mat.shape[0])), columns=['cell_names'])
    adata = sc.AnnData(TF_IDF_Pro(mat.T).T,obs=col_names, var=var_names,)
    ##this step used the former TF-IDF matrix
    adata.obs_names = col_names.cell_names
    adata.var_names = var_names.peak_names
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=30,n_pcs=30)
    sc.tl.umap(adata,min_dist=0.4)
    return adata.obsm['X_umap']

##import scipy
def GetSparseUnion(A, B = None):
    if B is None:
        B = A
    intersect = A.dot(B.T)
    A_sum = A.sum(axis=1).A1
    B_sum = B.sum(axis=1).A1
    xx,yy = np.meshgrid(A_sum,B_sum)
    union = ((xx+yy).T - intersect)
    return union

##import scipy
def CreateDoublets_Cluster_weighted(simulate_size, mat, distancemat, clusterlab):
#     print('simulating start',simulate_size)
    np.random.seed(1111)
    assert simulate_size != 0
#     mat = scipy.sparse.csr_matrix(mat.astype(np.uint8))
    clusterlab = clusterlab.astype(str)     
    cluster_uniq = list(set(clusterlab))
#     print(len(cluster_uniq), 'clusters')

    ####设置单个类大小权重dict
    cluster_weight = {}
    ####设置两类大小权重dict
    cluster_bi_weight = {}
    ###每个cluster的权重
    for i in range(len(cluster_uniq)):
        cluster_weight[str(cluster_uniq[i])] = clusterlab[clusterlab == cluster_uniq[i]].shape[0]
    ###归一化
    cw_st = MinMaxScaler(feature_range=(0,1)).fit_transform(np.array(list(cluster_weight.values())).reshape(-1,1)).squeeze()
    
    for i in range(len(cluster_uniq)):
        for j in range(i+1,len(cluster_uniq)):
            cluster_bi_weight[str(cluster_uniq[i])+'_'+str(cluster_uniq[j])] = 0
    
    idx = 0 
    for item in cluster_weight:
        cluster_weight[item] = cw_st[idx]
        idx += 1 
        
    for item in cluster_bi_weight:
        cluster_bi_weight[item] = cluster_weight[item.split('_')[0]] * cluster_weight[item.split('_')[1]]

#     single1set = []
#     single2set = []
#     simulatelabel = []
    L1_size = int(simulate_size/3*2)
    ####若聚类只有1，全部随机
    
    if len(cluster_uniq) == 1:
        L1_size = 0 
    for i in range(L1_size):

        ramtype = random_weight(i,cluster_bi_weight).split('_')
        
        idx1 = np.random.choice(len(clusterlab[clusterlab == ramtype[0]]),1)
        idx2 = np.random.choice(len(clusterlab[clusterlab == ramtype[1]]),1)
        try:
            single1set = scipy.sparse.vstack((single1set,mat[clusterlab == ramtype[0]][idx1]))
            single2set = scipy.sparse.vstack((single2set,mat[clusterlab == ramtype[1]][idx2]))
        except:
            single1set = mat[clusterlab == ramtype[0]][idx1]
            single2set = mat[clusterlab == ramtype[1]][idx2]
#         single1set.append(mat[clusterlab == ramtype[0]][idx1].A)
#         single2set.append(mat[clusterlab == ramtype[1]][idx2].A)。
#         doublet = scipy.sparse.csr_matrix(mat[clusterlab == ramtype[0]][idx1].A | mat[clusterlab == ramtype[1]][idx2].A)

    
    for i in range(simulate_size - L1_size):

        idxset = np.random.choice(mat.shape[0],2)
        single1set = scipy.sparse.vstack((single1set,mat[idxset[0]]))
        single2set = scipy.sparse.vstack((single2set,mat[idxset[1]]))
#         single1set.append(mat[idxset[0]].A)
#         single2set.append(mat[idxset[1]].A)
#         doublet = scipy.sparse.csr_matrix(mat[clusterlab == ramtype[0]][idx1].A | mat[clusterlab == ramtype[1]][idx2].A)


#     sim_set = GetSparseUnion(single1set,single2set)
    sim_set = (single1set + single2set).astype(np.uint8)

#     sim_set = scipy.sparse.csr_matrix(single1set|single2set)
#     simulatelabel = np.array(simulatelabel)
    print(simulate_size,'Doublets Construction Done!')
    return sim_set




def Splitmat(matsize,ncore):
    idx_set =[]
#     mat_size = mat.shape[0]
    size_scale = round(matsize/ncore)
    for i in range(ncore-1):
        idx_set.append(list(range(i*size_scale,(i+1)*size_scale)))
    idx_set.append(list(range((ncore-1)*size_scale, matsize)))
    return idx_set


##import csr
def ReshapeJac(mat_conc,label_conc,jac_raw,core):
#     print('reshape jaccard')
    A = mat_conc[label_conc != 1]  
    B = mat_conc[label_conc == 1]
    ###判断是不是simulated doublets
#     jac_extra = metrics.pairwise_distances(B.A, Y=mat_conc.A, metric='jaccard', n_jobs=core)
    jac_extra = CalJac(B,mat_conc)
    jac_new = np.vstack((jac_raw, jac_extra[:,0:A.shape[0]]))
    jac_new = np.hstack((jac_new,jac_extra.T))

#     for i in range(len(jac_new)):
#         assert jac_new[i][i] == 0 
    return jac_new

def PCoA(mat_jaccard,npc,showfig=False):

#     mat_jaccard = Jaccard_Index(mat)
    OrdinationResults = skbio.stats.ordination.pcoa(mat_jaccard, method='eigh', number_of_dimensions=0, inplace=False)
    return np.array(OrdinationResults.samples)[:,0:npc]

###validate when labels are known
def concat_set(iteration, jac,mat, sim_rate, label_refer):
#     print('concat start')
#     label_refer {0:singlets;
#              0.05:doublets reject
#              0.1-0.9:doublets predicted
#              1: simulated doublets}
    ####先对未被分类的细胞做聚类，基于原始的jaccard；label_refer为0或0.05 
    mat_unlabel = mat[np.any([label_refer == 0, label_refer == 0.05],axis = 0)]
    jac_unlabel = jac[np.any([label_refer == 0, label_refer == 0.05],axis = 0)][:,np.any([label_refer == 0, label_refer == 0.05],axis = 0)]
    cl_refer = label_refer[np.any([label_refer == 0, label_refer == 0.05],axis = 0)]
    cl, cd = GetCluster(mat_unlabel)
    #######在这一步决定simsize，为未分类样本数的10% - 40% 
    sim_size = int(sim_rate * mat_unlabel.shape[0])
    time_start = time.time()
    simDOU = CreateDoublets_Cluster_weighted(sim_size, mat_unlabel,cd,cl)
    mat_concat = scipy.sparse.vstack((mat,simDOU)).astype(np.uint8)
    label_refer[label_refer==0.05] = 0
    label_concat = np.hstack((label_refer,np.array([1]*sim_size))).astype(np.float64)

    return mat_concat,label_concat

##the size of sim_label corresponds to the sim_size

def get_knn_graph(mat, k, n_tree):
    from annoy import AnnoyIndex
    npc = mat.shape[1]
    ##cells
    ncell = mat.shape[0]
    ##peaks
#     annoy_index.unload() 
    ###存储npc 维的向量，metric ："euclidean"
    annoy_index = AnnoyIndex(npc, metric='euclidean')

    ###为索引i添加特征向量（peak）
    for i in range(ncell):
        annoy_index.add_item(i, list(mat[i,:]))
   
    ####建立 n_trees 的森林。查询时，树越多，精度越高。
    annoy_index.build(n_tree) 

    knn = []
    distance = []
    for iCell in range(ncell): 
        neighs = annoy_index.get_nns_by_item(iCell, k + 1)[1:]
        knn.append(neighs)
        nei_dis = []
        for neigh in neighs:
            nei_dis.append(annoy_index.get_distance(iCell,neigh))
        distance.append(nei_dis)
    knn = np.array(knn)
    distance = np.array(distance)
    return knn, distance 

def CalKNNScore(knn,distance,min_dis,max_dis,label):
    score_knn = np.zeros(len(knn))
    ###cell size
    for i in range(len(knn)):
#         neigh_labels = labelmat_concat[knn[i]]
        neigh_refer = label[knn[i]]
        score_knn[i] = np.dot((max_dis-distance[i])/(max_dis-min_dis) , neigh_refer)
    return score_knn

def Calsigmathresh(score1,score2):
    kernel1 = stats.gaussian_kde(score1)
    kernel2 = stats.gaussian_kde(score2)
    X = np.mgrid[np.min(score2):np.max(score2):0.1]
    idx = np.argwhere(np.diff(np.sign(kernel1(X)-kernel2(X)))).flatten()
#     scoremean = np.mean(score1)
    scorestd = np.std(score1)
    sigma = 2 * scorestd 
    thresh = X[idx] + sigma
    if len(thresh) != 0:
        return X[idx],kernel1(X)[idx],thresh[0]
    else:
        return 0, 0, np.mean(score2) - sigma

# def CalKL(a,b,scalesize):
#     minvalue = min(a.min(),b.min())
#     maxvalue = max(a.max(),b.max())
#     scalevalue = maxvalue-minvalue
#     proba, probb = np.zeros(scalesize),np.zeros(scalesize)
#     for i in range(scalesize):
#         proba[i] = (len(a[np.all([a>=minvalue+i*(scalevalue/scalesize), a<=minvalue+(i+1)*(scalevalue/scalesize)],axis=0)]))/len(a)+ 0.00000001
#         probb[i] = (len(b[np.all([b>=minvalue+i*(scalevalue/scalesize), b<=minvalue+(i+1)*(scalevalue/scalesize)],axis=0)]))/len(b)+ 0.00000001
#     KLvalue = scipy.stats.entropy(proba,probb)
#     return KLvalue

def Ramsplit(lenset,core):
    random.seed(1111)
    idxset = []
    setall = set(range(lenset))
    eachset = round(lenset/core)
    for i in range(core-1):
        tmpset = random.sample(list(setall),eachset)
        idxset.append(tmpset)
        setall = setall-set(tmpset)
    idxset.append(list(setall))
    return idxset




# def PlotDensity(score_raw,score_sim,score_detected,labelmat_det):
#     # score_raw, score_sim, score_detected = CallDoublet(mat, 800, jaccard_raw_mat, cores, 5, 40, 30,label_refer_mat)
# # import seaborn as sns
#     fig = plt.figure(1,figsize=(20,10), dpi= 100)

#     plt.title('Density')

#     plt.xlabel('Score')
#     plt.ylabel('Probability')
    
#     label_uniq = ['sim','detected','raw','crosspoint','sigma','quantitle']

#     score_detected = score_detected.values.squeeze()
#     sns.set() 

# #     sns.kdeplot(score_sim, shade=True, color="#F9F871",alpha=.5)
#     #background-image: linear-gradient(to right bottom, #5628e4, #fa0094, #ff2736, #faa200, #a8eb12);
#     sns.distplot(score_sim, color="yellow")
#     sns.distplot(score_detected,color="red")
# #     sns.distplot(score_detected[labelmat_det!='DOUB'], color="#F3C5FF")
#     sns.distplot(score_raw.values.reshape(-1),color="grey")

#     thresh = Cal3sigmathresh(score_raw.values.reshape(-1),score_sim)
#     plt.plot(thresh[0],thresh[1], 'ro')###交点
#     plt.axvline(thresh[2],ls=":",lw=2,c="blue") ##3sigma
#     plt.axvline(np.quantile(score_sim,0.6,interpolation='higher'),ls="--",lw=2,c="black")##60%分位

#     ax = fig.gca()
#     ax.patch.set_facecolor("#E6F4F1") 
#     ax.patch.set_alpha(0.2) 
#     ax.grid(color='r',
#             linestyle='--',
#             linewidth=1,
#             alpha=0.3)
#     for label in ax.xaxis.get_ticklabels():
#         label.set_rotation(30)
#     handles,labels = ax.get_legend_handles_labels()

#     ax.legend(handles, labels = label_uniq, loc='upper right', bbox_to_anchor=(1, 1),borderaxespad=0)
# #     pp = PdfPages('iters.pdf')
# #     pp.savefig(fig)
# #     pp.close()
#     plt.show()
#     return fig

    
# def PlotPCoA(PCoA,label_refer,label_base):

#     fig = plt.figure(1,figsize=(5,5), dpi= 50)

#     plt.title('PCoA')

#     plt.xlabel('PC1')
#     plt.ylabel('PC2')
#     xValue = np.array(PCoA[:,0])
#     yValue = np.array(PCoA[:,1])

#     label_uniq = ['Single','TP','FP','Undec','Sim']
    
#     ####label_refer为预测标签，0为single；1为sim；0-1为doublet
#     ####label_base为真实标签，celltype;DOUB;DOUB_sim
#     label_refer[label_refer==0.05] = 0
#     plt.scatter(xValue[np.all([label_refer==0, label_base !='DOUB'],axis =0)], yValue[np.all([label_refer==0, label_base !='DOUB'],axis =0)], s=1, c = '#FFE6D6', marker='o')
#     plt.scatter(xValue[np.all([label_refer!=0, label_refer!=1, label_base =='DOUB'],axis =0)], yValue[np.all([label_refer!=0, label_refer!=1, label_base =='DOUB'],axis =0)], s=1, c='#C45462', marker='o')
#     plt.scatter(xValue[np.all([label_refer!=0, label_refer!=1, label_base !='DOUB'],axis =0)], yValue[np.all([label_refer!=0, label_refer!=1, label_base !='DOUB'],axis =0)], s=1, c='#F49675', marker='o')
#     plt.scatter(xValue[np.all([label_refer==0, label_base =='DOUB'],axis =0)], yValue[np.all([label_refer==0, label_base =='DOUB'],axis =0)], s=1, c='#00B4FF', marker='o')
#     plt.scatter(xValue[label_refer==1], yValue[label_refer==1], s=1, c='#2F4858', marker='o')

#     ax = fig.gca()
    
#     ax.patch.set_facecolor("#FFFFFF") 
#     ax.patch.set_alpha(0.3) 
#     ax.grid(color='r',
#             linestyle='--',
#             linewidth=1,
#             alpha=0.3)
#     for label in ax.xaxis.get_ticklabels():
#         label.set_rotation(30)
#     handles,labels = ax.get_legend_handles_labels()
#     ax.legend(handles, labels = label_uniq, loc='upper right', bbox_to_anchor=(1.1, 1),borderaxespad=0)
#     plt.show()

    
    
    
def CallDoublet_PCA(iteration, mat, sim, jac, core, npc, k, n_tree,label_refer,labelmat):
    time1 = time.time()
    mat_concat, label_concat = concat_set(iteration, jac, mat,sim,label_refer)
    ####label_concat只包含{0:unlabel; 0-1:predicted doublets; 1:simulated doublets}
    ###mat_concat所有细胞+sim
    ###label_concat 带标记的raw(0 /3)+ sim(1) 
    
    distance_inuse = GetUMAP(mat_concat)
    ###k为构建knn图时所用的邻居节点default
    knn,distance = get_knn_graph(distance_inuse, k, n_tree)
    
    ####3 parts: 1. 未标记unlabel； 2.已标记detected; 3.simulated 
    cells_unlabel_knn = knn[label_concat==0]
    cells_sim_knn = knn[label_concat==1]
    cells_detected_knn = knn[np.all([label_concat!=0, label_concat !=1],axis =0)]
    
    cells_unlabel_distance = distance[label_concat==0]
    cells_sim_distance = distance[label_concat==1]
    
    cells_detected_distance = distance[np.all([label_concat!=0, label_concat !=1],axis =0)]
    ###上一轮unlabel的idx
    cell_unlabel_idx = np.where(label_concat==0)[0]
    ###上一轮detected doublet的idx
    cell_detected_idx = np.where((label_concat!=0)&(label_concat!=1))[0]
    
    label_raw = label_concat[label_concat!=1]
    score_unlabel = CalKNNScore(cells_unlabel_knn,cells_unlabel_distance,np.min(distance),np.max(distance),label_concat)
    score_sim = CalKNNScore(cells_sim_knn,cells_sim_distance,np.min(distance),np.max(distance),label_concat)
    
    score_unlabel = pd.DataFrame(score_unlabel)
    score_unlabel.index = cell_unlabel_idx
    
 
    ###返回pd格式 的unlabelled cell score
    if len(cells_detected_knn) == 0:
        score_detected = pd.DataFrame(np.array([]))
    else:
        score_detected = CalKNNScore(cells_detected_knn,cells_detected_distance,np.min(distance),np.max(distance),label_concat)
        score_detected = pd.DataFrame(score_detected)
        score_detected.index = cell_detected_idx 

    return score_unlabel,score_sim,score_detected

def CallDoublet_PCoA(iteration, mat, sim, jac, core, npc, k, n_tree,label_refer,labelmat):
    mat_concat, label_concat = concat_set(iteration, jac, mat,sim,label_refer)
    ####label_concat只包含{0:unlabel; 0-1:predicted doublets; 1:simulated doublets}
    ###mat_concat所有细胞+sim
    ####基于重构的mat，计算原始细胞与新的doublet之间的关系
    ###只需计算新的simulated 与 raw 的jaccard
    mat_jaccard = ReshapeJac(mat_concat, label_concat, jac, core)
#     mat_jaccard = CalJac(mat_concat)
    ##PCoA_ased
    mat_inuse = PCoA(mat_jaccard,npc)

    ###k为构建knn图时所用的邻居节点default
    knn,distance = get_knn_graph(mat_inuse, k, n_tree)
    
    ####3 parts: 1. 未标记unlabel； 2.已标记detected; 3.simulated 
    cells_unlabel_knn = knn[label_concat==0]
    cells_sim_knn = knn[label_concat==1]
    cells_detected_knn = knn[np.all([label_concat!=0, label_concat !=1],axis =0)]
    
    cells_unlabel_distance = distance[label_concat==0]
    cells_sim_distance = distance[label_concat==1]
    
    cells_detected_distance = distance[np.all([label_concat!=0, label_concat !=1],axis =0)]
    ###上一轮unlabel的idx
    cell_unlabel_idx = np.where(label_concat==0)[0]
    ###上一轮detected doublet的idx
    cell_detected_idx = np.where((label_concat!=0)&(label_concat!=1))[0]
    
    label_raw = label_concat[label_concat!=1]
    score_unlabel = CalKNNScore(cells_unlabel_knn,cells_unlabel_distance,np.min(distance),np.max(distance),label_concat)
    score_sim = CalKNNScore(cells_sim_knn,cells_sim_distance,np.min(distance),np.max(distance),label_concat)
    
    score_unlabel = pd.DataFrame(score_unlabel)
    score_unlabel.index = cell_unlabel_idx
    
    ###返回pd格式 的unlabelled cell score
    if len(cells_detected_knn) == 0:
        score_detected = pd.DataFrame(np.array([]))
    else:
        score_detected = CalKNNScore(cells_detected_knn,cells_detected_distance,np.min(distance),np.max(distance),label_concat)
        score_detected = pd.DataFrame(score_detected)
        score_detected.index = cell_detected_idx
    time5 = time.time()
    print('Scoring: ',time5-time4)
    return score_unlabel,score_sim,score_detected



def GetStrategy(mat,jac):
#     PCA_umap = GetUMAP(scipy.sparse.csr_matrix(mat))
    PCA_umap = GetUMAP(mat)
#     PCA_umap = MinMaxScaler().fit_transform(PCA_umap) 
    PCA_umap = StandardScaler().fit_transform(PCA_umap) 
#     distance = metrics.pairwise_distances(PCA_umap, Y=None, metric='euclidean', n_jobs=64)
    ##adjusted
    Parameters ={'eps':0.1,'min_samples':int(mat.shape[0]*0.005)}
    db1 = DBSCAN(eps = Parameters['eps'], min_samples = Parameters['min_samples'], metric = 'euclidean').fit(PCA_umap) 
    cl_PCA = db1.labels_

    try:
        chs1 = calinski_harabasz_score(PCA_umap,cl_PCA)
    except:
        chs1 = 0.0001
    PCoAX = PCoA(jac, 30)
    reducer = umap.UMAP(random_state=42, n_components = 2,min_dist=0.5)
    UMAP_PCoA = reducer.fit_transform(PCoAX)
    UMAP_PCoA = StandardScaler().fit_transform(UMAP_PCoA) 
    db2 = DBSCAN(eps = Parameters['eps'], min_samples = Parameters['min_samples'], metric = 'euclidean').fit(UMAP_PCoA) 
    cl_PCoA = db2.labels_
    try:
        chs2 = calinski_harabasz_score(PCA_umap,cl_PCoA)
    except:
        chs2 = 0.0001
    if chs1 > chs2:
        divides = chs1/chs2
    else:
        divides= chs2/chs1
    if divides<1.5:
        strategy = 'PCA'
    else:
        strategy = 'PCoA'
    print(strategy)
    return strategy

###define main function here
class KNNIter(object):
    def __init__(self, rawmat, strategy = None, core = None, sim_rate= None, nPC=None, neigbors=None, nTree=None, label= None, exprate = None):
        if isinstance(rawmat, anndata.AnnData):
            self.rawmat_ann = rawmat
            rawmat = rawmat.X 
        else:
            self.rawmat_ann = None
            rawmat = rawmat   
        ##preprocess1, select the peaks that are open in more than 1% cells
        ## transfer to scipy
        rawmat = scipy.sparse.csr_matrix(rawmat).astype(np.uint8)
        rawmat_pf = rawmat[:,(rawmat.sum(axis = 0)> 0.01* rawmat.shape[0]).A.squeeze()]
        print(rawmat_pf.shape)
        self.rawmat = rawmat
        self.fmat = rawmat_pf
        
        ###csr matrix
        self.core_max = mp.cpu_count()
        
        if core is None:
            self.core = self.core_max
        elif core > self.core_max:
            self.core = self.core_max
        else:
            self.core = core
        if sim_rate is None:
            self.sim_rate = 0.3
        else:
            self.sim_size = sim_size
        if nPC is None:
            self.nPC = 5
        else:
            self.nPC = nPC
        if neigbors is None:
            self.neigbors = 40
        else:
            self.neigbors = neigbors
        if nTree is None:
            self.nTree = 30
        else:
            self.nTree = nTree   
        if label is None:
            self.label = np.array(range(self.rawmat.shape[0]))
        else:
            self.label = label
        if exprate is None:
            self.exprate = 0.1
        else:
            self.exprate = exprate
        ##jaccard distance mat of raw set
        time0 = time.time()
#         self.jac = metrics.pairwise_distances(self.rawmat.A, Y=None, metric='jaccard', n_jobs=self.core)
        self.jac = CalJac(self.rawmat)

        time1 = time.time()
        print('raw jaccard cal: ',time1-time0)
        if strategy is None:
            print('evaluating the input data!')
            self.strategy = GetStrategy(self.rawmat,self.jac)
            time2 = time.time()
            print('Evaluate time: ',time2-time1)
        else:
            self.strategy = strategy
    def IterCall(self, mat=None, strategy=None, core=None, simrate=None, npc=None, k=None, n_tree=None, labelmat=None, exprate=None):
# #     def IterCall(mat=self.rawmat, strategy=self.strategy, core=self.core, simrate=self.sim_rate, npc=self.nPC, k=self.neigbors, n_tree=self.nTree, labelmat=self.label, exprate=self.exprate):
#     def IterCall(self):
        if mat is None:
            mat = self.fmat 
        if isinstance(mat, np.ndarray):
            mat = scipy.sparse.csr_matrix(mat)
        if strategy is None:
            strategy = self.strategy
        if core is None:
            core = self.core_max
        if simrate is None:
            simrate = self.sim_rate
        if npc is None:
            npc = self.nPC
        if k is None:
            k = self.neigbors
        if n_tree is None:
            n_tree = self.nTree
        if labelmat is None:
            labelmat = self.label
        if exprate is None:
            exprate = self.exprate

        time_start = time.time()
        ###原始数据的jaccard
        jaccard_raw = self.jac
        ###参与KNN计算
        label_refer = np.zeros(mat.shape[0]).astype(np.float64)
        ###控制迭代
        loopproba = 10
        iteration = 1
        ###保存每一轮的分数
        score_pool = np.array([]).reshape(-1,mat.shape[0])
        random.seed(2222)
        ####
        while random.sample([1]*loopproba+[0]*(10-loopproba),1)[0]:
            print('Iteration',iteration, 'start')
            # time_iter_start = time.time()
            if strategy == 'PCA':
                score_unlabel, score_sim, score_detected = CallDoublet_PCA(iteration, mat, simrate, jaccard_raw, core, npc, k, n_tree,label_refer,labelmat)

            else:
                score_unlabel, score_sim, score_detected = CallDoublet_PCoA(iteration, mat, simrate, jaccard_raw, core, npc, k, n_tree,label_refer,labelmat)
    
            ##当前轮次的分数
            score_merge = np.zeros(score_pool.shape[1])
            ###填入当前轮 unlabel的doublet scores,
            score_merge[score_unlabel.index] = score_unlabel.values.squeeze()
            ###填入当前轮 unlabel的doublet scores
            score_merge[score_detected.index] = score_detected.values.squeeze()
            score_pool = np.vstack((score_pool,score_merge))

            thresh_accept = Calsigmathresh(score_unlabel.values.reshape(-1),score_sim)[2] 
            ####交点+2sigma（自适应）
            
#             try:
#                 reject_idx = []
# #                 thresh_reject = np.quantile(score_detected,0.1,interpolation='higher')
# #                 reject_idx = score_detected[score_detected[0]<thresh_reject].index

#             except:
#                 reject_idx = []
#             #####option reject criteria : the lower quantile of sim
            
#             unlabel_neig = stac_nei[score_unlabel.index]
#             predict_idx = score_unlabel[np.all([score_unlabel[0]>thresh_accept,unlabel_neig<int(0.6*self.neigbors)],axis=0)].index
            predict_idx = score_unlabel[score_unlabel[0]>thresh_accept].index


            ####没有高于阈值的，开启退火
            if len(predict_idx) == 0:
                if random.sample([1]*loopproba+[0]*(10-loopproba),1)[0]:
                    iteration += 1
                    loopproba -= int(np.ceil(5/iteration))
                    continue
                else:
                    break
            else:
                scalar2 = MinMaxScaler(feature_range=(0.1,0.9))
                predict_idx_score = scalar2.fit_transform(np.array(score_unlabel[score_unlabel[0]>thresh_accept])).reshape(1,-1)
#                 predict_idx_score = scalar.fit_transform(np.array(score_unlabel[np.all([score_unlabel[0]>thresh_accept,unlabel_neig<int(0.6*self.neigbors)],axis=0)])).reshape(1,-1)

                 ###这一轮迭代的结果传给下一轮
                label_refer[predict_idx] = predict_idx_score
 
            
            #####labelrefer不为0(从未被预测为doublets)和0.05(被预测为doublets后又被reject)
            pred_doub = label_refer[np.all([label_refer!=0,label_refer!=0.05],axis=0)]
#             print(len(pred_doub),'cells have been detected')
            ####加上这轮预测的减去这轮拒绝后的数量达到expect，开启退火
            if len(pred_doub) > int(exprate * mat.shape[0]*0.9):
                loopproba -= int(np.ceil(5/iteration))  
            # print('Iteraion', iteration, 'elapse', time.time()-time_iter_start)
            iteration += 1
        
        proba_final = np.mean(score_pool,axis = 0)
        min_max_scaler = MinMaxScaler( )
        proba_final = min_max_scaler.fit_transform(proba_final.reshape(-1, 1)).squeeze()
        thresh_final = np.quantile(proba_final,1-exprate,interpolation= 'higher')
        print('Done')
        time_final = time.time()
        print('ALL time:', time_final - time_start)
        if self.rawmat_ann is None:
            return labelmat[proba_final > thresh_final],proba_final
        else:
            pred_results = np.zeros(proba_final.shape[0])
            pred_results[proba_final > thresh_final] = 1
            self.rawmat_ann.obs['PredDBL'] = pred_results
            self.rawmat_ann.obs['DBLscore'] = proba_final
            return self.rawmat_ann
        
