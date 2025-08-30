# /home/users/nus/e1503317/projects/dmeng/zhlin/SMOBench/Utils/SMOBench_clustering.py
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri

def universal_clustering(adata, n_clusters, used_obsm, method='kmeans', key='clusters', 
                         use_pca=False, n_comps=20, start=0.1, end=3.0, increment=0.01):
    """
    Universal clustering function supporting four methods: kmeans, mclust, leiden, and louvain.

    Parameters:
    adata: AnnData object containing the embedding features to be clustered
    n_clusters: Target number of clusters
    used_obsm: Key in adata.obsm where the embedding features are stored (e.g., 'merged_emb')
    method: Clustering method, optional values: 'kmeans', 'mclust', 'leiden', 'louvain'
    key: Key in adata.obs where the clustering results will be stored
    use_pca: Whether to perform PCA dimensionality reduction on the embedding features
    n_comps: Number of dimensions after PCA reduction (only effective when use_pca=True)
    start: Starting value for resolution search (only effective for leiden/louvain)
    end: Ending value for resolution search (only effective for leiden/louvain)
    increment: Step size for resolution search (only effective for leiden/louvain)

    Returns:
    adata: AnnData object containing the clustering results
    eg:adata.obs['leiden']: Leiden clustering results
    """
    if use_pca:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_comps, random_state=0)
        adata.obsm[f'{used_obsm}_pca'] = pca.fit_transform(adata.obsm[used_obsm])
        used_rep = f'{used_obsm}_pca'
    else:
        used_rep = used_obsm
    
    if method == 'kmeans':
        kmeans = KMeans(n_clusters=n_clusters, random_state=0)
        adata.obs[key] = kmeans.fit_predict(adata.obsm[used_rep]).astype(str)
    
    elif method == 'mclust':
        adata = mclust_R(adata, num_cluster=n_clusters, used_obsm=used_rep)
        adata.obs[key] = adata.obs['mclust'].astype(str)
    
    elif method in ['leiden', 'louvain']:
        sc.pp.neighbors(adata, n_neighbors=50, use_rep=used_rep)
        
        res = search_res(
            adata, 
            n_clusters=n_clusters, 
            method=method, 
            start=start, 
            end=end, 
            increment=increment
        )
        
        if method == 'leiden':
            sc.tl.leiden(adata, random_state=0, resolution=res)
            adata.obs[key] = adata.obs['leiden'].astype(str)
        else:
            sc.tl.louvain(adata, random_state=0, resolution=res)
            adata.obs[key] = adata.obs['louvain'].astype(str)
    
    else:
        raise ValueError(f"不支持的聚类方法: {method}，可选方法为kmeans、mclust、leiden、louvain")
    
    return adata


def mclust_R(adata, num_cluster, used_obsm='emb', random_seed=2020):
    """调用R的mclust进行聚类（修复类别类型转换错误）"""
    np.random.seed(random_seed)
    robjects.r.library("mclust")
    numpy2ri.activate()  # 启用numpy到R数组的转换
    
    # 设置随机种子
    robjects.r['set.seed'](random_seed)
    # 调用R的Mclust函数
    res = robjects.r['Mclust'](numpy2ri.numpy2rpy(adata.obsm[used_obsm]), num_cluster)
    # 提取聚类结果
    mclust_labels = np.array(res[-2])  # res[-2]是聚类标签
    
    adata.obs['mclust'] = pd.Categorical(mclust_labels.astype(int))
    return adata


def search_res(adata, n_clusters, method='leiden', start=0.1, end=3.0, increment=0.01):
    print(f"Searching resolution for {method}...")
    found = False
    for res in sorted(np.arange(start, end, increment), reverse=True):
        if method == 'leiden':
            sc.tl.leiden(adata, random_state=0, resolution=res, key_added='tmp_clust')
            cluster_count = adata.obs['tmp_clust'].nunique()
        else:
            sc.tl.louvain(adata, random_state=0, resolution=res, key_added='tmp_clust')
            cluster_count = adata.obs['tmp_clust'].nunique()
        
        print(f"resolution={res:.3f}, cluster number={cluster_count}")
        if cluster_count == n_clusters:
            found = True
            break
    
    if not found:
        raise ValueError(f"未找到满足聚类数{n_clusters}的分辨率，请调整搜索范围或步长")
    
    return res
