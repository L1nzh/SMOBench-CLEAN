import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import warnings
warnings.filterwarnings('ignore')

def universal_clustering(adata, n_clusters, used_obsm, method='kmeans', key='clusters', 
                         use_pca=False, n_comps=20, start=0.1, end=3.0, increment=0.01,
                         random_state=2024):
    """
    Universal clustering function supporting multiple clustering methods for SMOBench
    
    Parameters:
        adata: AnnData object containing integration embeddings
        n_clusters: Target number of clusters
        used_obsm: Key for embedding in adata.obsm (e.g., 'merged_emb', 'spatial_emb')
        method: Clustering method ('kmeans', 'mclust', 'leiden', 'louvain')
        key: Output key for storing clusters in adata.obs
        use_pca: Whether to apply PCA preprocessing to embeddings
        n_comps: Number of PCA components (only used if use_pca=True)
        start: Starting resolution for leiden/louvain search
        end: Ending resolution for leiden/louvain search  
        increment: Resolution search step size
        random_state: Random seed for reproducibility
    
    Returns:
        adata: Updated AnnData object with clustering results in adata.obs[key]
    """
    
    # Validate inputs
    if used_obsm not in adata.obsm.keys():
        raise ValueError(f"Embedding key '{used_obsm}' not found in adata.obsm. Available keys: {list(adata.obsm.keys())}")
    
    if method not in ['kmeans', 'mclust', 'leiden', 'louvain']:
        raise ValueError(f"Unsupported method '{method}'. Choose from: kmeans, mclust, leiden, louvain")
    
    # Set random seeds for reproducibility
    np.random.seed(random_state)
    
    # Apply PCA preprocessing if requested
    if use_pca:
        print(f"Applying PCA preprocessing: {adata.obsm[used_obsm].shape[1]} -> {n_comps} dimensions")
        pca = PCA(n_components=n_comps, random_state=random_state)
        adata.obsm[f'{used_obsm}_pca'] = pca.fit_transform(adata.obsm[used_obsm])
        embedding_key = f'{used_obsm}_pca'
    else:
        embedding_key = used_obsm
    
    # Perform clustering based on method
    if method == 'kmeans':
        print(f"Running K-means clustering with {n_clusters} clusters...")
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
        cluster_labels = kmeans.fit_predict(adata.obsm[embedding_key])
        adata.obs[key] = cluster_labels.astype(str)
        
    elif method == 'mclust':
        print(f"Running mclust clustering with {n_clusters} clusters...")
        try:
            adata = _mclust_R(adata, num_cluster=n_clusters, used_obsm=embedding_key, random_seed=random_state)
            adata.obs[key] = adata.obs['mclust'].astype(str)
        except Exception as e:
            print(f"mclust failed: {e}")
            print("Falling back to K-means clustering...")
            kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
            cluster_labels = kmeans.fit_predict(adata.obsm[embedding_key])
            adata.obs[key] = cluster_labels.astype(str)
    
    elif method in ['leiden', 'louvain']:
        print(f"Running {method} clustering with {n_clusters} clusters...")
        
        # Build neighborhood graph if not exists
        print("Building neighborhood graph...")
        sc.pp.neighbors(adata, use_rep=embedding_key, n_neighbors=50, random_state=random_state)
        
        # Search for optimal resolution
        resolution = _search_resolution(
            adata, 
            n_clusters=n_clusters, 
            method=method, 
            start=start, 
            end=end, 
            increment=increment,
            random_state=random_state
        )
        
        # Apply clustering with found resolution
        if method == 'leiden':
            sc.tl.leiden(adata, resolution=resolution, random_state=random_state, key_added=key)
        else:  # louvain
            sc.tl.louvain(adata, resolution=resolution, random_state=random_state, key_added=key)
        
        # Convert to string for consistency
        adata.obs[key] = adata.obs[key].astype(str)
    
    print(f"Clustering complete. Found {adata.obs[key].nunique()} clusters (target: {n_clusters})")
    return adata


def _mclust_R(adata, num_cluster, used_obsm='emb', random_seed=2024):
    """
    R mclust clustering implementation via rpy2
    """
    try:
        import rpy2.robjects as robjects
        from rpy2.robjects import numpy2ri
        
        # Load R mclust library
        robjects.r.library("mclust")
        
        # Enable numpy to R conversion
        numpy2ri.activate()
        
        # Set R random seed
        robjects.r['set.seed'](random_seed)
        
        # Get mclust function
        rmclust = robjects.r['Mclust']
        
        # Run clustering
        embedding_data = adata.obsm[used_obsm]
        res = rmclust(numpy2ri.numpy2rpy(embedding_data), num_cluster)
        
        # Extract clustering results
        mclust_labels = np.array(res[-2])  # Cluster assignments are in res[-2]
        
        # Store results
        adata.obs['mclust'] = mclust_labels.astype(int).astype('category')
        
        return adata
        
    except ImportError:
        raise ImportError("rpy2 not installed. Install with: conda install -c conda-forge rpy2")
    except Exception as e:
        raise RuntimeError(f"R mclust clustering failed: {e}")


def _search_resolution(adata, n_clusters, method='leiden', start=0.1, end=3.0, increment=0.01, random_state=2024):
    """
    Search for optimal resolution to achieve target cluster number
    """
    print(f"Searching resolution for {method} clustering (target: {n_clusters} clusters)...")
    
    # Search from high to low resolution for more stable results
    resolutions = sorted(np.arange(start, end, increment), reverse=True)
    
    for resolution in resolutions:
        # Apply clustering with current resolution
        if method == 'leiden':
            sc.tl.leiden(adata, resolution=resolution, random_state=random_state, key_added='tmp_search')
            n_found = adata.obs['tmp_search'].nunique()
        else:  # louvain
            sc.tl.louvain(adata, resolution=resolution, random_state=random_state, key_added='tmp_search')
            n_found = adata.obs['tmp_search'].nunique()
        
        print(f"Resolution {resolution:.3f}: {n_found} clusters")
        
        if n_found == n_clusters:
            print(f"Found optimal resolution: {resolution:.3f}")
            return resolution
        
        # Clean up temporary results
        if 'tmp_search' in adata.obs.columns:
            del adata.obs['tmp_search']
    
    # If exact match not found, return resolution closest to target
    print(f"Exact cluster number not found. Using resolution {resolution:.3f} with {n_found} clusters")
    return resolution


def batch_clustering(adata, n_clusters, used_obsm, methods=['leiden', 'louvain', 'kmeans', 'mclust'], 
                     prefix='', **kwargs):
    """
    Apply multiple clustering methods to the same embedding
    
    Parameters:
        adata: AnnData object
        n_clusters: Target cluster number
        used_obsm: Embedding key in adata.obsm  
        methods: List of clustering methods to apply
        prefix: Prefix for output keys (e.g., 'SpatialGlue_')
        **kwargs: Additional arguments passed to universal_clustering
    
    Returns:
        adata: AnnData with clustering results for all methods
    """
    
    for method in methods:
        key = f"{prefix}{method}" if prefix else method
        print(f"\n=== Applying {method} clustering ===")
        
        try:
            adata = universal_clustering(
                adata, 
                n_clusters=n_clusters,
                used_obsm=used_obsm, 
                method=method,
                key=key,
                **kwargs
            )
            print(f"✅ {method} clustering completed -> adata.obs['{key}']")
        except Exception as e:
            print(f"❌ {method} clustering failed: {e}")
    
    return adata


def evaluate_clustering_quality(adata, cluster_key, embedding_key=None):
    """
    Calculate clustering quality metrics
    
    Parameters:
        adata: AnnData object with clustering results
        cluster_key: Key for cluster labels in adata.obs
        embedding_key: Key for embedding in adata.obsm (optional)
    
    Returns:
        dict: Dictionary of quality metrics
    """
    
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    # Basic cluster statistics
    cluster_labels = adata.obs[cluster_key]
    n_clusters = cluster_labels.nunique()
    cluster_sizes = cluster_labels.value_counts().sort_index()
    
    metrics = {
        'n_clusters': n_clusters,
        'min_cluster_size': cluster_sizes.min(),
        'max_cluster_size': cluster_sizes.max(),
        'mean_cluster_size': cluster_sizes.mean(),
        'std_cluster_size': cluster_sizes.std()
    }
    
    # Silhouette score if embedding provided
    if embedding_key and embedding_key in adata.obsm.keys():
        try:
            silhouette = silhouette_score(adata.obsm[embedding_key], cluster_labels.astype(int))
            metrics['silhouette_score'] = silhouette
        except Exception as e:
            print(f"Could not calculate silhouette score: {e}")
            metrics['silhouette_score'] = np.nan
    
    return metrics


# Convenience functions for specific use cases
def spatialglue_clustering(adata, n_clusters, **kwargs):
    """Clustering wrapper for SpatialGlue results"""
    return batch_clustering(adata, n_clusters, used_obsm='spatial_emb', prefix='SpatialGlue_', **kwargs)

def spamosaic_clustering(adata, n_clusters, **kwargs):
    """Clustering wrapper for SpaMosaic results"""  
    return batch_clustering(adata, n_clusters, used_obsm='merged_emb', prefix='SpaMosaic_', **kwargs)

def praga_clustering(adata, n_clusters, **kwargs):
    """Clustering wrapper for PRAGA results"""
    return batch_clustering(adata, n_clusters, used_obsm='PRAGA_emb', prefix='PRAGA_', **kwargs)
