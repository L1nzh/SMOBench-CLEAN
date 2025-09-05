import os
import torch
import pandas as pd
import scanpy as sc
import numpy as np
import argparse
import time
import sys
import re
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans

# Add project root directory to module search path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.append(project_root)

# Add SpaMultiVAE to path
spamultivae_path = os.path.join(project_root, "Methods/SpaMultiVAE")
sys.path.append(spamultivae_path)

from spaMultiVAE import SPAMULTIVAE
from preprocess import normalize, geneSelection
from Utils.SMOBench_clustering import universal_clustering

# Custom SPAMULTIVAE class with robust jitter handling
class RobustSPAMULTIVAE(SPAMULTIVAE):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Increase jitter for numerical stability
        self.svgp.jitter = 1e-4  # Increased from default 1e-8


def parse_dataset_info(args):
    """
    Extract dataset_name and subset_name from RNA_path or save_path
    Support two modes:
    1. Manual specification --dataset Human_Lymph_Nodes/A1
    2. Auto extraction from paths
    """
    if hasattr(args, 'dataset') and args.dataset:
        parts = args.dataset.strip('/').split('/')
        if len(parts) == 2:
            return parts[0], parts[1]
        elif len(parts) == 1:
            return parts[0], "Unknown"
    
    # Auto parse RNA_path
    match = re.search(r'Dataset/([^/]+)/([^/]+)/([^/]+)/adata_RNA\.h5ad', args.RNA_path)
    if match:
        return match.group(2), match.group(3)
    return "Unknown", "Unknown"


def main(args):
    device = torch.device(args.device if torch.cuda.is_available() else 'cpu')
    print('device:', device)
    
    # === Load Data ===
    adata_rna = sc.read_h5ad(args.RNA_path)
    adata_rna.var_names_make_unique()

    if args.ADT_path:
        adata_adt = sc.read_h5ad(args.ADT_path)
        modality = 'ADT'
        modality_name = 'Proteome'
        adata_adt.var_names_make_unique()
        second_adata = adata_adt
    elif args.ATAC_path:
        adata_atac = sc.read_h5ad(args.ATAC_path)
        modality = 'ATAC'
        modality_name = 'Epigenome'
        adata_atac.var_names_make_unique()
        second_adata = adata_atac
    else:
        raise ValueError("Either ADT_path or ATAC_path must be provided.")

    print(f"Processing {args.data_type}: RNA + {modality_name} integration...")

    # === Ensure both datasets have same cells ===
    common_obs = adata_rna.obs_names.intersection(second_adata.obs_names)
    print(f"Common cells: {len(common_obs)}, RNA cells: {adata_rna.n_obs}, {modality_name} cells: {second_adata.n_obs}")
    adata_rna = adata_rna[common_obs].copy()
    second_adata = second_adata[common_obs].copy()
    print(f"After intersection - RNA: {adata_rna.n_obs}, {modality_name}: {second_adata.n_obs}")

    # === Extract data matrices and spatial coordinates ===
    x1 = adata_rna.X.toarray() if hasattr(adata_rna.X, 'toarray') else adata_rna.X
    x2 = second_adata.X.toarray() if hasattr(second_adata.X, 'toarray') else second_adata.X
    
    # Get spatial coordinates
    if 'spatial' in adata_rna.obsm:
        loc = adata_rna.obsm['spatial'].copy()
    elif 'spatial' in second_adata.obsm:
        loc = second_adata.obsm['spatial'].copy()
    else:
        raise ValueError("No spatial coordinates found in obsm['spatial']")
    
    # Ensure data types
    x1 = x1.astype('float64')
    x2 = x2.astype('float64')
    loc = loc.astype('float64')
    
    print(f"Data shapes - RNA: {x1.shape}, {modality_name}: {x2.shape}, Location: {loc.shape}")

    # === Preprocessing ===
    print("Preprocessing data...")
    
    # Auto batch size setting
    if x1.shape[0] <= 1024:
        batch_size = 128
    elif x1.shape[0] <= 2048:
        batch_size = 256
    else:
        batch_size = 512
    
    # Scale spatial coordinates
    scaler = MinMaxScaler()
    loc_range = 20.
    loc = scaler.fit_transform(loc) * loc_range
    
    # Create inducing points using K-means for better stability
    n_inducing = min(200, max(100, loc.shape[0] // 20))  # Adaptive number of inducing points
    print(f"Using {n_inducing} inducing points based on dataset size")
    
    from sklearn.cluster import KMeans
    loc_kmeans = KMeans(n_clusters=n_inducing, n_init=10, random_state=args.seed).fit(loc)
    initial_inducing_points = loc_kmeans.cluster_centers_
    print(f"Inducing points shape: {initial_inducing_points.shape}")
    
    # Create and normalize AnnData objects
    adata1 = sc.AnnData(x1, dtype="float64")
    adata1 = normalize(adata1,
                      size_factors=True,
                      normalize_input=True,
                      logtrans_input=True)

    adata2 = sc.AnnData(x2, dtype="float64")
    adata2 = normalize(adata2,
                      size_factors=False,
                      normalize_input=True,
                      logtrans_input=True)

    # For protein background prior calculation
    if modality == 'ADT':
        adata2_no_scale = sc.AnnData(x2, dtype="float64")
        adata2_no_scale = normalize(adata2_no_scale,
                          size_factors=False,
                          normalize_input=False,
                          logtrans_input=True)

        # Fit GMM model to the protein counts for background prior
        gm = GaussianMixture(n_components=2, covariance_type="diag", n_init=20).fit(adata2_no_scale.X)
        back_idx = np.argmin(gm.means_, axis=0)
        protein_log_back_mean = np.log(np.expm1(gm.means_[back_idx, np.arange(adata2_no_scale.n_vars)]))
        protein_log_back_scale = np.sqrt(gm.covariances_[back_idx, np.arange(adata2_no_scale.n_vars)])
        print("protein_back_mean shape", protein_log_back_mean.shape)
    else:
        # For ATAC data, use zero background
        protein_log_back_mean = np.zeros(adata2.n_vars)
        protein_log_back_scale = np.ones(adata2.n_vars)

    # === Initialize SpaMultiVAE model ===
    print("Initializing SpaMultiVAE model...")
    model = RobustSPAMULTIVAE(
        gene_dim=adata1.n_vars, 
        protein_dim=adata2.n_vars, 
        GP_dim=2, 
        Normal_dim=18,
        encoder_layers=[128, 64], 
        gene_decoder_layers=[128], 
        protein_decoder_layers=[128],
        gene_noise=0, 
        protein_noise=0, 
        encoder_dropout=0, 
        decoder_dropout=0,
        fixed_inducing_points=True, 
        initial_inducing_points=initial_inducing_points,
        fixed_gp_params=False, 
        kernel_scale=10.0,  # Reduced kernel scale for better stability 
        N_train=adata1.n_obs, 
        KL_loss=0.025,
        dynamicVAE=True,
        init_beta=10, 
        min_beta=4, 
        max_beta=25, 
        protein_back_mean=protein_log_back_mean, 
        protein_back_scale=protein_log_back_scale, 
        dtype=torch.float64,
        device=args.device
    )
    
    print(str(model))
    
    # === Train SpaMultiVAE model ===
    print("Training SpaMultiVAE model...")
    start_time = time.time()
    
    model.train_model(
        pos=loc, 
        gene_ncounts=adata1.X, 
        gene_raw_counts=adata1.raw.X, 
        gene_size_factors=adata1.obs.size_factors,
        protein_ncounts=adata2.X, 
        protein_raw_counts=adata2.raw.X,
        lr=5e-3, 
        weight_decay=1e-6, 
        batch_size=batch_size, 
        num_samples=1,
        train_size=0.95, 
        maxiter=5000, 
        patience=200, 
        save_model=False
    )
    
    end_time = time.time()
    print('Training time:', end_time - start_time)

    # === Get embeddings ===
    print("Extracting latent embeddings...")
    final_latent = model.batching_latent_samples(X=loc, gene_Y=adata1.X, protein_Y=adata2.X, batch_size=batch_size)
    
    # === Spatial resolution enhancement (optional) ===
    if args.enhance_resolution:
        print("Enhancing spatial resolution...")
        neigh = NearestNeighbors(n_neighbors=2).fit(loc)
        nearest_dist = neigh.kneighbors(loc, n_neighbors=2)[0]
        small_distance = np.median(nearest_dist[:,1])/4
        
        loc_new1 = loc.copy()
        loc_new2 = loc.copy()
        loc_new3 = loc.copy()
        loc_new4 = loc.copy()
        
        loc_new1[:,0] = loc_new1[:,0] - small_distance
        loc_new1[:,1] = loc_new1[:,1] + small_distance
        loc_new2[:,0] = loc_new2[:,0] + small_distance
        loc_new2[:,1] = loc_new2[:,1] + small_distance
        loc_new3[:,0] = loc_new3[:,0] - small_distance
        loc_new3[:,1] = loc_new3[:,1] - small_distance
        loc_new4[:,0] = loc_new4[:,0] + small_distance
        loc_new4[:,1] = loc_new4[:,1] - small_distance
        
        loc_enhance = np.concatenate((loc_new1, loc_new2, loc_new3, loc_new4, loc), axis=0)
        
        enhanced_latent, _, _ = model.batching_predict_samples(
            X_test=loc_enhance, 
            X_train=loc, 
            gene_Y_train=adata1.X, 
            protein_Y_train=adata2.X, 
            batch_size=batch_size, 
            n_samples=25
        )
        
        # Create enhanced AnnData with both original and enhanced data
        enhanced_obs_names = [f"{name}_enhanced_{i}" for i in range(4) for name in common_obs] + list(common_obs)
        
        # Create enhanced spatial coordinates (convert back from scaled)
        loc_enhance_orig = loc_enhance / loc_range
        loc_enhance_orig = scaler.inverse_transform(loc_enhance_orig)
        
        adata_enhanced = sc.AnnData(enhanced_latent, dtype="float64")
        adata_enhanced.obs_names = enhanced_obs_names
        adata_enhanced.obsm['SpaMultiVAE'] = enhanced_latent
        adata_enhanced.obsm['spatial'] = loc_enhance_orig
        
        # Use enhanced data for clustering
        adata = adata_enhanced.copy()
        
    else:
        # === Build Result AnnData ===
        # Use the first dataset (RNA) as the base
        adata = adata_rna.copy()
        adata.obsm['SpaMultiVAE'] = final_latent
    
    # === Parse Dataset Info ===
    dataset_name, subset_name = parse_dataset_info(args)
    print(f"Detected dataset: {dataset_name}, subset: {subset_name}")

    # === Plot Save Path ===
    plot_base_dir = "Results/plot"
    method_name = args.method if args.method else "SpaMultiVAE"
    plot_dir = os.path.join(plot_base_dir, method_name, dataset_name, subset_name)
    os.makedirs(plot_dir, exist_ok=True)
    print(f"Plot images will be saved to: {plot_dir}")

    # === Clustering and Visualization ===
    # tools = ['mclust', 'louvain', 'leiden', 'kmeans']
    tools = ['mclust']
    for tool in tools:
        adata = universal_clustering(
            adata,
            n_clusters=args.cluster_nums,
            used_obsm='SpaMultiVAE',
            method=tool,
            key=tool,
            use_pca=False
        )
        
        # Flip spatial coordinates for visualization if available
        if 'spatial' in adata.obsm.keys():
            adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

        fig, ax_list = plt.subplots(1, 2, figsize=(7, 3))
        
        # Generate UMAP for visualization
        sc.pp.neighbors(adata, use_rep='SpaMultiVAE', n_neighbors=30)
        sc.tl.umap(adata)

        # Plot UMAP and spatial
        sc.pl.umap(adata, color=tool, ax=ax_list[0], title=f'{method_name}-{tool}', s=20, show=False)
        if 'spatial' in adata.obsm.keys():
            sc.pl.embedding(adata, basis='spatial', color=tool, ax=ax_list[1], title=f'{method_name}-{tool}', s=20, show=False)
        else:
            # If no spatial coordinates, plot UMAP again
            sc.pl.umap(adata, color=tool, ax=ax_list[1], title=f'{method_name}-{tool} (no spatial)', s=20, show=False)

        plt.tight_layout(w_pad=0.3)
        plt.savefig(
            os.path.join(plot_dir, f'clustering_{tool}_umap_spatial.png'),
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()

    # === Save AnnData ===
    save_dir = os.path.dirname(args.save_path)
    os.makedirs(save_dir, exist_ok=True)
    adata.write(args.save_path)
    print(adata)
    print('Saving results to...', args.save_path)


if __name__ == "__main__":
    # Set environment variables for R and threading
    os.environ['R_HOME'] = '/home/zhenghong/miniconda3/envs/smobench/lib/R'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    print("Starting SpaMultiVAE integration...")
    parser = argparse.ArgumentParser(description='Run SpaMultiVAE integration')
    parser.add_argument('--data_type', type=str, default='10x', help='Data type, e.g. 10x, SPOTS, MISAR, simulation')
    parser.add_argument('--RNA_path', type=str, required=True, help='Path to RNA adata')
    parser.add_argument('--ADT_path', type=str, default='', help='Path to ADT adata')
    parser.add_argument('--ATAC_path', type=str, default='', help='Path to ATAC adata')
    parser.add_argument('--save_path', type=str, required=True, help='Path to save integrated adata')
    parser.add_argument('--seed', type=int, default=2024, help='Random seed')
    parser.add_argument('--device', type=str, default='cuda', help='Device to use, e.g. cuda or cpu')

    parser.add_argument('--method', type=str, default='SpaMultiVAE', help='Method name for plotting')
    parser.add_argument('--dataset', type=str, default='', help='Dataset name, e.g. Human_Lymph_Nodes/A1. If not provided, auto-extracted from paths.')
    parser.add_argument('--enhance_resolution', action='store_true', help='Whether to enhance spatial resolution')

    parser.add_argument('--cluster_nums', type=int, help='Number of clusters')

    args = parser.parse_args()
    main(args)