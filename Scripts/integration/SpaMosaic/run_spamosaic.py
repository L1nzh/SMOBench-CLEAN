#!/usr/bin/env python3
"""
Unified SpaMosaic Integration Script
Supports vertical, horizontal, and mosaic integration tasks
Based on SMOBench standards and SpaMosaic official tutorials
"""

import os
import sys
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import pynvml
import time
import warnings
from pathlib import Path
from os.path import join

# Add SMOBench utils to path
current_dir = Path(__file__).parent.absolute()
utils_path = current_dir.parent.parent.parent / "Utils"
sys.path.insert(0, str(utils_path))

# Import SMOBench utilities
from SMOBench_clustering import universal_clustering

# SpaMosaic imports
try:
    from spamosaic.framework import SpaMosaic
    from spamosaic.preprocessing2 import RNA_preprocess, ADT_preprocess, Epigenome_preprocess
    import spamosaic.utils as utls
except ImportError as e:
    print(f"Error importing SpaMosaic: {e}")
    print("Please ensure SpaMosaic is properly installed")
    sys.exit(1)

# Suppress warnings
warnings.filterwarnings('ignore')
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Unified SpaMosaic Integration')
    
    # Data arguments
    parser.add_argument('--dataset', type=str, required=True,
                        help='Dataset name (e.g., HLN, MISAR, Mouse_Brain)')
    parser.add_argument('--samples', type=str, nargs='+', required=True,
                        help='Sample names (e.g., A1 D1 for HLN)')
    parser.add_argument('--data_dir', type=str, required=True,
                        help='Base data directory path')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory path')
    
    # Integration type and modalities
    parser.add_argument('--task', type=str, choices=['vertical', 'horizontal', 'mosaic'], 
                        required=True, help='Integration task type')
    parser.add_argument('--modalities', type=str, nargs='+', default=['rna', 'adt'],
                        choices=['rna', 'adt', 'atac'],
                        help='Modalities to integrate (default: rna adt)')
    
    # Model parameters
    parser.add_argument('--gpu', type=int, default=0, help='GPU device ID')
    parser.add_argument('--lr', type=float, default=0.01, help='Learning rate')
    parser.add_argument('--epochs', type=int, default=100, help='Training epochs')
    parser.add_argument('--intra_knn', type=int, default=2, help='Intra-KNN neighbors')
    parser.add_argument('--inter_knn', type=int, default=2, help='Inter-KNN neighbors')
    parser.add_argument('--w_g', type=float, default=0.8, help='Graph weight')
    parser.add_argument('--seed', type=int, default=2024, help='Random seed')
    
    # Preprocessing parameters
    parser.add_argument('--n_hvg', type=int, default=5000, help='Number of highly variable genes')
    parser.add_argument('--n_peak', type=int, default=50000, help='Number of peaks for ATAC')
    parser.add_argument('--batch_corr', action='store_true', default=True,
                        help='Apply batch correction')
    
    # Clustering parameters
    parser.add_argument('--clustering_methods', type=str, nargs='+',
                        default=['leiden', 'louvain', 'kmeans', 'mclust'],
                        help='Clustering methods to use')
    parser.add_argument('--n_clusters', type=int, default=None,
                        help='Number of clusters (auto-detect if not specified)')
    
    return parser.parse_args()

def load_data(data_dir, dataset, samples, modalities):
    """Load data based on dataset and samples"""
    input_dict = {mod: [] for mod in modalities}
    
    print(f"Loading data for dataset: {dataset}")
    print(f"Samples: {samples}")
    print(f"Modalities: {modalities}")
    
    for sample in samples:
        sample_dir = join(data_dir, sample)
        print(f"Loading from: {sample_dir}")
        
        for modality in modalities:
            if modality == 'rna':
                file_pattern = 'adata_RNA.h5ad'
            elif modality == 'adt':
                file_pattern = 'adata_ADT.h5ad'
            elif modality == 'atac':
                file_pattern = 'adata_ATAC.h5ad'
            
            file_path = join(sample_dir, file_pattern)
            
            if os.path.exists(file_path):
                print(f"  Loading {modality} from {file_path}")
                adata = sc.read_h5ad(file_path)
                # Ensure unique observation names
                adata.obs_names_make_unique()
                adata.var_names_make_unique()
                input_dict[modality].append(adata)
            else:
                print(f"  Warning: {file_path} not found")
                input_dict[modality].append(None)
    
    # Remove empty modalities
    input_dict = {k: v for k, v in input_dict.items() if any(x is not None for x in v)}
    
    return input_dict

def align_features(input_dict):
    """Align features across samples for each modality"""
    print("Aligning features across samples...")
    
    for modality, adata_list in input_dict.items():
        if not adata_list or len(adata_list) < 2:
            continue
            
        valid_adatas = [ad for ad in adata_list if ad is not None]
        if len(valid_adatas) < 2:
            continue
            
        if modality in ['rna']:
            # Gene alignment for RNA
            print(f"  Aligning genes for {modality}")
            common_genes = set(valid_adatas[0].var_names)
            for adata in valid_adatas[1:]:
                common_genes.intersection_update(adata.var_names)
            common_genes = sorted(common_genes)
            
            if not common_genes:
                raise ValueError(f"No common genes found for {modality}")
                
            print(f"    Found {len(common_genes)} common genes")
            for i, adata in enumerate(adata_list):
                if adata is not None:
                    input_dict[modality][i] = adata[:, common_genes]
                    
        elif modality in ['atac']:
            # Peak alignment for ATAC
            print(f"  Aligning peaks for {modality}")
            from Methods.SpaMosaic.utils import peak_sets_alignment
            input_dict[modality] = peak_sets_alignment(valid_adatas)

    return input_dict

def preprocess_data(input_dict, args):
    """Preprocess data according to modality"""
    print("Preprocessing data...")
    input_key = 'dimred_bc'
    
    for modality, adata_list in input_dict.items():
        if modality == 'rna':
            print("  Preprocessing RNA data")
            RNA_preprocess(adata_list, 
                         batch_corr=args.batch_corr,
                         favor='scanpy',
                         n_hvg=args.n_hvg,
                         batch_key='src',
                         key=input_key)
        elif modality == 'adt':
            print("  Preprocessing ADT data")
            ADT_preprocess(adata_list,
                         batch_corr=args.batch_corr,
                         batch_key='src',
                         key=input_key)
        elif modality == 'atac':
            print("  Preprocessing ATAC data")
            Epigenome_preprocess(adata_list,
                               batch_corr=args.batch_corr,
                               n_peak=args.n_peak,
                               batch_key='src',
                               key=input_key)
    
    return input_dict, input_key

def setup_gpu_monitoring(gpu_id):
    """Setup GPU monitoring"""
    try:
        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(gpu_id)
        info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        return handle, info.used / (1024 * 1024)
    except:
        print("Warning: GPU monitoring not available")
        return None, 0

def train_spamosaic(input_dict, input_key, args):
    """Train SpaMosaic model"""
    print("Training SpaMosaic model...")
    
    # Setup GPU monitoring
    gpu_handle, gpu_memory_before = setup_gpu_monitoring(args.gpu)
    
    # Create model
    device = f'cuda:{args.gpu}' if args.gpu >= 0 else 'cpu'
    model = SpaMosaic(
        modBatch_dict=input_dict,
        input_key=input_key,
        batch_key='src',
        intra_knn=args.intra_knn,
        inter_knn=args.inter_knn,
        w_g=args.w_g,
        seed=args.seed,
        device=device
    )
    
    # Training
    start_time = time.time()
    model.train(net='wlgcn', lr=args.lr, T=0.01, n_epochs=args.epochs)
    
    # Inference
    ad_embs = model.infer_emb(input_dict, emb_key='emb', final_latent_key='merged_emb')
    ad_mosaic = sc.concat(ad_embs)
    ad_mosaic = utls.get_umap(ad_mosaic, use_reps=['merged_emb'])
    
    # Calculate metrics
    training_time = time.time() - start_time
    gpu_memory_after = 0
    if gpu_handle:
        try:
            info = pynvml.nvmlDeviceGetMemoryInfo(gpu_handle)
            gpu_memory_after = info.used / (1024 * 1024)
        except:
            pass
    
    return model, ad_embs, ad_mosaic, training_time, gpu_memory_before, gpu_memory_after

def perform_clustering(ad_mosaic, args):
    """Perform clustering using SMOBench universal clustering"""
    print("Performing clustering...")
    
    clustering_results = {}
    
    for method in args.clustering_methods:
        print(f"  Clustering with {method}")
        try:
            # Use SMOBench universal clustering
            result_adata = universal_clustering(
                adata=ad_mosaic.copy(),
                embedding_key='merged_emb',
                method=method,
                n_clusters=args.n_clusters,
                random_state=args.seed
            )
            
            # Get cluster labels
            if method == 'mclust':
                cluster_key = 'mclust'
            else:
                cluster_key = method
                
            if cluster_key in result_adata.obs.columns:
                clustering_results[method] = result_adata.obs[cluster_key].values
            else:
                print(f"    Warning: {cluster_key} not found in results")
                
        except Exception as e:
            print(f"    Error with {method}: {e}")
            continue
    
    return clustering_results

def save_results(model, ad_embs, ad_mosaic, clustering_results, 
                training_time, gpu_memory_before, gpu_memory_after, 
                input_dict, args):
    """Save all results"""
    print("Saving results...")
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save individual embeddings
    embed_dir = join(args.output_dir, 'embeddings')
    os.makedirs(embed_dir, exist_ok=True)
    
    for k in model.mod_list:
        for bi in range(model.n_batches):
            if input_dict[k][bi] is not None:
                emb = input_dict[k][bi].obsm['emb']
                np.save(join(embed_dir, f'{k}_batch{bi}_embedding.npy'), emb)
    
    # Save merged embeddings
    for bi, ad in enumerate(ad_embs):
        np.save(join(embed_dir, f'batch{bi}_merged_embedding.npy'), ad.obsm['merged_emb'])
    
    # Save clustering results
    cluster_dir = join(args.output_dir, 'clustering')
    os.makedirs(cluster_dir, exist_ok=True)
    
    for method, labels in clustering_results.items():
        np.save(join(cluster_dir, f'{method}_labels.npy'), labels)
    
    # Save UMAP
    if 'X_umap' in ad_mosaic.obsm:
        np.save(join(args.output_dir, 'umap.npy'), ad_mosaic.obsm['X_umap'])
    
    # Save integrated AnnData
    ad_mosaic.write_h5ad(join(args.output_dir, f'SpaMosaic_{args.dataset}_integrated.h5ad'))
    
    # Save metrics
    metrics = {
        'training_time': training_time,
        'gpu_memory_before': gpu_memory_before,
        'gpu_memory_after': gpu_memory_after,
        'dataset': args.dataset,
        'task': args.task,
        'modalities': args.modalities,
        'samples': args.samples,
        'n_epochs': args.epochs,
        'lr': args.lr,
        'seed': args.seed
    }
    
    with open(join(args.output_dir, 'metrics.txt'), 'w') as f:
        for key, value in metrics.items():
            f.write(f"{key}: {value}\n")
    
    print(f"Results saved to: {args.output_dir}")

def main():
    """Main execution function"""
    args = parse_arguments()
    
    print("="*60)
    print(f"SpaMosaic Integration: {args.dataset}")
    print(f"Task: {args.task}")
    print(f"Modalities: {args.modalities}")
    print(f"Samples: {args.samples}")
    print("="*60)
    
    # Set random seed
    np.random.seed(args.seed)
    
    # Load data
    input_dict = load_data(args.data_dir, args.dataset, args.samples, args.modalities)
    
    if not input_dict:
        print("Error: No data loaded")
        return
    
    # Align features if needed
    if args.task in ['horizontal', 'mosaic'] and len(args.samples) > 1:
        input_dict = align_features(input_dict)
    
    # Add sample identifiers for multi-sample integration
    if len(args.samples) > 1:
        for modality, adata_list in input_dict.items():
            for i, adata in enumerate(adata_list):
                if adata is not None:
                    adata.obs_names = f"s{i}_" + adata.obs_names
    
    # Preprocess data
    input_dict, input_key = preprocess_data(input_dict, args)
    
    # Train model
    model, ad_embs, ad_mosaic, training_time, gpu_before, gpu_after = train_spamosaic(
        input_dict, input_key, args)
    
    # Perform clustering
    clustering_results = perform_clustering(ad_mosaic, args)
    
    # Save results
    save_results(model, ad_embs, ad_mosaic, clustering_results,
                training_time, gpu_before, gpu_after, input_dict, args)
    
    print("SpaMosaic integration completed successfully!")

if __name__ == '__main__':
    main()