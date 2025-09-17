#!/usr/bin/env python3
"""
Convert swruan (SpaMosaic) results to SMOBench standard format
This script converts SpaMosaic embeddings and applies 4 clustering methods
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import re
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import sys

# Add project root to path
project_root = os.path.abspath(os.path.dirname(__file__))
sys.path.append(project_root)

from Utils.SMOBench_clustering import universal_clustering


# Dataset mapping for all available datasets (excluding Fusion)
DATASET_MAPPING = {
    'HLN': {
        'A1': 'A1',
        'D1': 'D1',
    },
    'HT': {
        'slice1': 'S1',
        'slice2': 'S2', 
        'slice3': 'S3',
    },
    'MISAR': {
        'E11': 'E11',
        'E13': 'E13',
        'E15': 'E15',
        'E18': 'E18',
    },
    'MISARS2': {
        'E11': 'E11',
        'E13': 'E13',
        'E15': 'E15',
        'E18': 'E18',
    },
    'Mouse_Thymus': {
        'Dataset3': 'Mouse_Thymus1',
        'Dataset4': 'Mouse_Thymus2', 
        'Dataset5': 'Mouse_Thymus3',
        'Dataset6': 'Mouse_Thymus4',
    },
    'Mouse_Brain': {
        'Dataset7': 'Mouse_Brain_ATAC',
        'Dataset8': 'Mouse_Brain_H3K27ac',
        'Dataset9': 'Mouse_Brain_H3K27me3',
        'Dataset10': 'Mouse_Brain_H3K4me3',
    }
}

def parse_running_time(time_file_path):
    """Parse running time from txt file"""
    try:
        with open(time_file_path, 'r') as f:
            content = f.read()
        
        # Extract running time in seconds
        time_match = re.search(r'Running Time:\s*([\d.]+)\s*seconds', content)
        if time_match:
            return float(time_match.group(1))
        else:
            print(f"Warning: Could not parse running time from {time_file_path}")
            return 0.0
    except Exception as e:
        print(f"Warning: Error reading {time_file_path}: {e}")
        return 0.0


def load_original_adata(dataset_name, subset_name, data_type='woGT'):
    """Load original AnnData to get metadata and coordinates"""
    
    # Determine modality and data_type based on dataset
    if dataset_name == 'HLN':
        modality = 'RNA_ADT'
        data_type = 'withGT'
        original_dataset = 'Human_Lymph_Nodes'
    elif dataset_name == 'HT':
        modality = 'RNA_ADT'
        data_type = 'withGT'
        original_dataset = 'Human_Tonsils'
    elif dataset_name == 'MISAR':
        modality = 'RNA_ATAC'
        data_type = 'withGT'
        original_dataset = 'Mouse_Embryos_S1'
    elif dataset_name == 'MISARS2':
        modality = 'RNA_ATAC'
        data_type = 'withGT'
        original_dataset = 'Mouse_Embryos_S2'
    elif 'Brain' in dataset_name:
        modality = 'RNA_ATAC'
        data_type = 'woGT'
        original_dataset = dataset_name
    elif 'Thymus' in dataset_name:
        modality = 'RNA_ADT'
        data_type = 'woGT'
        original_dataset = dataset_name
    else:
        print(f"Warning: Unknown dataset type: {dataset_name}")
        return None
    
    # Construct original data path
    rna_path = f"Dataset/{data_type}/{modality}/{original_dataset}/{subset_name}/adata_RNA.h5ad"
    
    if not os.path.exists(rna_path):
        print(f"Warning: Original RNA data not found at {rna_path}")
        return None
    
    try:
        adata = sc.read_h5ad(rna_path)
        print(f"Loaded original data: {adata.n_obs} cells, {adata.n_vars} genes")
        return adata
    except Exception as e:
        print(f"Error loading original data: {e}")
        return None


def convert_swruan_to_adata(swruan_dir, output_dir, method_name="spamosaic"):
    """Convert swruan results to AnnData format"""
    
    results_processed = 0
    results_skipped = 0
    
    # Process each dataset category
    for dataset_category, dataset_map in DATASET_MAPPING.items():
        dataset_dir = os.path.join(swruan_dir, dataset_category)
        
        if not os.path.exists(dataset_dir):
            print(f"Warning: Dataset directory not found: {dataset_dir}")
            continue
            
        # Process each dataset
        for swruan_dataset, smobench_subset in dataset_map.items():
            swruan_dataset_dir = os.path.join(dataset_dir, swruan_dataset)
            
            if not os.path.exists(swruan_dataset_dir):
                print(f"Warning: swruan dataset directory not found: {swruan_dataset_dir}")
                continue
            
            print(f"\n=== Processing {dataset_category}/{swruan_dataset} -> {dataset_category}/{smobench_subset} ===")
            
            # File paths
            merged_emb_path = os.path.join(swruan_dataset_dir, "batch0_merged_embedding.npy")
            rna_emb_path = os.path.join(swruan_dataset_dir, "rna_batch0_embedding.npy")
            adt_emb_path = os.path.join(swruan_dataset_dir, "adt_batch0_embedding.npy")
            leiden_path = os.path.join(swruan_dataset_dir, "leiden_y_pred_label.npy")
            time_path = os.path.join(swruan_dataset_dir, "running_time_and_memory.txt")
            
            # Check required files
            required_files = [merged_emb_path, leiden_path, time_path]
            missing_files = [f for f in required_files if not os.path.exists(f)]
            
            if missing_files:
                print(f"Skipping {swruan_dataset}: Missing files {missing_files}")
                results_skipped += 1
                continue
            
            try:
                # Load embeddings
                merged_embedding = np.load(merged_emb_path)
                leiden_labels = np.load(leiden_path, allow_pickle=True)
                train_time = parse_running_time(time_path)
                
                print(f"Loaded embeddings: {merged_embedding.shape}")
                print(f"Loaded leiden labels: {leiden_labels.shape}, unique labels: {len(np.unique(leiden_labels))}")
                print(f"Training time: {train_time:.2f} seconds")
                
                # Load original data for metadata
                original_adata = load_original_adata(dataset_category, smobench_subset)
                
                if original_adata is None:
                    print(f"Warning: Could not load original data, creating minimal AnnData")
                    # Create minimal AnnData
                    n_cells = merged_embedding.shape[0]
                    adata = sc.AnnData(X=np.zeros((n_cells, 1)))
                    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
                    adata.var_names = ["dummy_gene"]
                else:
                    # Use original data as base
                    adata = original_adata.copy()
                    
                    # Check cell count consistency
                    if adata.n_obs != merged_embedding.shape[0]:
                        print(f"Warning: Cell count mismatch - original: {adata.n_obs}, embedding: {merged_embedding.shape[0]}")
                        # Take the intersection or handle as needed
                        min_cells = min(adata.n_obs, merged_embedding.shape[0])
                        adata = adata[:min_cells].copy()
                        merged_embedding = merged_embedding[:min_cells]
                        leiden_labels = leiden_labels[:min_cells]
                        print(f"Adjusted to {min_cells} cells")
                
                # Store SpaMosaic results
                adata.obsm[method_name] = merged_embedding
                adata.obs['leiden'] = pd.Categorical(leiden_labels.astype(str))
                adata.uns['train_time'] = train_time
                
                # Store individual modality embeddings if available
                if os.path.exists(rna_emb_path):
                    rna_embedding = np.load(rna_emb_path)
                    if rna_embedding.shape[0] == adata.n_obs:
                        adata.obsm[f'{method_name}_rna'] = rna_embedding
                
                if os.path.exists(adt_emb_path):
                    adt_embedding = np.load(adt_emb_path)
                    if adt_embedding.shape[0] == adata.n_obs:
                        adata.obsm[f'{method_name}_adt'] = adt_embedding
                elif dataset_category in ['MISAR', 'MISARS2', 'Mouse_Brain']:
                    # For ATAC datasets, check for atac embedding
                    atac_emb_path = os.path.join(swruan_dataset_dir, "atac_batch0_embedding.npy")
                    if os.path.exists(atac_emb_path):
                        atac_embedding = np.load(atac_emb_path)
                        if atac_embedding.shape[0] == adata.n_obs:
                            adata.obsm[f'{method_name}_atac'] = atac_embedding
                
                # Clean embeddings for clustering (handle inf/nan values)
                embeddings_clean = adata.obsm[method_name].copy()
                if np.any(~np.isfinite(embeddings_clean)):
                    print("Warning: Found infinite or NaN values in embeddings. Cleaning...")
                    embeddings_clean[np.isinf(embeddings_clean)] = np.sign(embeddings_clean[np.isinf(embeddings_clean)]) * 1e10
                    embeddings_clean[np.isnan(embeddings_clean)] = 0
                    adata.obsm[method_name] = embeddings_clean
                    print(f"Cleaned embeddings shape: {embeddings_clean.shape}")
                
                # Generate UMAP for visualization
                print("Generating UMAP...")
                sc.pp.neighbors(adata, use_rep=method_name, n_neighbors=30)
                sc.tl.umap(adata)
                
                # Apply all 4 clustering methods to SpaMosaic embeddings  
                # Use standard SMOBench cluster numbers, not swruan's leiden result
                clustering_methods = ['mclust', 'leiden', 'louvain', 'kmeans']
                
                # Determine correct cluster numbers based on dataset
                cluster_nums_map = {
                    # withGT datasets
                    'HLN_A1': 10, 'HLN_D1': 11,
                    'HT_S1': 4, 'HT_S2': 5, 'HT_S3': 5,
                    'MISAR_E11': 8, 'MISAR_E13': 12, 'MISAR_E15': 12, 'MISAR_E18': 14,
                    'MISARS2_E11': 13, 'MISARS2_E13': 14, 'MISARS2_E15': 15, 'MISARS2_E18': 16,
                    # woGT datasets  
                    'Mouse_Thymus': 8, 'Mouse_Spleen': 5, 'Mouse_Brain': 18
                }
                
                # Generate dataset key for lookup
                if dataset_category == 'HLN':
                    dataset_key = f"HLN_{smobench_subset}"
                elif dataset_category == 'HT':
                    dataset_key = f"HT_{smobench_subset}"
                elif dataset_category == 'MISAR':
                    dataset_key = f"MISAR_{smobench_subset}"
                elif dataset_category == 'MISARS2':
                    dataset_key = f"MISARS2_{smobench_subset}"
                elif dataset_category == 'Mouse_Thymus':
                    dataset_key = 'Mouse_Thymus'
                elif dataset_category == 'Mouse_Brain':
                    dataset_key = 'Mouse_Brain'
                else:
                    dataset_key = 'Mouse_Spleen'  # Mouse_slpeen
                
                n_clusters = cluster_nums_map.get(dataset_key, 8)  # Default to 8 if not found
                print(f"Applying all clustering methods with {n_clusters} clusters (SMOBench standard)...")
                
                for method in clustering_methods:
                    try:
                        # Skip leiden if already exists and has correct number of clusters
                        if method == 'leiden' and 'leiden' in adata.obs.columns:
                            existing_clusters = len(adata.obs['leiden'].unique())
                            if existing_clusters == n_clusters:
                                print(f"✅ {method} clustering already exists with {existing_clusters} clusters")
                                continue
                        
                        adata = universal_clustering(
                            adata,
                            n_clusters=n_clusters,
                            used_obsm=method_name,
                            method=method,
                            key=method,
                            use_pca=False
                        )
                        print(f"✅ {method} clustering completed")
                    except Exception as e:
                        print(f"❌ {method} clustering failed: {e}")
                        # Create placeholder if clustering fails
                        adata.obs[method] = pd.Categorical(['0'] * adata.n_obs)
                
                # Create output directory based on dataset type - use t_spamosaic folder
                folder_name = f"t_{method_name}"
                
                if dataset_category == 'HLN':
                    output_dataset_dir = os.path.join(output_dir, folder_name, "HLN", smobench_subset)
                    output_filename = f"t_{method_name}_HLN_{smobench_subset}.h5ad"
                elif dataset_category == 'HT':
                    output_dataset_dir = os.path.join(output_dir, folder_name, "HT", smobench_subset)
                    output_filename = f"t_{method_name}_HT_{smobench_subset}.h5ad"
                elif dataset_category == 'MISAR':
                    output_dataset_dir = os.path.join(output_dir, folder_name, "MISAR_S1", smobench_subset)
                    output_filename = f"t_{method_name}_MISAR_S1_{smobench_subset}.h5ad"
                elif dataset_category == 'MISARS2':
                    output_dataset_dir = os.path.join(output_dir, folder_name, "MISAR_S2", smobench_subset)
                    output_filename = f"t_{method_name}_MISAR_S2_{smobench_subset}.h5ad"
                elif dataset_category == 'Mouse_Thymus':
                    output_dataset_dir = os.path.join(output_dir, folder_name, "Mouse_Thymus", smobench_subset.replace('Mouse_Thymus', 'Thymus'))
                    output_filename = f"t_{method_name}_MT_{smobench_subset.replace('Mouse_Thymus', 'Thymus')}.h5ad"
                elif dataset_category == 'Mouse_Brain':
                    brain_type = smobench_subset.replace('Mouse_Brain_', '')
                    output_dataset_dir = os.path.join(output_dir, folder_name, "Mouse_Brain", brain_type)
                    output_filename = f"t_{method_name}_MB_{brain_type}.h5ad"
                else:
                    output_dataset_dir = os.path.join(output_dir, folder_name, dataset_category, smobench_subset)
                    output_filename = f"t_{method_name}_{dataset_category}_{smobench_subset}.h5ad"
                
                os.makedirs(output_dataset_dir, exist_ok=True)
                
                # Save AnnData
                output_path = os.path.join(output_dataset_dir, output_filename)
                adata.write(output_path)
                print(f"✅ Saved to: {output_path}")
                
                # Generate visualization plots
                plot_dir = os.path.join("Results/plot", folder_name, dataset_category, smobench_subset)
                os.makedirs(plot_dir, exist_ok=True)
                
                clustering_methods = ['leiden', 'mclust', 'louvain', 'kmeans']
                available_methods = [m for m in clustering_methods if m in adata.obs.columns]
                
                for tool in available_methods:
                    try:
                        # Flip spatial coordinates for visualization if available
                        if 'spatial' in adata.obsm.keys():
                            adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

                        fig, ax_list = plt.subplots(1, 2, figsize=(7, 3))

                        # Plot UMAP and spatial
                        sc.pl.umap(adata, color=tool, ax=ax_list[0], title=f'{method_name}-{tool}', s=20, show=False)
                        if 'spatial' in adata.obsm.keys():
                            sc.pl.embedding(adata, basis='spatial', color=tool, ax=ax_list[1], 
                                          title=f'{method_name}-{tool}', s=20, show=False)
                        else:
                            # If no spatial coordinates, plot UMAP again
                            sc.pl.umap(adata, color=tool, ax=ax_list[1], 
                                     title=f'{method_name}-{tool} (no spatial)', s=20, show=False)

                        plt.tight_layout(w_pad=0.3)
                        plot_path = os.path.join(plot_dir, f'clustering_{tool}_umap_spatial.png')
                        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                        plt.close()
                        print(f"Plot saved: {plot_path}")
                        
                    except Exception as e:
                        print(f"Error generating plot for {tool}: {e}")
                
                results_processed += 1
                
            except Exception as e:
                print(f"Error processing {swruan_dataset}: {e}")
                results_skipped += 1
                continue
    
    print(f"\n=== Summary ===")
    print(f"Processed: {results_processed}")
    print(f"Skipped: {results_skipped}")
    print(f"Total: {results_processed + results_skipped}")


def main():
    parser = argparse.ArgumentParser(description='Convert swruan results to SMOBench format')
    parser.add_argument('--swruan_dir', type=str, default='swruan/results/embedding',
                       help='Path to swruan results directory')
    parser.add_argument('--output_dir', type=str, default='Results/adata',
                       help='Output directory for converted results')
    parser.add_argument('--method_name', type=str, default='spamosaic',
                       help='Method name to use in output')
    
    args = parser.parse_args()
    
    # Set environment variables for R and threading
    os.environ['R_HOME'] = '/home/zhenghong/miniconda3/envs/smobench/lib/R'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    
    print("Starting SpaMosaic (swruan) results conversion...")
    print(f"Input directory: {args.swruan_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Method name: {args.method_name}")
    
    convert_swruan_to_adata(args.swruan_dir, args.output_dir, args.method_name)
    print("SpaMosaic conversion completed!")


if __name__ == "__main__":
    main()