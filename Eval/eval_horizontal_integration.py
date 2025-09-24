#!/usr/bin/env python3
"""
Horizontal Integration Evaluation Script for SMOBench
Evaluates horizontal integration methods using three dimensions: SC (Spatial Coherence), BVC (Biological Conservation), BER (Batch Effect Removal)
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# Add current directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

# Import evaluation functions
from src.demo import eval_horizontal_integration, calculate_horizontal_dataset_summary, save_evaluation_results

# Dataset classification (same as vertical integration)
WITHGT_DATASETS = ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2']
WOGT_DATASETS = ['Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain']

# Integration methods (same as vertical integration)
METHODS = ['CANDIES', 'COSMOS', 'PRAGA', 'PRESENT', 'SpaMV', 'SpaMosaic', 'SpatialGlue', 'SpaMultiVAE']

# Method-dataset compatibility for horizontal integration
HORIZONTAL_METHOD_DATASET_COMPATIBILITY = {
    # SpatialGlue: 垂直7+水平7 - 支持所有数据集
    'SpatialGlue': ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2', 'Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain'],
    
    # SpaMV: 垂直7+水平6，水平没有mouse brain
    'SpaMV': ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2', 'Mouse_Thymus', 'Mouse_Spleen'],
    
    # CANDIES: 垂直7+水平5，水平没有mouse thymus/mouse brain
    'CANDIES': ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2', 'Mouse_Spleen'],
    
    # SpaMosaic: 垂直7+水平7 - 支持所有数据集
    'SpaMosaic': ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2', 'Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain'],
    
    # PRAGA: 垂直7+水平3，水平没有mouse thymus/mouse brain/mouse embryo
    'PRAGA': ['HLN', 'HT', 'Mouse_Spleen'],
    
    # PRESENT: 垂直7+水平7 - 支持所有数据集
    'PRESENT': ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2', 'Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain'],
    
    # COSMOS: 垂直7+水平7 - 支持所有数据集
    'COSMOS': ['HLN', 'HT', 'MISAR_S1', 'MISAR_S2', 'Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain'],
    
    # SpaMultiVAE: 垂直3+水平3，垂直和水平只支持RNA_ADT
    'SpaMultiVAE': ['HLN', 'HT', 'Mouse_Thymus', 'Mouse_Spleen']
}

# Clustering methods
CLUSTERING_METHODS = ['leiden', 'louvain', 'kmeans', 'mclust']

def get_ground_truth_labels(dataset_name):
    """
    Load ground truth labels from original dataset for horizontal integration
    
    Parameters:
    -----------
    dataset_name : str
        Name of the dataset (e.g., 'HLN', 'HT')
        
    Returns:
    --------
    np.ndarray or None : Ground truth labels if available
    """
    if dataset_name not in WITHGT_DATASETS:
        return None
    
    # Map dataset names to original paths
    dataset_mapping = {
        'HLN': 'Human_Lymph_Nodes',
        'HT': 'Human_Tonsils', 
        'MISAR_S1': 'Mouse_Embryos_S1',
        'MISAR_S2': 'Mouse_Embryos_S2'
    }
    
    # Determine data type
    if dataset_name in ['HLN', 'HT']:
        data_type = 'RNA_ADT'
    else:
        data_type = 'RNA_ATAC'
    
    # For horizontal integration, we need to load all slices and concatenate
    original_dataset = dataset_mapping.get(dataset_name, dataset_name)
    dataset_path = f"/home/zhenghong/SMOBench-CLEAN/Dataset/withGT/{data_type}/{original_dataset}"
    
    if not os.path.exists(dataset_path):
        print(f"Warning: Dataset path not found: {dataset_path}")
        return None
    
    # Load all slices for the dataset
    all_gt_labels = []
    slice_dirs = [d for d in os.listdir(dataset_path) if os.path.isdir(os.path.join(dataset_path, d))]
    
    for slice_name in sorted(slice_dirs):
        gt_path = os.path.join(dataset_path, slice_name, "adata_RNA.h5ad")
        try:
            if os.path.exists(gt_path):
                adata_gt = sc.read_h5ad(gt_path)
                if 'Spatial_Label' in adata_gt.obs.columns:
                    all_gt_labels.extend(adata_gt.obs['Spatial_Label'].values)
                else:
                    print(f"Warning: No 'Spatial_Label' found in {gt_path}")
        except Exception as e:
            print(f"Warning: Could not load GT from {gt_path}: {e}")
    
    return np.array(all_gt_labels) if all_gt_labels else None

def get_batch_labels_for_horizontal(adata):
    """
    Extract batch labels for horizontal integration
    In horizontal integration, batch labels typically indicate which original slice/sample each cell came from
    
    Parameters:
    -----------
    adata : AnnData
        Integrated AnnData object
        
    Returns:
    --------
    np.ndarray or None : Batch labels for BER metrics
    """
    # Common batch label columns
    batch_columns = ['batch', 'slice', 'sample', 'orig.ident', 'batch_id']
    
    for col in batch_columns:
        if col in adata.obs.columns:
            print(f"Using '{col}' as batch labels for BER metrics")
            batch_values = adata.obs[col].values
            # Convert pandas Categorical to string array if needed
            if hasattr(batch_values, 'categories'):
                batch_values = batch_values.astype(str)
            return batch_values
    
    # If no explicit batch column, try to infer from obs_names or other metadata
    if 'slice_id' in adata.obs.columns:
        print("Using 'slice_id' as batch labels for BER metrics")
        slice_values = adata.obs['slice_id'].values
        if hasattr(slice_values, 'categories'):
            slice_values = slice_values.astype(str)
        return slice_values
    
    # Fallback: create artificial batch labels based on index patterns
    # This is a heuristic - real horizontal integration should have proper batch labels
    print("Warning: No batch labels found, creating artificial batch labels based on cell indices")
    n_cells = adata.n_obs
    n_batches = min(4, max(2, n_cells // 1000))  # Assume 2-4 batches
    batch_size = n_cells // n_batches
    
    batch_labels = []
    for i in range(n_batches):
        start_idx = i * batch_size
        end_idx = (i + 1) * batch_size if i < n_batches - 1 else n_cells
        batch_labels.extend([f'batch_{i}'] * (end_idx - start_idx))
    
    return np.array(batch_labels[:n_cells])

def evaluate_horizontal_integration_methods(
    input_base_dir="/home/zhenghong/SMOBench-CLEAN/Results/adata/horizontal_integration",
    output_base_dir="/home/zhenghong/SMOBench-CLEAN/Results/evaluation/horizontal_integration",
    test_mode=False
):
    """
    Evaluate all horizontal integration methods with three dimensions: SC, BVC, BER
    
    Parameters:
    -----------
    input_base_dir : str
        Base directory containing integration results
    output_base_dir : str
        Base directory for evaluation outputs
    test_mode : bool
        If True, process only a few files for testing
    """
    
    print("=== SMOBench Horizontal Integration Evaluation ===")
    print("Three-dimensional evaluation: SC (Spatial Coherence) + BVC (Biological Conservation) + BER (Batch Effect Removal)")
    print()
    
    os.makedirs(output_base_dir, exist_ok=True)
    
    all_results = []
    files_processed = 0
    max_test_files = 5 if test_mode else float('inf')
    
    # Process each method
    for method_name in METHODS:
        method_dir = os.path.join(input_base_dir, method_name)
        if not os.path.exists(method_dir):
            print(f"Method directory not found: {method_dir}")
            continue
        
        print(f"\n--- Processing Method: {method_name} ---")
        
        # Create output directory for this method
        method_output_dir = os.path.join(output_base_dir, method_name)
        os.makedirs(method_output_dir, exist_ok=True)
        
        # Process each dataset
        for dataset_name in os.listdir(method_dir):
            dataset_path = os.path.join(method_dir, dataset_name)
            if not os.path.isdir(dataset_path):
                continue
            
            # Check method-dataset compatibility for horizontal integration
            if method_name in HORIZONTAL_METHOD_DATASET_COMPATIBILITY:
                if dataset_name not in HORIZONTAL_METHOD_DATASET_COMPATIBILITY[method_name]:
                    print(f"  Skipping {dataset_name} (not supported by {method_name} in horizontal integration)")
                    continue
            
            print(f"\n  Dataset: {dataset_name}")
            
            # Load horizontal integration result file
            h5ad_file = f"{method_name}_{dataset_name}_horizontal.h5ad"
            h5ad_path = os.path.join(dataset_path, h5ad_file)
            
            if not os.path.exists(h5ad_path):
                print(f"    File not found: {h5ad_path}")
                continue
            
            if test_mode and files_processed >= max_test_files:
                print(f"    Test mode: stopping after {max_test_files} files")
                break
            
            try:
                print(f"    Loading: {h5ad_path}")
                adata = sc.read_h5ad(h5ad_path)
                
                # Get ground truth labels
                y_GT = get_ground_truth_labels(dataset_name)
                has_gt = y_GT is not None
                
                # Get batch labels for BER metrics
                batch_labels = get_batch_labels_for_horizontal(adata)
                
                print(f"    Data shape: {adata.shape}, GT available: {has_gt}, Batches: {len(np.unique(batch_labels)) if batch_labels is not None else 'None'}")
                
                # Process each clustering method
                dataset_results = []
                for clustering_method in CLUSTERING_METHODS:
                    # For horizontal integration, clustering columns are named directly (leiden, louvain, etc.)
                    cluster_key = clustering_method
                    
                    if cluster_key not in adata.obs.columns:
                        print(f"    Warning: {cluster_key} not found in adata.obs")
                        continue
                    
                    # Get required data
                    try:
                        # Get embeddings
                        # First try method-specific embedding (e.g., 'CANDIES', 'SpatialGlue', etc.)
                        method_embedding_key = method_name
                        if method_embedding_key in adata.obsm:
                            embeddings = adata.obsm[method_embedding_key]
                            print(f"    Using embeddings from: {method_embedding_key}")
                        else:
                            # Fallback to common embedding keys
                            embedding_keys = ['X_integrated', 'X_emb', 'X_pca', 'X_umap']
                            embeddings = None
                            for key in embedding_keys:
                                if key in adata.obsm:
                                    embeddings = adata.obsm[key]
                                    print(f"    Using embeddings from: {key}")
                                    break
                        
                        if embeddings is None:
                            print(f"    Warning: No embeddings found for {dataset_name}")
                            continue
                        
                        # Get clustering results
                        y_pred = adata.obs[cluster_key].values
                        
                        # Convert string labels to numeric
                        if y_pred.dtype == 'object':
                            unique_labels = sorted(set(y_pred))
                            label_map = {label: i for i, label in enumerate(unique_labels)}
                            y_pred = np.array([label_map[label] for label in y_pred])
                        
                        # Get spatial coordinates
                        if 'spatial' in adata.obsm:
                            spatial_coords = adata.obsm['spatial']
                        elif 'X_spatial' in adata.obsm:
                            spatial_coords = adata.obsm['X_spatial']
                        else:
                            print(f"    Warning: No spatial coordinates found, using first 2 embedding dimensions")
                            spatial_coords = embeddings[:, :2]
                        
                        # Create adjacency matrix (simplified)
                        n_cells = len(y_pred)
                        adj_matrix = sparse.eye(n_cells, format='csr')
                        
                        # Align ground truth if available
                        if y_GT is not None and len(y_GT) != len(y_pred):
                            min_len = min(len(y_GT), len(y_pred))
                            print(f"    Warning: GT length mismatch. Truncating to {min_len}")
                            y_GT = y_GT[:min_len]
                            y_pred = y_pred[:min_len]
                            embeddings = embeddings[:min_len]
                            spatial_coords = spatial_coords[:min_len]
                            batch_labels = batch_labels[:min_len] if batch_labels is not None else None
                            adj_matrix = adj_matrix[:min_len, :min_len]
                        
                        # Evaluate with three dimensions
                        print(f"    Evaluating {clustering_method} clustering...")
                        metrics = eval_horizontal_integration(
                            embeddings=embeddings,
                            adj_matrix=adj_matrix,
                            y_pred=y_pred,
                            y_GT=y_GT,
                            spatial_coords=spatial_coords,
                            batch_labels=batch_labels,
                            method_name=method_name,
                            dataset_name=dataset_name,
                            slice_name="horizontal",  # For horizontal integration
                            clustering_method=clustering_method
                        )
                        
                        dataset_results.append(metrics)
                        all_results.append(metrics)
                        
                        # Save individual results
                        save_evaluation_results(
                            metrics_dict=metrics,
                            output_dir=method_output_dir,
                            method_name=method_name,
                            dataset_name=dataset_name,
                            slice_name="horizontal",
                            clustering_method=clustering_method,
                            has_gt=has_gt
                        )
                        
                        print(f"      {clustering_method}: SC={metrics['SC_Score']:.3f}, BVC={metrics['BVC_Score']:.3f}, BER={metrics['BER_Score']:.3f}, Final={metrics['Final_Score']:.3f}")
                        
                    except Exception as e:
                        print(f"    Error processing {clustering_method}: {e}")
                        continue
                
                # Note: Dataset summaries will be generated by generate_final_results.py
                
                files_processed += 1
                
            except Exception as e:
                print(f"    Error processing {h5ad_path}: {e}")
                continue
        
        if test_mode and files_processed >= max_test_files:
            break
    
    print(f"\n--- Horizontal Integration Evaluation Complete ---")
    print(f"Total results processed: {len(all_results)}")
    print(f"Results saved to: {output_base_dir}")
    
    return all_results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate horizontal integration methods with SC + BVC + BER')
    parser.add_argument('--test', action='store_true', help='Run in test mode (process only a few files)')
    parser.add_argument('--input', type=str, 
                       default='/home/zhenghong/SMOBench-CLEAN/Results/adata/horizontal_integration',
                       help='Input directory containing integration results')
    parser.add_argument('--output', type=str,
                       default='/home/zhenghong/SMOBench-CLEAN/Results/evaluation/horizontal_integration', 
                       help='Output directory for evaluation results')
    
    args = parser.parse_args()
    
    results = evaluate_horizontal_integration_methods(
        input_base_dir=args.input,
        output_base_dir=args.output,
        test_mode=args.test
    )
    
    print(f"\nEvaluation complete! Processed {len(results)} method-dataset-clustering combinations.")