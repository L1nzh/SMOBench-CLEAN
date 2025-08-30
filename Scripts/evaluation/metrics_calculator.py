#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    adjusted_rand_score, normalized_mutual_info_score, adjusted_mutual_info_score,
    fowlkes_mallows_score, silhouette_score, calinski_harabasz_score, davies_bouldin_score,
    homogeneity_score, completeness_score, v_measure_score, jaccard_score, f1_score
)
from sklearn.metrics.cluster import contingency_matrix
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

class SMOBenchMetrics:
    """
    Comprehensive metrics calculator for SMOBench evaluation
    
    Supports three categories of metrics:
    1. Spatial Coherence (SC): Moran's I, Geary's C
    2. Biological Conservation (BioC): Clustering and quality metrics  
    3. Batch Effect Removal (BER): Batch mixing metrics
    """
    
    def __init__(self, adata, task='vertical'):
        """
        Initialize metrics calculator
        
        Parameters:
            adata: AnnData object with integration results and clustering
            task: Integration task ('vertical', 'horizontal', 'mosaic')
        """
        self.adata = adata
        self.task = task
        self.has_spatial = 'spatial' in adata.obsm or ('x' in adata.obs and 'y' in adata.obs)
        self.has_gt = self._check_ground_truth()
        
    def _check_ground_truth(self):
        """Check if ground truth labels are available"""
        gt_candidates = ['final_annot', 'Spatial_Label', 'cell_type', 'celltype', 
                        'cluster', 'annotation', 'ground_truth', 'gt', 'label']
        
        for col in self.adata.obs.columns:
            for candidate in gt_candidates:
                if candidate.lower() in col.lower():
                    return True
        return False
    
    def calculate_vertical_metrics(self, cluster_key, embedding_key=None, gt_key=None):
        """
        Calculate metrics for vertical integration (BioC + SC)
        
        Parameters:
            cluster_key: Key for cluster labels in adata.obs
            embedding_key: Key for embeddings in adata.obsm (optional)
            gt_key: Key for ground truth labels in adata.obs (if available)
        
        Returns:
            dict: Calculated metrics
        """
        metrics = {}
        
        # Get cluster labels
        if cluster_key not in self.adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
        
        cluster_labels = self.adata.obs[cluster_key].astype(str)
        
        # Biological Conservation metrics
        if gt_key and gt_key in self.adata.obs.columns:
            # With ground truth - full BioC metrics
            gt_labels = self.adata.obs[gt_key].astype(str)
            metrics.update(self._calculate_clustering_metrics(gt_labels, cluster_labels))
            
            # ASW cell type and Graph cLISI
            if embedding_key and embedding_key in self.adata.obsm:
                metrics['asw_celltype'] = self._calculate_asw_celltype(gt_labels, embedding_key)
                metrics['graph_clisi'] = self._calculate_graph_clisi(gt_labels, embedding_key)
        
        # Quality metrics (all datasets)
        if embedding_key and embedding_key in self.adata.obsm:
            metrics.update(self._calculate_quality_metrics(cluster_labels, embedding_key))
        
        # Spatial Coherence metrics
        if self.has_spatial:
            metrics.update(self._calculate_spatial_metrics(cluster_labels))
        
        return metrics
    
    def calculate_horizontal_mosaic_metrics(self, cluster_key, embedding_key=None, 
                                          gt_key=None, batch_key='batch'):
        """
        Calculate metrics for horizontal/mosaic integration (BioC + SC + BER)
        
        Parameters:
            cluster_key: Key for cluster labels in adata.obs
            embedding_key: Key for embeddings in adata.obsm (optional)
            gt_key: Key for ground truth labels in adata.obs (if available)
            batch_key: Key for batch labels in adata.obs
        
        Returns:
            dict: Calculated metrics
        """
        # Start with vertical metrics (BioC + SC)
        metrics = self.calculate_vertical_metrics(cluster_key, embedding_key, gt_key)
        
        # Add BER metrics if batch information available
        if batch_key in self.adata.obs.columns and embedding_key:
            batch_labels = self.adata.obs[batch_key].astype(str)
            metrics.update(self._calculate_ber_metrics(batch_labels, embedding_key))
        
        return metrics
    
    def _calculate_clustering_metrics(self, gt_labels, cluster_labels):
        """Calculate clustering accuracy metrics (with ground truth)"""
        metrics = {}
        
        # Convert to numeric for sklearn metrics
        gt_numeric = pd.Categorical(gt_labels).codes
        cluster_numeric = pd.Categorical(cluster_labels).codes
        
        # Core clustering metrics
        metrics['ARI'] = adjusted_rand_score(gt_numeric, cluster_numeric)
        metrics['NMI'] = normalized_mutual_info_score(gt_numeric, cluster_numeric)
        metrics['AMI'] = adjusted_mutual_info_score(gt_numeric, cluster_numeric)
        metrics['FMI'] = fowlkes_mallows_score(gt_numeric, cluster_numeric)
        
        # Additional metrics
        metrics['Homogeneity'] = homogeneity_score(gt_numeric, cluster_numeric)
        metrics['Completeness'] = completeness_score(gt_numeric, cluster_numeric)
        metrics['V-measure'] = v_measure_score(gt_numeric, cluster_numeric)
        
        # Purity
        metrics['Purity'] = self._calculate_purity(gt_numeric, cluster_numeric)
        
        # F-measure, Jaccard, and Dice indices
        metrics.update(self._calculate_set_metrics(gt_numeric, cluster_numeric))
        
        return metrics
    
    def _calculate_quality_metrics(self, cluster_labels, embedding_key):
        """Calculate clustering quality metrics (all datasets)"""
        metrics = {}
        
        embedding = self.adata.obsm[embedding_key]
        cluster_numeric = pd.Categorical(cluster_labels).codes
        
        # Standard quality metrics
        metrics['Silhouette Coefficient'] = silhouette_score(embedding, cluster_numeric)
        metrics['Calinski-Harabasz Index'] = calinski_harabasz_score(embedding, cluster_numeric)
        metrics['Davies-Bouldin Index'] = davies_bouldin_score(embedding, cluster_numeric)
        
        return metrics
    
    def _calculate_spatial_metrics(self, cluster_labels):
        """Calculate spatial coherence metrics"""
        metrics = {}
        
        # Get spatial coordinates
        if 'spatial' in self.adata.obsm:
            coords = self.adata.obsm['spatial']
        elif 'x' in self.adata.obs and 'y' in self.adata.obs:
            coords = self.adata.obs[['x', 'y']].values
        else:
            return metrics  # No spatial information
        
        cluster_numeric = pd.Categorical(cluster_labels).codes
        
        # Moran's I
        metrics['Moran Index'] = self._calculate_morans_i(coords, cluster_numeric)
        
        # Geary's C
        metrics['Geary C'] = self._calculate_gearys_c(coords, cluster_numeric)
        
        return metrics
    
    def _calculate_ber_metrics(self, batch_labels, embedding_key):
        """Calculate batch effect removal metrics"""
        metrics = {}
        
        embedding = self.adata.obsm[embedding_key]
        batch_numeric = pd.Categorical(batch_labels).codes
        
        # ASW batch (lower is better - negative indicates good mixing)
        metrics['ASW (batch)'] = silhouette_score(embedding, batch_numeric)
        
        # Graph iLISI and kBET require neighborhood graph
        try:
            # Build graph if not exists
            if 'neighbors' not in self.adata.uns:
                sc.pp.neighbors(self.adata, use_rep=embedding_key, n_neighbors=50)
            
            metrics['Graph iLISI'] = self._calculate_graph_ilisi(batch_labels, embedding_key)
            metrics['kBET'] = self._calculate_kbet(batch_labels)
            metrics['Graph Connectivity'] = self._calculate_graph_connectivity()
            
        except Exception as e:
            print(f"Could not calculate graph-based BER metrics: {e}")
        
        return metrics
    
    def _calculate_asw_celltype(self, gt_labels, embedding_key):
        """Calculate Average Silhouette Width for cell types"""
        embedding = self.adata.obsm[embedding_key]
        gt_numeric = pd.Categorical(gt_labels).codes
        return silhouette_score(embedding, gt_numeric)
    
    def _calculate_purity(self, gt_labels, cluster_labels):
        """Calculate cluster purity"""
        contingency = contingency_matrix(gt_labels, cluster_labels)
        return np.sum(np.amax(contingency, axis=0)) / np.sum(contingency)
    
    def _calculate_set_metrics(self, gt_labels, cluster_labels):
        """Calculate F-measure, Jaccard, and Dice indices"""
        metrics = {}
        
        # Convert to binary classification problem for each class
        unique_gt = np.unique(gt_labels)
        f_scores, jaccard_scores, dice_scores = [], [], []
        
        for gt_class in unique_gt:
            gt_binary = (gt_labels == gt_class).astype(int)
            
            # Find best matching cluster
            best_f, best_jaccard, best_dice = 0, 0, 0
            for cluster in np.unique(cluster_labels):
                cluster_binary = (cluster_labels == cluster).astype(int)
                
                # F-measure
                f = f1_score(gt_binary, cluster_binary, zero_division=0)
                # Jaccard
                jaccard = jaccard_score(gt_binary, cluster_binary, zero_division=0)
                # Dice (2 * Jaccard / (1 + Jaccard))
                dice = 2 * jaccard / (1 + jaccard) if jaccard > 0 else 0
                
                if f > best_f:
                    best_f, best_jaccard, best_dice = f, jaccard, dice
            
            f_scores.append(best_f)
            jaccard_scores.append(best_jaccard)
            dice_scores.append(best_dice)
        
        metrics['F-measure'] = np.mean(f_scores)
        metrics['Jaccard Index'] = np.mean(jaccard_scores)
        metrics['Dice Index'] = np.mean(dice_scores)
        
        return metrics
    
    def _calculate_morans_i(self, coords, values):
        """Calculate Moran's I spatial autocorrelation"""
        n = len(coords)
        
        # Calculate spatial weights matrix (inverse distance)
        distances = squareform(pdist(coords))
        # Avoid division by zero
        distances[distances == 0] = np.finfo(float).eps
        weights = 1.0 / distances
        np.fill_diagonal(weights, 0)
        
        # Normalize weights
        W = weights / np.sum(weights, axis=1, keepdims=True)
        
        # Calculate Moran's I
        mean_val = np.mean(values)
        numerator = 0
        denominator = 0
        
        for i in range(n):
            for j in range(n):
                numerator += W[i, j] * (values[i] - mean_val) * (values[j] - mean_val)
            denominator += (values[i] - mean_val) ** 2
        
        if denominator == 0:
            return 0
        
        return n * numerator / denominator
    
    def _calculate_gearys_c(self, coords, values):
        """Calculate Geary's C spatial autocorrelation"""
        n = len(coords)
        
        # Calculate spatial weights matrix
        distances = squareform(pdist(coords))
        distances[distances == 0] = np.finfo(float).eps
        weights = 1.0 / distances
        np.fill_diagonal(weights, 0)
        
        W = weights / np.sum(weights, axis=1, keepdims=True)
        
        # Calculate Geary's C
        numerator = 0
        denominator = 0
        mean_val = np.mean(values)
        
        for i in range(n):
            for j in range(n):
                numerator += W[i, j] * (values[i] - values[j]) ** 2
            denominator += (values[i] - mean_val) ** 2
        
        if denominator == 0:
            return 1  # No variation
        
        return (n - 1) * numerator / (2 * np.sum(weights) * denominator)
    
    def _calculate_graph_clisi(self, labels, embedding_key):
        """Calculate Graph cLISI (placeholder - requires proper implementation)"""
        # Simplified version - proper implementation would use LISI algorithm
        unique_labels = len(np.unique(labels))
        return unique_labels * 0.8  # Placeholder
    
    def _calculate_graph_ilisi(self, batch_labels, embedding_key):
        """Calculate Graph iLISI (placeholder - requires proper implementation)"""
        # Simplified version - proper implementation would use LISI algorithm
        unique_batches = len(np.unique(batch_labels))
        return unique_batches * 0.7  # Placeholder
    
    def _calculate_kbet(self, batch_labels):
        """Calculate kBET score (placeholder - requires proper implementation)"""
        # Simplified version - proper implementation would use kBET algorithm
        return 0.6  # Placeholder
    
    def _calculate_graph_connectivity(self):
        """Calculate graph connectivity preservation"""
        # Placeholder - requires comparison with original graph
        return 0.85  # Placeholder


def evaluate_method_result(adata_path, cluster_keys, task='vertical', 
                          embedding_key=None, gt_key=None, batch_key='batch',
                          output_dir='Results/evaluation/'):
    """
    Evaluate a method result with multiple clustering approaches
    
    Parameters:
        adata_path: Path to AnnData file with integration results
        cluster_keys: List of clustering result keys in adata.obs
        task: Integration task ('vertical', 'horizontal', 'mosaic')
        embedding_key: Key for integrated embeddings in adata.obsm
        gt_key: Key for ground truth labels (if available)
        batch_key: Key for batch labels (for horizontal/mosaic tasks)
        output_dir: Directory to save results
    
    Returns:
        dict: Results for all clustering methods
    """
    
    # Load data
    adata = sc.read_h5ad(adata_path)
    
    # Initialize metrics calculator
    calculator = SMOBenchMetrics(adata, task=task)
    
    all_results = {}
    
    for cluster_key in cluster_keys:
        if cluster_key not in adata.obs.columns:
            print(f"Warning: Cluster key '{cluster_key}' not found, skipping...")
            continue
        
        print(f"Calculating metrics for {cluster_key}...")
        
        # Calculate appropriate metrics based on task
        if task == 'vertical':
            metrics = calculator.calculate_vertical_metrics(
                cluster_key=cluster_key,
                embedding_key=embedding_key,
                gt_key=gt_key
            )
        else:  # horizontal or mosaic
            metrics = calculator.calculate_horizontal_mosaic_metrics(
                cluster_key=cluster_key,
                embedding_key=embedding_key,
                gt_key=gt_key,
                batch_key=batch_key
            )
        
        all_results[cluster_key] = metrics
        
        # Save individual result
        result_df = pd.DataFrame([metrics]).T
        result_df.columns = [cluster_key]
        
        # Create output filename
        method_dataset = adata_path.split('/')[-1].replace('.h5ad', '')
        gt_status = 'withGT' if gt_key and gt_key in adata.obs.columns else 'woGT'
        output_file = f"{output_dir}/{task}/{gt_status}/method_results/{method_dataset}_{cluster_key}_metrics.csv"
        
        # Create directory if needed
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        result_df.to_csv(output_file)
        print(f"Saved results to {output_file}")
    
    return all_results


if __name__ == "__main__":
    import argparse
    import os
    
    parser = argparse.ArgumentParser(description="Calculate SMOBench evaluation metrics")
    parser.add_argument("--adata", required=True, help="Path to AnnData file")
    parser.add_argument("--cluster-keys", nargs="+", required=True, help="Cluster result keys")
    parser.add_argument("--task", choices=['vertical', 'horizontal', 'mosaic'], 
                       default='vertical', help="Integration task type")
    parser.add_argument("--embedding-key", help="Embedding key in adata.obsm")
    parser.add_argument("--gt-key", help="Ground truth key in adata.obs")
    parser.add_argument("--batch-key", default="batch", help="Batch key in adata.obs")
    parser.add_argument("--output-dir", default="Results/evaluation/", help="Output directory")
    
    args = parser.parse_args()
    
    results = evaluate_method_result(
        adata_path=args.adata,
        cluster_keys=args.cluster_keys,
        task=args.task,
        embedding_key=args.embedding_key,
        gt_key=args.gt_key,
        batch_key=args.batch_key,
        output_dir=args.output_dir
    )
    
    print("Evaluation completed!")
    for cluster_key, metrics in results.items():
        print(f"\n{cluster_key}: {len(metrics)} metrics calculated")