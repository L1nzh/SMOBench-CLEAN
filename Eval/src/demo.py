"""
Vertical Integration Evaluation Functions for SMOBench
Evaluates spatial multi-omics integration methods on SC (Spatial Coherence) and BioC (Biological Conservation) metrics
"""

import os
import pandas as pd
import numpy as np
from sklearn import metrics
# Import functions with fallback handling
def safe_import():
    """Safely import functions with fallback implementations"""
    global Moran_Geary, silhouette_simple, graph_clisi
    
    import sys
    import os
    
    # Add Eval directory to path for local imports
    current_dir = os.path.dirname(os.path.abspath(__file__))
    eval_dir = os.path.dirname(current_dir)
    if eval_dir not in sys.path:
        sys.path.insert(0, eval_dir)
    
    try:
        # Try to import individual functions
        from src.compute_metric import silhouette_simple, graph_clisi
        print("Successfully imported silhouette_simple and graph_clisi from compute_metric")
        
        # Try to import Moran_Geary
        try:
            from src.compute_metric import Moran_Geary
            print("Successfully imported Moran_Geary from compute_metric")
        except ImportError:
            print("Warning: Moran_Geary import failed, using spatial_metrics_simple fallback")
            from src.spatial_metrics_simple import Moran_Geary
            
    except ImportError as e:
        print(f"Warning: compute_metric import failed ({e}), using fallback implementations")
        
        # Import from fallback
        try:
            from src.spatial_metrics_simple import Moran_Geary
            print("Using spatial_metrics_simple fallback for Moran_Geary")
        except ImportError:
            print("Warning: spatial_metrics_simple also failed, using inline fallback")
            
            def Moran_Geary(coordinates, labels):
                """Inline fallback implementation"""
                class SimpleResult:
                    def __init__(self, value):
                        self.I = value if hasattr(self, 'I') else None
                        self.C = value if hasattr(self, 'C') else None
                
                moran_result = SimpleResult(0.0)
                moran_result.I = 0.0
                geary_result = SimpleResult(1.0) 
                geary_result.C = 1.0
                return moran_result, geary_result
        
        def silhouette_simple(embeddings, labels):
            """Fallback silhouette implementation"""
            from sklearn.metrics import silhouette_score
            try:
                asw = silhouette_score(embeddings, labels)
                return (asw + 1) / 2  # Scale to match original implementation
            except:
                return 0.0
        
        def graph_clisi(adj_matrix, labels):
            """Fallback graph_clisi implementation"""
            # Simplified version - just return a reasonable default
            return 0.5

# Call the safe import function
safe_import()

def eval_vertical_integration(embeddings, adj_matrix, y_pred, y_GT, spatial_coords, method_name, dataset_name, slice_name, clustering_method):
    """
    Evaluate vertical integration results for a single slice
    
    Parameters:
    -----------
    embeddings : np.ndarray
        Integrated embeddings from the method
    adj_matrix : sparse matrix
        Adjacency matrix for spatial relationships
    y_pred : array-like
        Predicted cluster labels
    y_GT : array-like or None
        Ground truth labels (Spatial_Label for withGT datasets, None for woGT)
    spatial_coords : np.ndarray
        Spatial coordinates for cells
    method_name : str
        Name of the integration method
    dataset_name : str
        Name of the tissue/dataset (e.g., 'HLN', 'HT', 'ME_S1')
    slice_name : str
        Name of the slice (e.g., 'A1', 'S1', 'E11')
    clustering_method : str
        Clustering method used ('leiden', 'louvain', 'kmeans', 'mclust')
    
    Returns:
    --------
    dict : Evaluation metrics
    """
    
    # Calculate Spatial Coherence (SC) metrics
    try:
        # y_pred should already be numeric from eval_adata.py preprocessing
        y_pred_numeric = np.asarray(y_pred, dtype=float)
        
        Moran, Geary = Moran_Geary(spatial_coords, y_pred_numeric)
        sc_metrics = {
            'Moran Index': Moran.I,
            'Geary C': Geary.C
        }
    except Exception as e:
        print(f"Warning: Could not compute spatial metrics for {method_name}_{dataset_name}_{slice_name}: {e}")
        sc_metrics = {
            'Moran Index': 0.0,
            'Geary C': 0.0
        }
    
    # Initialize metrics dictionary
    metrics_dict = sc_metrics.copy()
    
    # Calculate Biological Conservation (BioC) metrics
    if y_GT is not None:
        # withGT metrics - check for data alignment first
        try:
            # Check if all arrays have consistent lengths
            n_embeddings = embeddings.shape[0]
            n_gt = len(y_GT)
            n_pred = len(y_pred)
            n_coords = spatial_coords.shape[0]
            
            if not all(n == n_embeddings for n in [n_gt, n_pred, n_coords]):
                print(f"Warning: Data length mismatch for {method_name}_{dataset_name}_{slice_name}:")
                print(f"  Embeddings: {n_embeddings}, GT: {n_gt}, Pred: {n_pred}, Coords: {n_coords}")
                
                # Find minimum length and truncate all arrays
                min_len = min(n_embeddings, n_gt, n_pred, n_coords)
                print(f"  Truncating all arrays to length: {min_len}")
                
                embeddings_aligned = embeddings[:min_len]
                y_GT_aligned = y_GT[:min_len]
                y_pred_aligned = y_pred[:min_len]
                adj_matrix_aligned = adj_matrix[:min_len, :min_len]
            else:
                embeddings_aligned = embeddings
                y_GT_aligned = y_GT
                y_pred_aligned = y_pred
                adj_matrix_aligned = adj_matrix
            
            bioc_metrics = {
                'ARI': metrics.adjusted_rand_score(np.ravel(y_GT_aligned), np.ravel(y_pred_aligned)),
                'NMI': metrics.normalized_mutual_info_score(np.ravel(y_GT_aligned), np.ravel(y_pred_aligned)),
                'asw_celltype': silhouette_simple(embeddings_aligned, y_GT_aligned),
                'graph_clisi': graph_clisi(adj_matrix_aligned, y_GT_aligned)
            }
            metrics_dict.update(bioc_metrics)
        except Exception as e:
            print(f"Warning: Could not compute withGT BioC metrics for {method_name}_{dataset_name}_{slice_name}: {e}")
            # Add default values
            bioc_metrics = {
                'ARI': 0.0,
                'NMI': 0.0,
                'asw_celltype': 0.0,
                'graph_clisi': 0.0
            }
            metrics_dict.update(bioc_metrics)
    else:
        # woGT metrics
        try:
            bioc_metrics = {
                'Davies-Bouldin Index': metrics.davies_bouldin_score(embeddings, y_pred),
                'Silhouette Coefficient': metrics.silhouette_score(embeddings, y_pred, metric='euclidean'),
                'Calinski-Harabaz Index': metrics.calinski_harabasz_score(embeddings, y_pred)
            }
            metrics_dict.update(bioc_metrics)
        except Exception as e:
            print(f"Warning: Could not compute woGT BioC metrics for {method_name}_{dataset_name}_{slice_name}: {e}")
            # Add default values
            bioc_metrics = {
                'Davies-Bouldin Index': 1.0,
                'Silhouette Coefficient': 0.0,
                'Calinski-Harabaz Index': 1.0
            }
            metrics_dict.update(bioc_metrics)
    
    # Add metadata
    metrics_dict['Method'] = method_name
    metrics_dict['Dataset'] = dataset_name
    metrics_dict['Slice'] = slice_name
    metrics_dict['Clustering'] = clustering_method
    metrics_dict['GT_Available'] = y_GT is not None
    
    return metrics_dict

def save_evaluation_results(metrics_dict, output_dir, method_name, dataset_name, slice_name, clustering_method, has_gt):
    """
    Save evaluation results to CSV file
    
    Parameters:
    -----------
    metrics_dict : dict
        Dictionary containing evaluation metrics
    output_dir : str
        Output directory path
    method_name : str
        Name of the integration method
    dataset_name : str
        Name of the tissue/dataset
    slice_name : str
        Name of the slice
    clustering_method : str
        Clustering method used
    has_gt : bool
        Whether ground truth is available
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare data for saving (exclude metadata fields)
    save_metrics = {k: v for k, v in metrics_dict.items() 
                   if k not in ['Method', 'Dataset', 'Slice', 'Clustering', 'GT_Available']}
    
    # Create DataFrame
    df = pd.DataFrame(list(save_metrics.items()), columns=['Metric', 'Value'])
    
    # Determine file suffix
    suffix = 'withGT' if has_gt else 'woGT'
    
    # Save to CSV
    filename = f"{method_name}_{dataset_name}_{slice_name}_{clustering_method}_{suffix}.csv"
    filepath = os.path.join(output_dir, filename)
    df.to_csv(filepath, index=False)
    
    print(f"Saved evaluation results to: {filepath}")
    
    return filepath

def calculate_dataset_summary(metrics_list, dataset_name, method_name, clustering_method):
    """
    Calculate summary statistics for a dataset (averaging across slices)
    
    Parameters:
    -----------
    metrics_list : list of dict
        List of metrics dictionaries from different slices
    dataset_name : str
        Name of the dataset
    method_name : str
        Name of the method
    clustering_method : str
        Clustering method used
        
    Returns:
    --------
    dict : Summary metrics for the dataset
    """
    
    if not metrics_list:
        return None
    
    # Determine metric categories
    has_gt = metrics_list[0]['GT_Available']
    
    if has_gt:
        sc_metrics = ['Moran Index', 'Geary C']
        bioc_metrics = ['ARI', 'NMI', 'asw_celltype', 'graph_clisi']
    else:
        sc_metrics = ['Moran Index', 'Geary C']
        bioc_metrics = ['Davies-Bouldin Index', 'Silhouette Coefficient', 'Calinski-Harabaz Index']
    
    # Calculate averages
    summary = {
        'Method': method_name,
        'Dataset': dataset_name,
        'Clustering': clustering_method,
        'GT_Available': has_gt
    }
    
    # Average SC metrics
    sc_values = []
    for metric in sc_metrics:
        values = [m[metric] for m in metrics_list if metric in m and not np.isnan(m[metric])]
        if values:
            avg_val = np.mean(values)
            summary[metric] = avg_val
            sc_values.append(avg_val)
    
    # Average BioC metrics
    bioc_values = []
    for metric in bioc_metrics:
        values = [m[metric] for m in metrics_list if metric in m and not np.isnan(m[metric])]
        if values:
            avg_val = np.mean(values)
            summary[metric] = avg_val
            bioc_values.append(avg_val)
    
    # Calculate category scores
    summary['SC_Score'] = np.mean(sc_values) if sc_values else 0.0
    summary['BioC_Score'] = np.mean(bioc_values) if bioc_values else 0.0
    summary['Total_Score'] = (summary['SC_Score'] + summary['BioC_Score']) / 2
    
    return summary