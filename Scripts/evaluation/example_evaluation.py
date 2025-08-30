#!/usr/bin/env python3
"""
Example evaluation script demonstrating SMOBench evaluation workflow

This script shows how to:
1. Load integration results
2. Apply universal clustering  
3. Calculate comprehensive evaluation metrics
4. Save results in organized structure
"""

import scanpy as sc
import pandas as pd
from Utils.SMOBench_clustering import batch_clustering
from Scripts.evaluation.metrics_calculator import SMOBenchMetrics
import os

def run_example_evaluation():
    """
    Example evaluation workflow for a SpatialGlue result on Human Tonsils S1
    """
    print("=== SMOBench Evaluation Example ===\n")
    
    # Step 1: Load integration result (this would be the output from SpatialGlue)
    print("1. Loading integration result...")
    # In practice, this would be: adata = sc.read_h5ad('path/to/spatialglue/result.h5ad')
    # For demo, we'll create a mock result
    print("   (Using mock data for demonstration)")
    
    # Step 2: Apply universal clustering with multiple methods
    print("\n2. Applying universal clustering...")
    
    # Mock clustering results for demonstration
    cluster_methods = ['leiden', 'louvain', 'kmeans', 'mclust']
    print(f"   Applying clustering methods: {cluster_methods}")
    print("   - leiden: Resolution search -> 6 clusters found")
    print("   - louvain: Resolution search -> 6 clusters found") 
    print("   - kmeans: K-means with k=6")
    print("   - mclust: Gaussian mixture with k=6")
    
    # Step 3: Calculate evaluation metrics
    print("\n3. Calculating evaluation metrics...")
    
    # Example for vertical integration with ground truth
    print("   Task: Vertical Integration")
    print("   Ground Truth: Available (withGT)")
    print("   Calculating 19 comprehensive metrics...")
    
    # Mock metrics results matching your examples
    sample_metrics_withGT = {
        'leiden': {
            'asw_celltype': 0.511,
            'graph_clisi': 0.886,
            'ARI': 0.104,
            'NMI': 0.274,
            'FMI': 0.313,
            'Silhouette Coefficient': 0.147,
            'Calinski-Harabasz Index': 893.968,
            'Davies-Bouldin Index': 1.673,
            'Purity': 0.750,
            'AMI': 0.272,
            'Homogeneity': 0.450,
            'Completeness': 0.197,
            'V-measure': 0.274,
            'F-measure': 0.241,
            'Jaccard Index': 0.137,
            'Dice Index': 0.241,
            'Moran Index': 0.448,
            'Geary C': 0.554
        }
    }
    
    sample_metrics_woGT = {
        'leiden': {
            'Silhouette Coefficient': 0.147,
            'Calinski-Harabasz Index': 893.968,
            'Davies-Bouldin Index': 1.673,
            'Moran Index': 0.448,
            'Geary C': 0.554
        }
    }
    
    # Step 4: Display results
    print("\n4. Results Summary:")
    print("\n   With Ground Truth (19 metrics):")
    for metric, value in sample_metrics_withGT['leiden'].items():
        print(f"     {metric:25}: {value:.3f}")
    
    print("\n   Without Ground Truth (5 metrics):")
    for metric, value in sample_metrics_woGT['leiden'].items():
        print(f"     {metric:25}: {value:.3f}")
    
    # Step 5: Save results
    print("\n5. Saving results...")
    
    # Create output directories
    output_dirs = [
        'Results/evaluation/vertical/withGT/method_results/',
        'Results/evaluation/vertical/woGT/method_results/'
    ]
    
    for dir_path in output_dirs:
        os.makedirs(dir_path, exist_ok=True)
        print(f"   Created directory: {dir_path}")
    
    # Save sample results
    for cluster_method in ['leiden', 'louvain', 'kmeans', 'mclust']:
        # With GT results
        withgt_df = pd.DataFrame([sample_metrics_withGT['leiden']]).T
        withgt_df.columns = [f'SpatialGlue_HT_S1_{cluster_method}']
        withgt_file = f'Results/evaluation/vertical/withGT/method_results/SpatialGlue_HT_S1_{cluster_method}_metrics.csv'
        withgt_df.to_csv(withgt_file, header=False)
        
        # Without GT results  
        wogt_df = pd.DataFrame([sample_metrics_woGT['leiden']]).T
        wogt_df.columns = [f'SpatialGlue_HT_S1_{cluster_method}']
        wogt_file = f'Results/evaluation/vertical/woGT/method_results/SpatialGlue_HT_S1_{cluster_method}_metrics.csv'
        wogt_df.to_csv(wogt_file, header=False)
        
        print(f"   Saved: {withgt_file}")
        print(f"   Saved: {wogt_file}")
    
    print("\n=== Evaluation Complete ===")
    print(f"✅ Processed 1 method × 1 dataset × 4 clustering approaches")
    print(f"✅ Generated 8 evaluation result files")
    print(f"✅ Ready for aggregation and benchmarking")

def show_evaluation_structure():
    """Show the evaluation framework structure"""
    print("\n=== SMOBench Evaluation Framework ===")
    print("""
Evaluation Categories:
1. Spatial Coherence (SC)
   - Moran's I: Spatial autocorrelation
   - Geary's C: Local spatial association

2. Biological Conservation (BioC)  
   - Clustering Accuracy (with GT): ARI, NMI, AMI, FMI, etc.
   - Quality Metrics (all): Silhouette, Calinski-Harabasz, Davies-Bouldin
   - Additional BioC (with GT): ASW cell type, Graph cLISI

3. Batch Effect Removal (BER) - Horizontal/Mosaic only
   - ASW batch, Graph iLISI, kBET, Graph Connectivity

Task-Specific Evaluation:
• Vertical Integration: BioC + SC
  - With GT: 19 metrics (full evaluation)
  - Without GT: 5 metrics (quality + spatial only)

• Horizontal/Mosaic Integration: BioC + SC + BER  
  - With GT: 19+ metrics (includes BER metrics)
  - Without GT: 5+ metrics (includes BER metrics)
""")

def show_usage_examples():
    """Show practical usage examples"""
    print("\n=== Usage Examples ===")
    
    print("\n1. Evaluate single method result:")
    print("""
conda run -n smobench python Scripts/evaluation/metrics_calculator.py \\
  --adata Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1_clustered.h5ad \\
  --cluster-keys leiden louvain kmeans mclust \\
  --task vertical \\
  --embedding-key spatial_emb \\
  --gt-key final_annot \\
  --output-dir Results/evaluation/
""")
    
    print("\n2. Batch evaluation of all results:")
    print("""
for task in vertical horizontal mosaic; do
  for gt_status in withGT woGT; do
    conda run -n smobench python Scripts/evaluation/batch_evaluate.py \\
      --task $task \\
      --gt_status $gt_status \\
      --results_dir Results/adata/ \\
      --output_dir Results/evaluation/$task/$gt_status/
  done
done
""")
    
    print("\n3. Generate final benchmark results:")
    print("""
conda run -n smobench python Scripts/evaluation/generate_final_results.py \\
  --evaluation_dir Results/evaluation/ \\
  --output_path Results/SMOBench_final_benchmark.csv
""")

if __name__ == "__main__":
    run_example_evaluation()
    show_evaluation_structure()
    show_usage_examples()