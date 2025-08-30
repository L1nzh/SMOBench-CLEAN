# SMOBench Research Workflow

## Overview
This document outlines the systematic workflow for benchmarking spatial multi-omics integration methods in SMOBench. The workflow ensures reproducible, comprehensive evaluation across methods, datasets, and clustering approaches.

## Phase 1: Environment Setup and Data Preparation

### 1.1 Environment Initialization
```bash
# Activate conda environment
conda activate smobench

# Verify GPU availability
nvidia-smi

# Set CUDA device if needed
export CUDA_VISIBLE_DEVICES=0
```

### 1.2 Dataset Organization
The datasets are organized in the standardized structure:
```
Dataset/
├── withGT/          # Datasets with ground truth annotations
└── woGT/            # Datasets without ground truth
```

### 1.3 Data Verification
```python
# Verify data integrity and structure
conda run -n smobench python Dataset/data_info/generate_data_info.py
```

## Phase 2: Method Execution Pipeline

### 2.1 Individual Method Testing
For each integration method, follow this sequence:

1. **Parameter Configuration**: Set method-specific hyperparameters
2. **Data Loading**: Load RNA and secondary modality (ADT/ATAC)
3. **Preprocessing**: Apply var_names_make_unique() and method-specific preprocessing
4. **Integration**: Execute spatial multi-omics integration
5. **Result Storage**: Save integrated embeddings in AnnData format

### 2.2 Batch Execution Framework
Use the PBS template system for systematic evaluation:

```bash
# Generate PBS scripts from template
python Scripts/generate_pbs_scripts.py

# Submit batch jobs
for pbs in run_*.pbs; do qsub $pbs; done
```

### 2.3 Method-Specific Execution

#### SpatialGlue
```bash
conda run -n smobench python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
  --data_type 10x \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad \
  --save_path Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1.h5ad
```

#### SpaMosaic
```bash
conda run -n smobench python Scripts/integration/SpaMosaic/run_SpaMosaic.py \
  --dataset Human_Tonsils \
  --sample S1 \
  --output_dir Results/adata/SpaMosaic/
```

## Phase 3: Universal Clustering Application

### 3.1 Clustering Interface
Apply multiple clustering methods to each integration result:

```python
from Utils.SMOBench_clustering import universal_clustering

# Load integration result
adata = sc.read_h5ad('Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1.h5ad')

# Apply all clustering methods
for method in ['leiden', 'louvain', 'kmeans', 'mclust']:
    adata = universal_clustering(
        adata, 
        n_clusters=6, 
        used_obsm='spatial_emb',  # or method-specific embedding key
        method=method, 
        key=f'{method}_clusters'
    )

# Save updated result
adata.write('Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1_clustered.h5ad')
```

### 3.2 Parameter Optimization
For leiden/louvain methods:
- Automatic resolution search to match target cluster numbers
- Range: 0.1 to 3.0 with 0.01 increments
- Early stopping when target cluster count is achieved

## Phase 4: Evaluation and Metrics Calculation

### 4.1 Ground Truth Evaluation (withGT datasets)
```python
conda run -n smobench python Scripts/evaluation/eval_with_gt.py \
  --adata_path Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1_clustered.h5ad \
  --gt_column 'final_annot' \
  --output_path Results/metrics/SpatialGlue_HT_S1_with_gt.csv
```

Metrics calculated:
- ARI (Adjusted Rand Index)
- NMI (Normalized Mutual Information)  
- FMI (Fowlkes-Mallows Index)
- Silhouette Score

### 4.2 Unsupervised Evaluation (woGT datasets)
```python
conda run -n smobench python Scripts/evaluation/eval_without_gt.py \
  --adata_path Results/adata/SpatialGlue/MS/Spleen1/SpatialGlue_MS_Spleen1_clustered.h5ad \
  --output_path Results/metrics/SpatialGlue_MS_Spleen1_without_gt.csv
```

Metrics calculated:
- Silhouette Score
- Moran's I (spatial autocorrelation)
- Spatial Consistency Index

### 4.3 Batch Processing
```bash
# Evaluate all results systematically
conda run -n smobench python Scripts/evaluation/batch_evaluate.py \
  --results_dir Results/adata/ \
  --output_dir Results/metrics/
```

## Phase 5: Result Aggregation and Analysis

### 5.1 Metric Standardization
Within each dataset:
1. Calculate min/max values for each metric across all methods
2. Apply min-max normalization: `(value - min) / (max - min)`
3. Handle directionality (higher-is-better vs lower-is-better metrics)

### 5.2 Cross-Dataset Aggregation
```python
conda run -n smobench python Scripts/evaluation/generate_final_results.py \
  --metrics_dir Results/metrics/ \
  --output_path Results/final_benchmark.csv
```

Aggregation steps:
1. **Method-Dataset Averaging**: Mean across clustering methods for each method-dataset pair
2. **Batch Averaging**: Mean across sample replicates (S1, S2) within tissue types  
3. **Task-Specific Averaging**: Separate averages for vertical, horizontal, mosaic tasks
4. **Overall Ranking**: Weighted average across all evaluation criteria

### 5.3 Statistical Analysis
- Significance testing between methods (paired t-test, Wilcoxon)
- Effect size calculation (Cohen's d)
- Confidence intervals for performance metrics
- Robustness analysis across hyperparameter variations

## Phase 6: Visualization and Reporting

### 6.1 Spatial Visualizations
```python
conda run -n smobench python Scripts/visualization/plot_spatial_results.py \
  --adata_path Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1_clustered.h5ad \
  --output_dir Results/plots/SpatialGlue/Human_Tonsils/S1/
```

Generated plots:
- Spatial scatter plots colored by cluster assignment
- UMAP embeddings of integrated representations
- Side-by-side comparison: predicted vs ground truth (if available)

### 6.2 Performance Heatmaps
```python
conda run -n smobench python Scripts/visualization/plot_benchmark_heatmap.py \
  --results_path Results/final_benchmark.csv \
  --output_path Results/plots/benchmark_heatmap.png
```

### 6.3 Method Comparison Plots
- Box plots of performance distributions
- Radar charts for multi-metric comparison
- Scatter plots of metric correlations
- Runtime and memory usage comparison

## Phase 7: Quality Control and Reproducibility

### 7.1 Result Validation
- Cross-check results with known method implementations
- Verify metric calculations with reference implementations
- Test on synthetic datasets with known ground truth

### 7.2 Documentation Updates
- Update status.md with current progress
- Document any parameter modifications or workflow changes
- Record computational resource usage and runtime statistics

### 7.3 Version Control
```bash
# Commit major workflow milestones
git add .
git commit -m "Complete benchmark evaluation for [method/dataset]"
git push origin main
```

## Phase 8: Advanced Analysis (Optional)

### 8.1 Method Failure Analysis
- Identify datasets where methods consistently underperform
- Analyze correlation between dataset characteristics and method performance
- Document method-specific limitations and optimal use cases

### 8.2 Ensemble Methods
- Combine predictions from multiple integration methods
- Evaluate ensemble performance vs individual methods
- Identify complementary method combinations

### 8.3 Hyperparameter Sensitivity
- Systematic variation of key hyperparameters
- Performance robustness analysis
- Optimal hyperparameter identification for each tissue type

## Error Handling and Troubleshooting

### Common Issues
1. **GPU Memory**: Use `nvidia-smi` to check availability, clear cache with `torch.cuda.empty_cache()`
2. **R Interface**: Ensure rpy2 and mclust are properly installed for mclust clustering
3. **Data Format**: Verify AnnData objects have required obsm keys for embeddings
4. **Missing Dependencies**: Use conda to install method-specific requirements

### Recovery Procedures
- Checkpoint intermediate results after each major processing step
- Implement automatic retry logic for transient failures
- Maintain backup copies of processed datasets and results