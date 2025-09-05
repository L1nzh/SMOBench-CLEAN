# SMOBench Evaluation Results Structure

## Directory Organization

```
Results/evaluation/
├── vertical/           # Cross-modality, same batch integration
│   ├── withGT/        # Datasets with ground truth annotations
│   │   ├── method_results/     # Individual method-dataset-clustering results
│   │   │   ├── SpatialGlue_HT_S1_leiden_metrics.csv
│   │   │   ├── SpatialGlue_HT_S1_louvain_metrics.csv
│   │   │   ├── SpaMosaic_HT_S1_kmeans_metrics.csv
│   │   │   └── ...
│   │   └── aggregated/        # Aggregated results across methods
│   │       ├── vertical_withGT_summary.csv
│   │       └── vertical_withGT_rankings.csv
│   └── woGT/          # Datasets without ground truth
│       ├── method_results/
│       └── aggregated/
├── horizontal/         # Same modality, cross-batch integration
│   ├── withGT/
│   └── woGT/
└── mosaic/            # Mixed modality and batch integration
    ├── withGT/
    └── woGT/
```

## Evaluation Metrics Categories

### 1. Spatial Coherence (SC)
- **Moran's I**: Spatial autocorrelation of cluster assignments
- **Geary's C**: Local spatial autocorrelation measure

### 2. Biological Conservation (BioC)
#### Clustering Metrics (with GT)
- **ARI**: Adjusted Rand Index
- **NMI**: Normalized Mutual Information
- **AMI**: Adjusted Mutual Information
- **FMI**: Fowlkes-Mallows Index
- **Purity**: Cluster purity measure
- **Homogeneity**: All samples in cluster belong to same class
- **Completeness**: All samples from same class in same cluster
- **V-measure**: Harmonic mean of homogeneity and completeness
- **Jaccard Index**: Intersection over union of cluster assignments
- **Dice Index**: Sørensen-Dice coefficient
- **F-measure**: Harmonic mean of precision and recall

#### Quality Metrics (all datasets)
- **Silhouette Coefficient**: Cluster cohesion and separation
- **Calinski-Harabasz Index**: Ratio of between-cluster to within-cluster dispersion
- **Davies-Bouldin Index**: Average similarity between clusters
- **ASW (cell type)**: Average silhouette width for cell types (with GT)
- **Graph cLISI**: Local inverse Simpson's index for clustering (with GT)

### 3. Batch Effect Removal (BER) - Not used in vertical integration
- **ASW (batch)**: Average silhouette width for batches
- **Graph iLISI**: Local inverse Simpson's index for batch mixing
- **kBET**: k-nearest neighbor batch effect test
- **Graph Connectivity**: Preservation of local neighborhood structure

## File Naming Convention

### Individual Results
`{Method}_{Dataset}_{Clustering}_{GTStatus}_metrics.csv`

Examples:
- `SpatialGlue_HT_S1_leiden_withGT_metrics.csv`
- `SpaMosaic_MB_ATAC_kmeans_woGT_metrics.csv`

### Aggregated Results
`{Task}_{GTStatus}_{AggregationType}.csv`

Examples:
- `vertical_withGT_summary.csv` - All method-dataset combinations
- `vertical_withGT_rankings.csv` - Method rankings by metric
- `vertical_woGT_normalized.csv` - Min-max normalized scores

## Evaluation Pipeline

1. **Method Execution**: Run integration method on dataset
2. **Clustering**: Apply multiple clustering algorithms (leiden, louvain, kmeans, mclust)
3. **Metric Calculation**: Compute appropriate metrics based on task and GT availability
4. **Result Storage**: Save individual results to method_results/
5. **Aggregation**: Combine results and compute rankings in aggregated/