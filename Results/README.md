# SMOBench Results

This folder contains the outputs from SMOBench spatial multi-omics integration methods.

## Folder Structure

```
Results/
├── adata/                              # Processed AnnData results
│   ├── SpatialGlue/                   # SpatialGlue integration results
│   │   ├── HLN/                       # Human Lymph Nodes
│   │   │   ├── A1/                    # Sample A1
│   │   │   │   └── SpatialGlue_HLN_A1.h5ad
│   │   │   └── D1/                    # Sample D1
│   │   ├── HT/                        # Human Tonsils  
│   │   │   ├── S1/, S2/, S3/          # Multiple samples
│   │   ├── Mouse_Thymus/              # Mouse thymus samples
│   │   ├── Mouse_Spleen/              # Mouse spleen samples
│   │   ├── Mouse_Brain/               # Mouse brain samples
│   │   └── MISAR_S1/, MISAR_S2/       # Mouse embryo samples
│   ├── PRAGA/                         # PRAGA integration results
│   │   └── [same structure as above]
│   ├── COSMOS/                        # COSMOS integration results  
│   │   └── [same structure as above]
│   ├── PRESENT/                       # PRESENT integration results
│   │   └── [same structure as above]
│   ├── SpaMultiVAE/                   # SpaMultiVAE integration results
│   │   └── [same structure as above]
│   └── SpaMosaic/                     # SpaMosaic integration results
│       └── [same structure as above]
├── plot/                              # Visualization outputs
│   ├── SpatialGlue/                   # SpatialGlue plots
│   │   ├── Human_Lymph_Nodes/         
│   │   │   ├── A1/
│   │   │   │   ├── clustering_mclust_umap_spatial.png
│   │   │   │   ├── clustering_leiden_umap_spatial.png
│   │   │   │   ├── clustering_louvain_umap_spatial.png
│   │   │   │   └── clustering_kmeans_umap_spatial.png
│   │   │   └── D1/
│   │   ├── Human_Tonsils/
│   │   │   ├── S1/, S2/, S3/
│   │   ├── Mouse_Thymus/
│   │   ├── Mouse_Spleen/  
│   │   ├── Mouse_Brain/
│   │   └── Mouse_Embryos_S1/, Mouse_Embryos_S2/
│   ├── PRAGA/                         # PRAGA plots
│   │   └── [same structure as above]
│   ├── COSMOS/                        # COSMOS plots
│   │   └── [same structure as above]  
│   ├── PRESENT/                       # PRESENT plots
│   │   └── [same structure as above]
│   ├── SpaMultiVAE/                   # SpaMultiVAE plots
│   │   └── [same structure as above]
│   └── SpaMosaic/                     # SpaMosaic plots
│       └── [same structure as above]
└── evaluation/                        # Evaluation metrics and reports
    ├── clustering_results.csv         # Clustering performance metrics
    ├── integration_metrics.csv        # Integration quality metrics
    └── method_comparison.csv          # Cross-method comparison
```

## File Types

### AnnData Results (`*.h5ad`)
- **Location**: `Results/adata/{Method}/{Dataset}/{Sample}/`
- **Content**: Integrated spatial multi-omics data with embeddings
- **Key Fields**:
  - `adata.obsm['{Method}']`: Integration embeddings
  - `adata.obs['{clustering_method}']`: Clustering labels
  - `adata.obsm['spatial']`: Spatial coordinates
  - `adata.obsm['X_umap']`: UMAP coordinates for visualization

### Visualization Plots (`*.png`)
- **Location**: `Results/plot/{Method}/{Dataset}/{Sample}/`
- **Types**:
  - `clustering_{method}_umap_spatial.png`: UMAP + spatial plots
  - Methods include: mclust, leiden, louvain, kmeans
- **Format**: High-resolution (300 DPI) publication-ready figures

### Evaluation Metrics (`*.csv`)
- **Location**: `Results/evaluation/`
- **Content**: Quantitative performance metrics for method comparison
- **Metrics**: ARI, NMI, Silhouette score, spatial coherence, etc.

## Data Formats

### Integration Results Structure
Each `.h5ad` file contains:
- **X**: Original or processed expression matrix
- **obs**: Cell metadata with clustering results
- **obsm**: Embeddings and coordinates
  - `{Method}`: Method-specific integration embeddings
  - `spatial`: Spatial coordinates
  - `X_umap`: UMAP visualization coordinates
- **var**: Feature metadata

### Clustering Results
For each clustering method (mclust, leiden, louvain, kmeans):
- Stored in `adata.obs['{method}']` as categorical labels
- Consistent cluster numbering across methods
- Used for downstream evaluation and comparison

## Generating Results

To regenerate these results, run the integration scripts:

```bash
# Run individual methods
bash Scripts/integration/SpatialGlue/run.sh
bash Scripts/integration/PRAGA/run.sh
bash Scripts/integration/COSMOS/run.sh
bash Scripts/integration/PRESENT/run.sh
bash Scripts/integration/SpaMultiVAE/run.sh

# Or use the task-based system
python Scripts/run_tasks.py --config tasks_sample.csv
```

## Note

Result files are excluded from the repository due to large file sizes. Users should generate their own results using the provided scripts and datasets.