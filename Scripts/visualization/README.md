# SMOBench Visualization Module

This module provides comprehensive visualization tools for SMOBench integration results, focusing on UMAP and spatial plotting capabilities.

## Overview

The visualization module consists of three main components:
- **`plot_umap_spatial.py`**: Core visualization script for individual files
- **`run_visualization.sh`**: Shell wrapper for easy command-line usage
- **`batch_visualization.py`**: Batch processing for multiple files

## Quick Start - Generate All UMAP Plots

To create UMAP plots for all your integration results at once:

```bash
# Navigate to SMOBench root directory
cd /path/to/SMOBench-CLEAN

# Generate all UMAP plots automatically
conda run -n smobench python Scripts/visualization/batch_visualization.py \
  --results_dir Results/adata \
  --auto_detect_clustering \
  --summary_file umap_generation_summary.csv
```

This will:
- Process all `.h5ad` files in `Results/adata/`
- Auto-detect integration methods (SpatialGlue, PRAGA, SpaMosaic, etc.)
- Auto-detect embedding keys and clustering methods
- Generate UMAP + spatial plots for all clustering methods
- Save plots to `Results/plot/Method/Dataset/Subset/`
- Create a summary report

## Individual File Visualization

### Using Shell Script (Recommended)

```bash
# Basic usage
bash Scripts/visualization/run_visualization.sh \
  -i Results/adata/SpatialGlue/Human_Lymph_Nodes/A1/adata_integrated.h5ad \
  -m SpatialGlue \
  -e SpatialGlue

# With custom clustering methods
bash Scripts/visualization/run_visualization.sh \
  -i Results/adata/PRAGA/MISAR/sample1/adata_integrated.h5ad \
  -m PRAGA \
  -e PRAGA \
  -c "mclust,leiden"

# With custom settings
bash Scripts/visualization/run_visualization.sh \
  -i Results/adata/SpaMosaic/simulation/data1/adata_integrated.h5ad \
  -m SpaMosaic \
  -e merged_emb \
  --point-size 15 \
  --n-neighbors 50
```

### Using Python Script Directly

```bash
conda run -n smobench python Scripts/visualization/plot_umap_spatial.py \
  --adata_path Results/adata/SpatialGlue/Human_Lymph_Nodes/A1/adata_integrated.h5ad \
  --method SpatialGlue \
  --embedding_key SpatialGlue \
  --clustering_methods "mclust,leiden,louvain,kmeans"
```

## Batch Processing Options

### 1. Process All Results

```bash
# Process everything in Results/adata/
conda run -n smobench python Scripts/visualization/batch_visualization.py \
  --results_dir Results/adata \
  --auto_detect_clustering
```

### 2. Process Specific Methods

```bash
# Only SpatialGlue results
conda run -n smobench python Scripts/visualization/batch_visualization.py \
  --results_dir Results/adata \
  --method SpatialGlue \
  --embedding_key SpatialGlue \
  --clustering_methods "mclust,leiden,louvain,kmeans"
```

### 3. Process from File List

```bash
# Create file list
find Results/adata -name "*SpatialGlue*.h5ad" > spatialglue_files.txt

# Process specific files
conda run -n smobench python Scripts/visualization/batch_visualization.py \
  --file_list spatialglue_files.txt \
  --auto_detect_clustering
```

## Command-Line Arguments

### Shell Script Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-i, --input` | Input AnnData file path | Required |
| `-m, --method` | Integration method name | Required |
| `-e, --embedding` | Embedding key in adata.obsm | Required |
| `-d, --dataset` | Dataset name | Auto-detected |
| `-s, --subset` | Subset/sample name | Auto-detected |
| `-c, --clustering` | Clustering methods (comma-separated) | `mclust,leiden,louvain,kmeans` |
| `-o, --output` | Output plot directory | `Results/plot` |
| `--point-size` | Point size for plots | `20` |
| `--n-neighbors` | Number of neighbors for UMAP | `30` |
| `--no-spatial-flip` | Don't flip Y coordinates | `false` |
| `--force-recompute` | Force UMAP recomputation | `false` |
| `--compare-emb` | Compare multiple embeddings | None |

### Batch Script Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--results_dir` | Directory with AnnData files | Required* |
| `--file_list` | Text file with file paths | Required* |
| `--pattern` | Glob pattern for file search | `**/*.h5ad` |
| `--method` | Force specific method name | Auto-detected |
| `--embedding_key` | Force specific embedding key | Auto-detected |
| `--clustering_methods` | Clustering methods | Auto-detected |
| `--auto_detect_clustering` | Auto-detect clustering methods | `false` |
| `--output_dir` | Output plot directory | `Results/plot` |
| `--summary_file` | Summary report file | `batch_visualization_summary.csv` |
| `--stop_on_error` | Stop on first error | `false` |

*Either `--results_dir` or `--file_list` is required.

## Output Structure

Plots are saved to the following directory structure:
```
Results/plot/
├── SpatialGlue/
│   ├── Human_Lymph_Nodes/
│   │   ├── A1/
│   │   │   ├── clustering_mclust_umap_spatial.png
│   │   │   ├── clustering_leiden_umap_spatial.png
│   │   │   ├── clustering_louvain_umap_spatial.png
│   │   │   └── clustering_kmeans_umap_spatial.png
│   │   └── A2/
│   └── MISAR/
├── PRAGA/
│   ├── Human_Lymph_Nodes/
│   └── MISAR/
└── SpaMosaic/
    └── simulation/
```

## Generated Plot Types

For each clustering method, the script generates:

1. **UMAP + Spatial Plot** (`clustering_{method}_umap_spatial.png`):
   - Left panel: UMAP projection colored by clusters
   - Right panel: Spatial coordinates colored by clusters

2. **UMAP Only** (`clustering_{method}_umap.png`):
   - Generated when spatial coordinates are not available

3. **Embedding Comparison** (`embedding_{key}_umap.png`):
   - When using `--compare_embeddings` option

## Supported Integration Methods

The visualization module automatically detects and supports:

| Method | Default Embedding Key | Common Alt Keys |
|--------|----------------------|-----------------|
| SpatialGlue | `SpatialGlue` | `spatial_emb` |
| PRAGA | `PRAGA` | `PRAGA_emb` |
| SpaMosaic | `merged_emb` | `SpaMosaic` |
| COSMOS | `COSMOS` | `cosmos_emb` |
| PRESENT | `PRESENT` | `present_emb` |
| SpaMultiVAE | `SpaMultiVAE` | `spa_emb` |

## Supported Clustering Methods

- **mclust**: R-based Gaussian mixture model clustering
- **leiden**: Leiden community detection algorithm
- **louvain**: Louvain community detection algorithm
- **kmeans**: K-means clustering

## Troubleshooting

### Common Issues

1. **"Embedding key not found"**
   ```bash
   # Check available keys in your AnnData file
   python -c "import scanpy as sc; adata = sc.read_h5ad('your_file.h5ad'); print(list(adata.obsm.keys()))"
   ```

2. **"No clustering results found"**
   ```bash
   # Check available clustering columns
   python -c "import scanpy as sc; adata = sc.read_h5ad('your_file.h5ad'); print([col for col in adata.obs.columns if col in ['mclust', 'leiden', 'louvain', 'kmeans']])"
   ```

3. **UMAP computation fails**
   ```bash
   # Use force recompute option
   bash Scripts/visualization/run_visualization.sh \
     -i your_file.h5ad -m Method -e embedding_key --force-recompute
   ```

4. **Spatial coordinates look wrong**
   ```bash
   # Try without Y-axis flipping
   bash Scripts/visualization/run_visualization.sh \
     -i your_file.h5ad -m Method -e embedding_key --no-spatial-flip
   ```

### Getting Help

```bash
# Show shell script help
bash Scripts/visualization/run_visualization.sh --help

# Show batch script help
python Scripts/visualization/batch_visualization.py --help

# Show Python script help
python Scripts/visualization/plot_umap_spatial.py --help
```

## Examples for Common Use Cases

### Visualize SpatialGlue Results
```bash
bash Scripts/visualization/run_visualization.sh \
  -i Results/adata/SpatialGlue/Human_Lymph_Nodes/A1/adata_integrated.h5ad \
  -m SpatialGlue \
  -e SpatialGlue
```

### Visualize PRAGA Results
```bash
bash Scripts/visualization/run_visualization.sh \
  -i Results/adata/PRAGA/MISAR/sample1/adata_integrated.h5ad \
  -m PRAGA \
  -e PRAGA
```

### Batch Process All Results with Custom Settings
```bash
conda run -n smobench python Scripts/visualization/batch_visualization.py \
  --results_dir Results/adata \
  --auto_detect_clustering \
  --point_size 15 \
  --n_neighbors 50 \
  --summary_file custom_umap_summary.csv
```

### Compare Multiple Embeddings
```bash
bash Scripts/visualization/run_visualization.sh \
  -i Results/adata/SpaMosaic/simulation/data1/adata_integrated.h5ad \
  -m SpaMosaic \
  -e merged_emb \
  --compare-emb "spatial_emb,merged_emb,latent_emb"
```

## Integration with SMOBench Workflow

This visualization module is designed to work seamlessly with the SMOBench workflow:

1. **After Integration**: Run integration methods (SpatialGlue, PRAGA, etc.)
2. **Generate Plots**: Use this visualization module to create UMAP plots
3. **Evaluation**: Proceed with evaluation metrics and analysis

The module follows SMOBench conventions:
- Uses `conda run -n smobench` for execution
- Reads from standard `Results/adata/` directory
- Saves to standard `Results/plot/` directory
- Supports all SMOBench clustering methods
- Handles spatial coordinate conventions