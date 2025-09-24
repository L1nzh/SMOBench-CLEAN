# Results Directory

Integration outputs and evaluation results from spatial multi-omics benchmarking.

**Note**: This directory is not included in the GitHub repository due to large file sizes. Download from Google Drive: https://drive.google.com/drive/u/1/folders/11zYh27BK9QuqU7zObApCYSzSEMqHS0G6

## Directory Structure

### adata/
**Vertical Integration Results:**
- `vertical_integration/`: AnnData files with integration embeddings
- Format: `{Method}_{Dataset}_{Sample}.h5ad`
- Contains: Integration embeddings, clustering results, spatial coordinates

**Horizontal Integration Results:**
- `horizontal_integration/`: Cross-sample integration outputs  
- Format: `{Method}_{Dataset}_horizontal.h5ad`
- Contains: Batch-corrected embeddings, multi-sample clustering

### evaluation/
**Vertical Integration Evaluation:**
- `vertical_integration/`: Individual method evaluation results
- `vertical_integration/final_results/`: Aggregated summary tables

**Horizontal Integration Evaluation:**  
- `horizontal_integration/`: Individual method evaluation results
- `horizontal_integration/final_results/`: Aggregated summary tables

### plots/
Visualization outputs including UMAP plots, spatial plots, and performance comparisons.

## File Formats

**Integration Results (AnnData):**
- `adata.obsm[METHOD]`: Integration embeddings
- `adata.obs[clustering]`: Clustering assignments
- `adata.obsm['spatial']`: Spatial coordinates

**Evaluation Results (CSV):**
- Individual metrics per method-dataset-clustering combination
- Summary tables with overall performance scores

## Data Access

To use this framework:
1. Download the complete Results folder from Google Drive
2. Place it in the SMOBench-CLEAN root directory
3. The evaluation scripts will automatically find the results