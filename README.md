## Framework Overview

**Integration Tasks:**
- **Vertical Integration**: Cross-modality integration within the same sample (RNA+ADT, RNA+ATAC)
- **Horizontal Integration**: Cross-sample integration with batch effect removal
- **Mosaic Integration**: Mixed modality and batch integration

**Evaluation Dimensions:**
- **Spatial Coherence (SC)**: Spatial clustering quality
- **Biological Conservation (BioC)**: Clustering accuracy and biological validity
- **Batch Effect Removal (BER)**: Cross-batch mixing metrics (horizontal/mosaic only)

**Methods Evaluated:**
CANDIES, COSMOS, PRAGA, PRESENT, SpaMV, SpaMosaic, SpatialGlue, SpaMultiVAE

**Clustering Methods:**
Leiden, Louvain, K-means, Mclust

## Project Structure

- **Dataset/**: Multi-modal spatial omics datasets (withGT/woGT) - *Download from Google Drive*
- **Methods/**: Integration method implementations
- **Scripts/**: Execution workflows for integration and evaluation
- **Eval/**: Evaluation framework and metrics calculation
- **Results/**: Integration outputs and evaluation results - *Download from Google Drive*
- **Utils/**: Shared utilities and clustering interface

## Quick Start

1. **Download Data**: Get datasets and results from Google Drive link below
2. **Integration**: Run methods using scripts in `Scripts/`
3. **Evaluation**: Use evaluation scripts in `Eval/`

## Data Access

**Note**: Due to large file sizes, datasets and integration results are not included in this repository.

- **Complete Data Package**: https://drive.google.com/drive/u/1/folders/11zYh27BK9QuqU7zObApCYSzSEMqHS0G6
  - Multi-modal spatial omics datasets (Dataset folder)
  - Integration results in AnnData format (Results folder)
  - All processed evaluation outputs

