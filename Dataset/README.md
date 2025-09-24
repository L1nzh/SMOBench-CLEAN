# Dataset Directory

Multi-modal spatial omics datasets for benchmarking spatial integration methods.

## Dataset Structure

### With Ground Truth (withGT/)
Datasets with cell type annotations for supervised evaluation.

**RNA_ADT (Protein)**
- **Human_Lymph_Nodes** (HLN): 2 samples with spatial protein data
- **Human_Tonsils** (HT): 3 samples with spatial protein data  
- **Mouse_Thymus**: 4 samples with spatial protein data
- **Mouse_Spleen**: 2 samples with spatial protein data

**RNA_ATAC (Chromatin)**
- **Mouse_Embryos_S1** (MISAR_S1): 4 timepoints with chromatin accessibility
- **Mouse_Embryos_S2** (MISAR_S2): 4 timepoints with chromatin accessibility

### Without Ground Truth (woGT/)
Datasets without cell type annotations for unsupervised evaluation.

**RNA_ADT (Protein)**
- **Mouse_Thymus**: 4 samples
- **Mouse_Spleen**: 2 samples

**RNA_ATAC (Chromatin)**  
- **Mouse_Brain**: 4 chromatin modalities (ATAC, H3K27ac, H3K27me3, H3K4me3)

## Data Format

Each dataset contains:
- **adata_RNA.h5ad**: RNA expression matrix with spatial coordinates
- **adata_ADT.h5ad**: Protein expression matrix (for ADT datasets)
- **adata_ATAC.h5ad**: Chromatin accessibility matrix (for ATAC datasets)
- **Spatial_Label**: Ground truth cell type annotations (withGT only)

## Download

Complete datasets available at:
https://drive.google.com/drive/u/1/folders/11zYh27BK9QuqU7zObApCYSzSEMqHS0G6