# SMOBench Dataset

This folder contains the spatial multi-omics datasets used in the SMOBench benchmark framework.

## Download Data

The complete dataset can be downloaded from Google Drive:

**[Download SMOBench Dataset (SMOBench_Data.tar.gz)](https://drive.google.com/file/d/1XPmGicNOeaKjMnMGqvMPD3HLIBy8ostX/view?usp=drive_link)**

After downloading, extract the data using:
```bash
tar -xzf SMOBench_Data.tar.gz
```

## Dataset Structure

```
Dataset/
├── data_info/                          # Dataset information files
├── withGT/                            # Datasets with ground truth spatial labels
│   ├── RNA_ADT/                      # RNA + Protein datasets
│   │   ├── Human_Lymph_Nodes/        # Human lymph node samples
│   │   │   ├── A1/                   # Sample A1
│   │   │   │   ├── adata_RNA.h5ad    # RNA expression data
│   │   │   │   └── adata_ADT.h5ad    # Protein expression data
│   │   │   └── D1/                   # Sample D1
│   │   └── Human_Tonsils/            # Human tonsil samples
│   │       ├── S1/, S2/, S3/         # Multiple samples
│   ├── RNA_ATAC/                     # RNA + Chromatin accessibility datasets
│   │   ├── Mouse_Embryos_S1/         # Mouse embryo samples S1
│   │   │   ├── E11/, E13/, E15/, E18/ # Different developmental stages
│   │   └── Mouse_Embryos_S2/         # Mouse embryo samples S2
│   └── 3M_Simulation/                # Simulated 3 million cell dataset
│       ├── adata_RNA.h5ad           # Simulated RNA data
│       └── adata_ADT.h5ad           # Simulated protein data
└── woGT/                            # Datasets without ground truth (for evaluation)
    ├── RNA_ADT/                     # RNA + Protein datasets
    │   ├── Mouse_Thymus/            # Mouse thymus samples
    │   │   ├── Mouse_Thymus1/, Mouse_Thymus2/, Mouse_Thymus3/, Mouse_Thymus4/
    │   └── Mouse_Spleen/            # Mouse spleen samples
    │       ├── Mouse_Spleen1/, Mouse_Spleen2/
    └── RNA_ATAC/                    # RNA + Chromatin accessibility datasets
        └── Mouse_Brain/             # Mouse brain samples
            ├── Mouse_Brain_ATAC/
            ├── Mouse_Brain_H3K4me3/
            ├── Mouse_Brain_H3K27ac/
            └── Mouse_Brain_H3K27me3/
```

## Data Format

All datasets are stored in **AnnData format** (`.h5ad` files) compatible with scanpy and spatial omics analysis tools.

### Key Features:

- **withGT datasets**: Contains ground truth spatial domain labels in `adata.obs['spatial_label']`
- **woGT datasets**: Used for benchmarking without ground truth labels
- **Spatial coordinates**: Available in `adata.obsm['spatial']`
- **Multi-modal data**: Separate files for different omics modalities (RNA, ADT, ATAC)

### Dataset Types:

1. **10x Visium**: Human lymph nodes and tonsils (RNA + protein)
2. **Stereo-CITE-seq**: Mouse thymus (RNA + protein)  
3. **SPOTS**: Mouse spleen (RNA + protein)
4. **MISAR**: Mouse embryos (RNA + chromatin accessibility)
5. **Spatial-epigenome-transcriptome**: Mouse brain (RNA + chromatin accessibility)
6. **Simulation**: Large-scale simulated dataset (3M cells)

## Usage

Each dataset includes:
- Gene/protein/peak expression matrices
- Spatial coordinates
- Cell type annotations (where available)
- Quality metrics and metadata

The datasets support various spatial multi-omics integration tasks:
- **Vertical integration**: Cross-modality integration (RNA + ADT/ATAC)
- **Horizontal integration**: Cross-batch integration  
- **Mosaic integration**: Mixed modality and batch integration

## Citation

If you use this dataset, please cite the SMOBench paper and the original data sources as specified in the documentation.