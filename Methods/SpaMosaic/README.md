# SpaMosaic Integration Method

This directory contains the SpaMosaic integration method implementation for SMOBench.

## Directory Structure

```
Methods/SpaMosaic/
├── README.md          # This file
├── spamosaic/         # SpaMosaic core implementation
└── backup/            # Original dataset-specific scripts (archived)
```

## New Unified Implementation

The SpaMosaic method has been reorganized to use a unified script that supports all integration types:

- **Location**: `Scripts/integration/SpaMosaic/`
- **Main Script**: `run_spamosaic.py`
- **Execution Script**: `run.sh`

## Usage

### Basic Usage

```bash
# Navigate to the script directory
cd Scripts/integration/SpaMosaic/

# Vertical integration (RNA+ADT) for HLN dataset
./run.sh HLN vertical "rna adt" A1 D1

# Horizontal integration (RNA only) for Mouse_Brain
./run.sh Mouse_Brain horizontal rna sample1 sample2 sample3

# Mosaic integration (RNA+ATAC) for MISAR
./run.sh MISAR mosaic "rna atac" E11 E13 E15 E18
```

### Advanced Usage

```bash
# Custom parameters
./run.sh -g 1 -c 10 --epochs 200 --lr 0.005 HLN vertical "rna adt" A1 D1

# Custom data and output directories
./run.sh -d /path/to/data -o /path/to/output HLN vertical "rna adt" A1

# Dry run to see the command
./run.sh --dry-run HLN vertical "rna adt" A1 D1
```

### Help

```bash
./run.sh --help
```

## Integration Types Supported

1. **Vertical Integration**: Same batch, cross-modality (e.g., RNA + ADT from same sample)
2. **Horizontal Integration**: Same modality, cross-batch (e.g., RNA from multiple samples)  
3. **Mosaic Integration**: Mixed modality and batch integration

## Features

- **Unified Interface**: Single script handles all integration types
- **SMOBench Compatible**: Follows SMOBench standards and uses universal clustering
- **Flexible Parameters**: Customizable GPU, learning rate, epochs, clustering methods
- **Comprehensive Output**: Embeddings, clustering results, metrics, and integrated AnnData
- **Error Handling**: Robust error checking and informative messages
- **GPU Monitoring**: Tracks GPU memory usage during training

## Requirements

- Conda environment: `smobench`
- SpaMosaic package installed
- SMOBench Utils for universal clustering
- NVIDIA GPU for CUDA acceleration (optional)

## Output Structure

```
Results/adata/SpaMosaic/DATASET/
├── SpaMosaic_DATASET_integrated.h5ad   # Integrated AnnData object
├── umap.npy                            # UMAP coordinates
├── metrics.txt                         # Training metrics and parameters
├── embeddings/                         # Individual and merged embeddings
│   ├── rna_batch0_embedding.npy
│   ├── adt_batch0_embedding.npy
│   └── batch0_merged_embedding.npy
└── clustering/                         # Clustering results
    ├── leiden_labels.npy
    ├── louvain_labels.npy
    ├── kmeans_labels.npy
    └── mclust_labels.npy
```

## Migration from Old Scripts

The original dataset-specific scripts have been moved to the `backup/` directory. The new unified implementation provides:

1. **Better Maintainability**: Single script instead of multiple dataset-specific scripts
2. **Consistent Interface**: Standardized parameters and outputs
3. **Enhanced Features**: GPU monitoring, comprehensive clustering, better error handling
4. **SMOBench Integration**: Compatible with SMOBench evaluation framework

## Configuration

The script automatically detects and uses appropriate preprocessing based on:
- **RNA**: Highly variable gene selection, log normalization, PCA
- **ADT**: Centered log ratio normalization, PCA  
- **ATAC**: Peak alignment, dimensionality reduction

Batch correction is applied automatically for multi-sample integration tasks.