# SMOBench: Spatial Multi-Omics Integration Benchmark

A comprehensive benchmarking framework for evaluating spatial multi-omics integration methods.

## Overview

SMOBench provides standardized datasets, evaluation metrics, and integration pipelines for benchmarking spatial multi-omics integration methods. It supports both vertical integration (cross-modality, same sample) and horizontal integration (same modality, cross-sample) tasks.

## Dataset
Complete benchmark datasets available at: https://drive.google.com/file/d/1XPmGicNOeaKjMnMGqvMPD3HLIBy8ostX/view?usp=drive_link

## Quick Start

### Environment Setup
```bash
conda create -n smobench python=3.9
conda activate smobench
pip install scanpy pandas numpy scikit-learn torch
```

### Running Integration Methods

#### Vertical Integration (Single Sample, Multi-Modal)
```bash
# SpaMultiVAE
conda run -n smobench python Scripts/vertical_integration/SpaMultiVAE/run_spamultivae.py

# SpaMosaic  
conda run -n smobench python Scripts/vertical_integration/SpaMosaic/run_spamosaic.py

# SpatialGlue
conda run -n smobench python Scripts/vertical_integration/SpatialGlue/run_spatialglue.py

# PRAGA
conda run -n smobench python Scripts/vertical_integration/PRAGA/run_praga.py
```

#### Horizontal Integration (Multi-Sample Fusion)
```bash
# Run horizontal integration for specific method
conda run -n smobench python Scripts/horizontal_integration/[METHOD]/run_[method].py
```

### Evaluation
```bash
# Evaluate specific results
conda run -n smobench python Eval/eval_adata.py --method SpaMultiVAE --dataset HLN

# Generate comprehensive results
conda run -n smobench python Eval/generate_final_results.py

# Validate all results
python test.py
```

## Methods

### Implemented Methods

| Method | Type | Integration | Clustering | Status |
|--------|------|-------------|------------|--------|
| **SpaMultiVAE** | VAE | Vertical + Horizontal | ✓ | Complete |
| **SpaMosaic** | Mosaic | Vertical + Horizontal | ✓ | Complete |
| **SpatialGlue** | Graph | Vertical + Horizontal | ✓ | Complete |
| **PRAGA** | GAN | Vertical + Horizontal | ✓ | Complete |

## Datasets

### With Ground Truth (withGT)
- **HLN**: 2 samples (A1, D1) - RNA + ADT
- **HT**: 3 samples (S1, S2, S3) - RNA + ADT  
- **MISAR_S1**: 4 timepoints (E11, E13, E15, E18) - RNA + ATAC
- **MISAR_S2**: 4 timepoints (E11, E13, E15, E18) - RNA + ATAC

### Without Ground Truth (woGT)
- **Mouse_Thymus**: 4 samples - RNA + ADT
- **Mouse_Spleen**: 2 samples - RNA + ADT
- **Mouse_Brain**: 4 modalities (ATAC, H3K27ac, H3K27me3, H3K4me3) - RNA + Chromatin

## Project Structure

```
SMOBench-CLEAN/
├── Dataset/                     # Benchmark datasets
│   ├── withGT/                  # Datasets with ground truth
│   └── woGT/                    # Datasets without ground truth
├── Methods/                     # Integration method implementations
│   ├── SpaMultiVAE/            # VAE-based integration
│   ├── SpaMosaic/              # Mosaic integration
│   ├── SpatialGlue/            # Graph-based integration
│   └── PRAGA/                  # Adversarial integration
├── Scripts/                     # Execution scripts
│   ├── vertical_integration/    # Vertical integration workflows
│   ├── horizontal_integration/  # Horizontal integration workflows
│   └── data_preparation/        # Data preprocessing
├── Utils/                       # Utility functions
│   ├── SMOBench_clustering.py  # Universal clustering interface
│   └── evaluation_utils.py     # Evaluation metrics
├── Results/                     # Integration results
│   ├── adata/                  # AnnData format results
│   └── plots/                  # Visualization outputs
├── Eval/                       # Evaluation framework
│   ├── src/                    # Evaluation source code
│   └── generate_final_results.py # Results aggregation
└── test.py                     # Comprehensive validation
```

## Evaluation Metrics

### Spatial Coherence (SC)
- **Moran's I**: Spatial autocorrelation
- **Geary's C**: Spatial association

### Biological Conservation (BioC)
- **ARI/NMI**: Clustering accuracy (withGT)
- **Silhouette**: Cluster quality
- **Davies-Bouldin**: Cluster separation

## Results Format

All results are stored in AnnData format with standardized structure:

```python
adata.obsm['METHOD_NAME']    # Integration embeddings
adata.obsm['X_umap']         # UMAP coordinates
adata.obsm['spatial']        # Spatial coordinates
adata.obs['mclust']          # mclust clustering
adata.obs['leiden']          # Leiden clustering  
adata.obs['louvain']         # Louvain clustering
adata.obs['kmeans']          # K-means clustering
adata.uns['method']          # Method metadata
```

## Citation

If you use SMOBench in your research, please cite our work.

