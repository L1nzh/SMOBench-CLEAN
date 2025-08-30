# SMOBench Development Standards

## Mandatory Files
These files MUST be present and referenced in every conversation:
- `CLAUDE.md` - This file with development standards
- `docs/brief.md` - Project brief and research objectives
- `docs/workflow.md` - Research workflow and experimental processes
- `docs/status.md` - Current research status and progress
- `docs/setup.md` - Environment setup and installation guide

## Core Constraints

### Technology Stack
- **Language**: Python 3.9+ (conda environment: smobench)
- **ML Framework**: PyTorch with CUDA support, scanpy, AnnData
- **Spatial Omics**: SpatialGlue, SpaMosaic, PRAGA, COSMOS, PRESENT
- **Data Processing**: NumPy, pandas, h5py, scanpy
- **Clustering**: scikit-learn, leiden, louvain, mclust (R interface)

### Code Standards
- Use conda run -n smobench for all Python executions
- Write reproducible research code with fixed random seeds (default: 2024)
- Implement universal clustering interface for method comparison
- Maintain compatibility with AnnData/h5ad format

### Research Requirements
- Support three integration tasks: vertical, horizontal, mosaic
- Evaluate with/without ground truth scenarios
- Test multiple clustering methods (leiden, louvain, kmeans, mclust)
- Generate comprehensive evaluation metrics and visualizations

### Project Structure
```
SMOBench-CLEAN/
├── Dataset/                    # Organized datasets
│   ├── data_info/             # Dataset information files
│   ├── withGT/                # Datasets with ground truth
│   │   ├── 3M_Simulation/     # Multi-modal simulation data
│   │   ├── RNA_ADT/           # RNA + protein datasets
│   │   └── RNA_ATAC/          # RNA + chromatin datasets
│   └── woGT/                  # Datasets without ground truth
│       ├── RNA_ADT/           # Protein datasets for evaluation
│       └── RNA_ATAC/          # Chromatin datasets for evaluation
├── Methods/                   # Integration methods
│   ├── SpatialGlue/          # Graph-based integration
│   ├── SpaMosaic/            # Mosaic integration
│   ├── PRAGA/                # Adversarial integration
│   ├── COSMOS/               # Multi-modal integration
│   └── PRESENT/              # Spatial-aware integration
├── Scripts/                  # Organized execution scripts
│   ├── integration/          # Method execution scripts
│   ├── evaluation/           # Evaluation and metrics
│   └── visualization/        # Plotting and analysis
├── Utils/                    # Utilities and common functions
│   ├── SMOBench_clustering.py # Universal clustering interface
│   └── evaluation_utils.py   # Evaluation metrics
├── Results/                  # Experimental outputs
│   ├── adata/                # Processed AnnData results
│   └── plots/                # Visualization outputs
└── docs/                     # Documentation
```

### GPU Management Protocol
When encountering GPU memory issues:
1. **Check Availability**: `nvidia-smi` to identify free GPUs
2. **Switch GPU**: Use `CUDA_VISIBLE_DEVICES=N` where N is the free GPU index
3. **Clear Cache**: `python -c "import torch; torch.cuda.empty_cache()"`
4. **Monitor Usage**: Use `nvidia-smi` continuously with `watch -n 1 nvidia-smi`

### Integration Task Definitions
1. **Vertical Integration**: Same batch, cross-modality (e.g., RNA + ADT from same sample)
2. **Horizontal Integration**: Same modality, cross-batch (e.g., RNA from multiple batches)  
3. **Mosaic Integration**: Mixed modality and batch integration

### Data Processing Pipeline
1. **Load Data**: Use scanpy.read_h5ad() for AnnData objects
2. **Preprocess**: Apply var_names_make_unique() for duplicate genes
3. **Integration**: Run spatial multi-omics integration methods
4. **Clustering**: Apply universal clustering with multiple algorithms
5. **Evaluation**: Calculate metrics with/without ground truth
6. **Visualization**: Generate UMAP and spatial plots

### Clustering Standards
- **Methods**: leiden, louvain, kmeans, mclust
- **Resolution Search**: Automatic resolution tuning for leiden/louvain
- **Input**: Use integration embeddings from adata.obsm
- **Output**: Store results in adata.obs with method-specific keys
- **Interface**: Use Utils/SMOBench_clustering.py universal function

### Evaluation Metrics
#### With Ground Truth
- **ARI**: Adjusted Rand Index (clustering accuracy)
- **NMI**: Normalized Mutual Information 
- **FMI**: Fowlkes-Mallows Index
- **Silhouette Score**: Cluster cohesion quality

#### Without Ground Truth
- **Silhouette Score**: Unsupervised cluster quality
- **Moran's I**: Spatial autocorrelation of clusters
- **Spatial Consistency**: Neighborhood cluster agreement

### Development Environment
- **Conda Environment**: `smobench` with spatial omics dependencies
- **GPU**: NVIDIA with sufficient VRAM for large datasets
- **Storage**: Adequate space for h5ad files (~50GB+ for full dataset)
- **R Interface**: rpy2 for mclust clustering integration

### File Naming Conventions
- **Results**: `{Method}_{Dataset}_{Sample}.h5ad`
- **Plots**: `clustering_{method}_{plot_type}.png`
- **Metrics**: `{Method}_{Dataset}_metrics_{gt_status}.csv`
- **Scripts**: `run_{method}.py` and `run_{method}.sh`

### Script Execution Standards
- Use PBS template system to reduce redundancy
- Generate batch scripts from CSV parameter files
- Include comprehensive error handling and logging
- Support both single runs and batch processing
- Use relative paths with proper base directory handling

### Quality Assurance
1. **Reproducibility**: Fixed random seeds across all methods
2. **Validation**: Test on known datasets before new experiments
3. **Documentation**: Comment complex spatial processing steps
4. **Version Control**: Commit major changes with descriptive messages

### Method Integration Requirements
Each spatial integration method must:
1. Accept RNA and secondary modality (ADT/ATAC) as input
2. Output integrated embeddings in AnnData.obsm format
3. Preserve spatial coordinates and metadata
4. Support both 10x and spatial transcriptomics data formats
5. Handle variable gene names and feature preprocessing

## Mandatory Documentation
Each development session must reference these files and update them as needed. All experimental results must be logged with:
- Method parameters and hyperparameters
- Dataset characteristics and preprocessing steps
- Evaluation metrics for all clustering methods
- Computational resources and runtime information
- Reproducible execution commands with conda environment