#!/usr/bin/env python3

import os
import torch
import pandas as pd
import scanpy as sc
import argparse
import time
import sys
import re
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Add project root to path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.append(project_root)

# Add SpaMV to path
spamv_path = os.path.join(project_root, "Methods/SpaMV")
sys.path.append(spamv_path)

from SpaMV.spamv import SpaMV
from SpaMV.utils import preprocess_dc, clustering
from Utils.SMOBench_clustering import universal_clustering


def parse_dataset_info(args):
    """
    从 RNA_path 或 save_path 中提取 dataset_name 和 subset_name
    支持两种模式：
    1. 手动指定 --dataset Human_Lymph_Nodes/A1
    2. 自动从路径中提取
    """
    if hasattr(args, 'dataset') and args.dataset:
        parts = args.dataset.strip('/').split('/')
        if len(parts) == 2:
            return parts[0], parts[1]
        elif len(parts) == 1:
            return parts[0], "Unknown"
    
    # 自动解析 RNA_path - support both absolute and relative paths
    match = re.search(r'(SMOBench_Data|Dataset)/([^/]+)/([^/]+)/adata_RNA\.h5ad', args.RNA_path)
    if match:
        return match.group(2), match.group(3)
    return "Unknown", "Unknown"


def main(args):
    print(f"Starting SpaMV integration...")
    
    # === 读取数据 ===
    adata_omics1 = sc.read_h5ad(args.RNA_path)
    adata_omics1.var_names_make_unique()
    print(f"Loaded RNA data: {adata_omics1.shape}")

    if args.ADT_path:
        adata_omics2 = sc.read_h5ad(args.ADT_path)
        modality = 'ADT'
        omics_names = ['Transcriptome', 'Proteome']
        adata_omics2.var_names_make_unique()
        print(f"Loaded ADT data: {adata_omics2.shape}")
        # Use scaling for ADT data as per tutorial
        scale = True
    elif args.ATAC_path:
        adata_omics2 = sc.read_h5ad(args.ATAC_path)
        modality = 'ATAC'
        omics_names = ['Transcriptome', 'Epigenome']
        adata_omics2.var_names_make_unique()
        print(f"Loaded ATAC data: {adata_omics2.shape}")
        # No scaling for ATAC data as per tutorial
        scale = False
    else:
        raise ValueError("Either ADT_path or ATAC_path must be provided.")

    # === 使用 SpaMV 官方预处理（修改避免seurat_v3问题） ===
    print("Preprocessing with modified SpaMV preprocess_dc...")
    datasets = [adata_omics1, adata_omics2]
    
    # Use custom preprocessing to avoid seurat_v3 issue
    try:
        datasets = preprocess_dc(datasets, omics_names, scale=scale)
    except Exception as e:
        print(f"Official preprocess_dc failed: {e}")
        print("Using custom preprocessing...")
        # Custom preprocessing following SpaMV patterns
        from SpaMV.utils import clr_normalize_each_cell, lsi
        
        # Process RNA data
        sc.pp.filter_genes(datasets[0], min_cells=10)
        sc.pp.filter_cells(datasets[0], min_genes=200)
        sc.pp.normalize_total(datasets[0], target_sum=1e4)
        sc.pp.log1p(datasets[0])
        sc.pp.highly_variable_genes(datasets[0], n_top_genes=3000)  # Use default flavor after normalization
        datasets[0] = datasets[0][:, datasets[0].var.highly_variable]
        if scale:
            sc.pp.scale(datasets[0])
        sc.pp.pca(datasets[0], n_comps=min(50, datasets[0].n_vars - 1))
        datasets[0].obsm['embedding'] = datasets[0].obsm['X_pca']
        
        # Process second modality
        if modality == 'ADT':
            datasets[1] = clr_normalize_each_cell(datasets[1])
            if scale:
                sc.pp.scale(datasets[1])
            sc.pp.pca(datasets[1], n_comps=min(50, datasets[1].n_vars - 1))
            datasets[1].obsm['embedding'] = datasets[1].obsm['X_pca']
        elif modality == 'ATAC':
            lsi(datasets[1], use_highly_variable=False, n_components=min(50, datasets[1].n_vars - 1) + 1)
            datasets[1].obsm['embedding'] = datasets[1].obsm['X_lsi']

    # === 训练 SpaMV 模型 ===
    print("Training SpaMV model...")
    
    # 设置随机种子
    if hasattr(args, 'seed'):
        torch.manual_seed(args.seed)
        
    model = SpaMV(datasets, interpretable=False, max_epochs_stage1=args.max_epochs)
    
    start_time = time.time()
    model.train()
    end_time = time.time()
    print(f'Training time: {end_time - start_time:.2f} seconds')

    # === 获取嵌入 ===
    embedding = model.get_embedding()
    print(f"SpaMV embedding shape: {embedding.shape}")

    # === 构建结果 AnnData ===
    adata = datasets[0].copy()  # Use preprocessed data
    adata.obsm['SpaMV'] = embedding

    # === 解析数据集信息 ===
    dataset_name, subset_name = parse_dataset_info(args)
    print(f"Detected dataset: {dataset_name}, subset: {subset_name}")

    # === 图像保存路径 ===
    plot_base_dir = "Results/plot"  # Use relative path for portability
    method_name = args.method if args.method else "SpaMV"
    plot_dir = os.path.join(plot_base_dir, method_name, dataset_name, subset_name)
    os.makedirs(plot_dir, exist_ok=True)
    print(f"Plot images will be saved to: {plot_dir}")

    # === 聚类与可视化 ===
    # 使用universal_clustering支持多种聚类方法（mclust, leiden, louvain等）
    # 可以轻松扩展到其他聚类方法
    clustering_methods = ['mclust']  # mclust is working properly
    
    # 翻转Y坐标用于可视化
    if 'spatial' in adata.obsm:
        adata.obsm['spatial'][:, 1] = -1 * adata.obsm['spatial'][:, 1]

    # 计算UMAP（只需计算一次）
    print("Computing UMAP...")
    sc.pp.neighbors(adata, use_rep='SpaMV', n_neighbors=30)
    sc.tl.umap(adata)

    for method in clustering_methods:
        print(f"Performing {method} clustering with {args.cluster_nums} clusters...")
        adata = universal_clustering(
            adata,
            n_clusters=args.cluster_nums,
            used_obsm='SpaMV',
            method=method,
            key=method,
            use_pca=False
        )

        # 创建可视化
        fig, ax_list = plt.subplots(1, 2, figsize=(12, 5))
        
        sc.pl.umap(adata, color=method, ax=ax_list[0], title=f'{method_name}-{method} (UMAP)', 
                   s=20, show=False, legend_loc='right margin')
        
        if 'spatial' in adata.obsm:
            sc.pl.embedding(adata, basis='spatial', color=method, ax=ax_list[1], 
                           title=f'{method_name}-{method} (Spatial)', s=20, show=False, 
                           legend_loc='right margin')
        else:
            # 如果没有空间信息，显示第二个UMAP
            sc.pl.umap(adata, color=method, ax=ax_list[1], title=f'{method_name}-{method}', 
                       s=20, show=False, legend_loc='right margin')

        plt.tight_layout()
        plt.savefig(
            os.path.join(plot_dir, f'clustering_{method}_umap_spatial.png'),
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()
        print(f"Saved {method} clustering plot to: {plot_dir}/clustering_{method}_umap_spatial.png")

    # === 保存 AnnData ===
    save_dir = os.path.dirname(args.save_path)
    os.makedirs(save_dir, exist_ok=True)
    adata.write(args.save_path)
    print(f'Results saved to: {args.save_path}')
    print(f'Final AnnData shape: {adata.shape}')
    print("SpaMV integration completed successfully!")


if __name__ == "__main__":
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1' 
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    print("Starting SpaMV integration...")
    parser = argparse.ArgumentParser(description='Run SpaMV integration')
    
    # Data paths
    parser.add_argument('--RNA_path', type=str, required=True, help='Path to RNA adata')
    parser.add_argument('--ADT_path', type=str, default='', help='Path to ADT adata')
    parser.add_argument('--ATAC_path', type=str, default='', help='Path to ATAC adata')
    parser.add_argument('--save_path', type=str, required=True, help='Path to save integrated adata')
    
    # Model parameters
    parser.add_argument('--seed', type=int, default=2024, help='Random seed')
    parser.add_argument('--max_epochs', type=int, default=200, help='Maximum number of training epochs')
    
    # Clustering parameters
    parser.add_argument('--cluster_nums', type=int, required=True, help='Number of clusters')
    
    # Optional parameters
    parser.add_argument('--method', type=str, default='SpaMV', help='Method name for plotting')
    parser.add_argument('--dataset', type=str, default='', help='Dataset name, e.g. Human_Lymph_Nodes/A1')

    args = parser.parse_args()
    main(args)