# /home/users/nus/e1503317/projects/dmeng/zhlin/SMOBench/Scripts/integration/SpaMosaic/run_SpaMosaic.py

import os
import sys
import argparse
import numpy as np
import scanpy as sc
import time
import pynvml

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.append(project_root)
from spamosaic.framework import SpaMosaic
from spamosaic.preprocessing2 import RNA_preprocess, ADT_preprocess, Epigenome_preprocess
import spamosaic.utils as utls
from Utils.SMOBench_clustering import universal_clustering

os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'


def main(args):
    print("Starting SpaMosaic integration...")

    # === 1. 读取数据 ===
    input_dict = {}

    # 支持多个模态：rna, adt, atac（可扩展）
    if args.rna_paths:
        rna_adatas = [sc.read_h5ad(p.strip()) for p in args.rna_paths.split(",") if p.strip()]
        input_dict['rna'] = rna_adatas
    else:
        raise ValueError("At least one RNA path must be provided via --rna_paths.")

    if args.adt_paths:
        adt_adatas = [sc.read_h5ad(p.strip()) for p in args.adt_paths.split(",") if p.strip()]
        input_dict['adt'] = adt_adatas

    if args.atac_paths:
        atac_adatas = [sc.read_h5ad(p.strip()) for p in args.atac_paths.split(",") if p.strip()]
        input_dict['atac'] = atac_adatas

    if not input_dict:
        raise ValueError("No valid modality data provided.")

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output will be saved to: {args.output_dir}")

    # === 2. 预处理 ===
    input_key = 'dimred_bc'

    if 'rna' in input_dict:
        RNA_preprocess(
            input_dict['rna'],
            batch_corr=args.batch_corr,
            favor='scanpy',
            n_hvg=args.n_hvg,
            batch_key=args.batch_key,
            key=input_key
        )

    if 'adt' in input_dict:
        ADT_preprocess(
            input_dict['adt'],
            batch_corr=args.batch_corr,
            batch_key=args.batch_key,
            key=input_key
        )

    if 'atac' in input_dict:
        Epigenome_preprocess(
            input_dict['atac'],
            batch_corr=args.batch_corr,
            batch_key=args.batch_key,
            key=input_key
        )

    # === 3. 构建模型 ===
    model = SpaMosaic(
        modBatch_dict=input_dict,
        input_key=input_key,
        batch_key=args.batch_key,
        intra_knn=args.intra_knn,
        inter_knn=args.inter_knn,
        w_g=args.w_g,
        seed=args.seed,
        device=args.device
    )

    # GPU 内存监控初始化
    if args.device.startswith('cuda'):
        pynvml.nvmlInit()
        handle = pynvml.nvmlDeviceGetHandleByIndex(0)
        info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        gpu_memory_before = info.used / (1024 * 1024)  # MB
    else:
        gpu_memory_before = None

    # === 4. 训练模型 ===
    start_time = time.time()
    model.train(net=args.net, lr=args.lr, T=args.T, n_epochs=args.n_epochs)
    running_time = time.time() - start_time

    if args.device.startswith('cuda'):
        info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        gpu_memory_after = info.used / (1024 * 1024)
    else:
        gpu_memory_after = None

    # === 5. 推断嵌入 ===
    ad_embs = model.infer_emb(
        input_dict,
        emb_key='emb',
        final_latent_key='merged_emb'
    )
    ad_mosaic = sc.concat(ad_embs)
    ad_mosaic = utls.get_umap(ad_mosaic, use_reps=['merged_emb'])

    # === 6. 保存嵌入 ===
    # 保存各模态各批次的嵌入
    for modality in model.mod_list:
        if modality in input_dict:
            for batch_idx in range(model.n_batches):
                adata = input_dict[modality][batch_idx]
                if adata is not None:
                    emb = adata.obsm['emb']
                    np.save(
                        os.path.join(args.output_dir, f'{modality}_batch{batch_idx}_embedding.npy'),
                        emb
                    )

    # 保存合并后的嵌入
    for batch_idx, adata in enumerate(ad_embs):
        np.save(
            os.path.join(args.output_dir, f'batch{batch_idx}_merged_embedding.npy'),
            adata.obsm['merged_emb']
        )

    # 保存 UMAP
    np.save(
        os.path.join(args.output_dir, 'umap.npy'),
        ad_mosaic.obsm['X_umap']
    )

    # === 7. 聚类 ===
    if args.clustering_method == 'leiden':
        resolution = args.leiden_resolution
        sc.tl.leiden(ad_mosaic, resolution=resolution, key_added='leiden')
        cluster_labels = ad_mosaic.obs['leiden'].values
        np.save(os.path.join(args.output_dir, 'leiden_y_pred_label.npy'), cluster_labels)
        print(f"Leiden clustering done with resolution={resolution}")

    elif args.clustering_method == 'mclust':
        n_clusters = args.mclust_n_clusters
        utls.clustering(ad_mosaic, n_cluster=n_clusters, used_obsm='merged_emb', algo='mclust', key='mclust')
        utls.split_adata_ob(ad_embs, ad_mosaic, 'obs', 'mclust')
        cluster_labels = ad_mosaic.obs['mclust'].values
        np.save(os.path.join(args.output_dir, 'mclust_y_pred_label.npy'), cluster_labels)
        print(f"mclust clustering done with n_clusters={n_clusters}")

    else:
        print("No clustering method selected.")

    # === 8. 可视化（可选）===
    if args.plot:
        # UMAP + spatial plots
        utls.plot_basis(ad_mosaic, basis='merged_emb_umap', color=['mclust'] if args.clustering_method == 'mclust' else ['leiden'])
        for ad in ad_embs:
            utls.plot_basis(ad, 'spatial', 'mclust' if args.clustering_method == 'mclust' else 'leiden', s=70)
        print("Plots generated.")

    # === 9. 保存运行信息 ===
    with open(os.path.join(args.output_dir, 'running_time_and_memory.txt'), 'w') as f:
        f.write(f"Running Time: {running_time:.2f} seconds\n")
        if gpu_memory_before is not None:
            f.write(f"GPU Memory Usage Before: {gpu_memory_before:.2f} MB\n")
            f.write(f"GPU Memory Usage After: {gpu_memory_after:.2f} MB\n")
        else:
            f.write("GPU Memory Usage: N/A (CPU mode)\n")

    print("SpaMosaic integration completed.")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SpaMosaic multi-omics integration")

    # 数据输入路径（支持多个批次，用逗号分隔）
    parser.add_argument('--rna_paths', type=str, required=True, help='Comma-separated paths to RNA adata files')
    parser.add_argument('--adt_paths', type=str, default='', help='Comma-separated paths to ADT adata files')
    parser.add_argument('--atac_paths', type=str, default='', help='Comma-separated paths to ATAC adata files')

    # 输出路径
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save results')

    # 预处理参数
    parser.add_argument('--batch_corr', action='store_true', help='Enable batch correction in preprocessing')
    parser.add_argument('--batch_key', type=str, default='src', help='Batch key in .obs')
    parser.add_argument('--n_hvg', type=int, default=5000, help='Number of highly variable genes')

    # 模型参数
    parser.add_argument('--intra_knn', type=int, default=2, help='Intra-modality KNN')
    parser.add_argument('--inter_knn', type=int, default=2, help='Inter-modality KNN')
    parser.add_argument('--w_g', type=float, default=0.8, help='Weight for graph regularization')
    parser.add_argument('--net', type=str, default='wlgcn', help='Network type for training')
    parser.add_argument('--lr', type=float, default=0.01, help='Learning rate')
    parser.add_argument('--T', type=float, default=0.01, help='Temperature for Gumbel-Softmax')
    parser.add_argument('--n_epochs', type=int, default=100, help='Number of training epochs')
    parser.add_argument('--seed', type=int, default=1234, help='Random seed')
    parser.add_argument('--device', type=str, default='cuda:0', help='Device: cuda:0 or cpu')

    # 聚类参数
    parser.add_argument('--clustering_method', type=str, choices=['none', 'leiden', 'mclust'], default='leiden', help='Clustering method')
    parser.add_argument('--leiden_resolution', type=float, default=1.0, help='Resolution for Leiden clustering')
    parser.add_argument('--mclust_n_clusters', type=int, default=6, help='Number of clusters for mclust')

    # 其他
    parser.add_argument('--plot', action='store_true', help='Generate and save plots')

    args = parser.parse_args()
    main(args)