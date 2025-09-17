#!/usr/bin/env python3
"""
转换COSMOS结果到SMOBench标准格式

使用说明:
1. 确保在smobench conda环境中运行: conda activate smobench
2. 运行脚本: python convert_cosmos_results.py
3. 结果将保存在 Results/adata/ 目录下，文件名以 t_ 开头

功能说明:
- 读取yjxiao/COSMOS/results下的所有结果
- 提取embedding并计算UMAP坐标
- 运行4种聚类方法 (mclust, leiden, louvain, kmeans)
- 保存为标准SMOBench格式的h5ad文件
"""

import os
import sys
import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add project root to path
sys.path.append('/home/zhenghong/SMOBench-CLEAN')
from Utils.SMOBench_clustering import universal_clustering

def load_original_data(dataset_type, sample_name):
    """加载原始数据来获取metadata"""
    base_path = "/home/zhenghong/SMOBench-CLEAN/Dataset"
    
    # 映射COSMOS的数据集名称到我们的路径
    mapping = {
        'HLN': {
            'A1': f"{base_path}/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_RNA.h5ad",
            'D1': f"{base_path}/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_RNA.h5ad"
        },
        'HT': {
            'slice1': f"{base_path}/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad",
            'slice2': f"{base_path}/withGT/RNA_ADT/Human_Tonsils/S2/adata_RNA.h5ad", 
            'slice3': f"{base_path}/withGT/RNA_ADT/Human_Tonsils/S3/adata_RNA.h5ad"
        },
        'MISAR_S1': {
            'E11': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S1/E11/adata_RNA.h5ad",
            'E13': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S1/E13/adata_RNA.h5ad",
            'E15': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S1/E15/adata_RNA.h5ad",
            'E18': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S1/E18/adata_RNA.h5ad"
        },
        'MISAR_S2': {
            'E11': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S2/E11/adata_RNA.h5ad",
            'E13': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S2/E13/adata_RNA.h5ad",
            'E15': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S2/E15/adata_RNA.h5ad",
            'E18': f"{base_path}/withGT/RNA_ATAC/Mouse_Embryos_S2/E18/adata_RNA.h5ad"
        },
        'Mouse_Thymus': {
            'Dataset3_Mouse_Thymus1': f"{base_path}/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus1/adata_RNA.h5ad",
            'Dataset4_Mouse_Thymus2': f"{base_path}/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus2/adata_RNA.h5ad",
            'Dataset5_Mouse_Thymus3': f"{base_path}/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus3/adata_RNA.h5ad",
            'Dataset6_Mouse_Thymus4': f"{base_path}/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus4/adata_RNA.h5ad"
        },
        'Mouse_slpeen': {
            'Dataset1_Mouse_Spleen1': f"{base_path}/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen1/adata_RNA.h5ad",
            'Dataset2_Mouse_Spleen2': f"{base_path}/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen2/adata_RNA.h5ad"
        },
        'Mouse_Brain': {
            'Dataset7_Mouse_Brain_ATAC': f"{base_path}/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_RNA.h5ad",
            'Dataset8_Mouse_Brain_H3K4me3': f"{base_path}/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K4me3/adata_RNA.h5ad",
            'Dataset9_Mouse_Brain_H3K27ac': f"{base_path}/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K27ac/adata_RNA.h5ad",
            'Dataset10_Mouse_Brain_H3K27me3': f"{base_path}/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K27me3/adata_RNA.h5ad"
        }
    }
    
    if dataset_type in mapping and sample_name in mapping[dataset_type]:
        original_path = mapping[dataset_type][sample_name]
        if os.path.exists(original_path):
            return sc.read_h5ad(original_path)
    
    return None

def get_cluster_numbers():
    """获取每个数据集的正确聚类数量"""
    cluster_nums = {
        # withGT datasets - based on ground truth
        'HLN_A1': 10, 'HLN_D1': 11,
        'HT_slice1': 4, 'HT_slice2': 5, 'HT_slice3': 5,
        'MISAR_S1_E11': 8, 'MISAR_S1_E13': 12, 'MISAR_S1_E15': 12, 'MISAR_S1_E18': 14,
        'MISAR_S2_E11': 13, 'MISAR_S2_E13': 14, 'MISAR_S2_E15': 15, 'MISAR_S2_E18': 16,
        
        # woGT datasets - preset numbers
        'Mouse_Thymus': 8, 'Mouse_Spleen': 5, 'Mouse_Brain': 18
    }
    return cluster_nums

def convert_cosmos_sample(sample_path, output_dir):
    """转换单个COSMOS样本到SMOBench格式"""
    
    sample_dir = Path(sample_path)
    dataset_type = sample_dir.parent.name
    sample_name = sample_dir.name
    
    # 跳过Fusion数据，之后统一处理
    if 'Fusion' in sample_name:
        print(f"\n⏭️  跳过Fusion数据: {dataset_type}/{sample_name}")
        return None
    
    print(f"\n🔄 转换: {dataset_type}/{sample_name}")
    
    # 检查必需文件
    required_files = ['adata1.h5ad', 'adata2.h5ad', 'embedding_before_clustering.npy']
    for file_name in required_files:
        if not (sample_dir / file_name).exists():
            print(f"❌ 缺少必需文件: {file_name}")
            return None
    
    try:
        # 1. 读取COSMOS结果
        adata1 = sc.read_h5ad(sample_dir / 'adata1.h5ad')  # RNA数据
        adata2 = sc.read_h5ad(sample_dir / 'adata2.h5ad')  # ADT/ATAC数据
        embedding = np.load(sample_dir / 'embedding_before_clustering.npy')
        
        # 读取预测标签(如果存在)
        cosmos_labels = None
        if (sample_dir / 'Y_Pred_label.npy').exists():
            cosmos_labels = np.load(sample_dir / 'Y_Pred_label.npy')
        
        print(f"   📊 adata1: {adata1.shape}, adata2: {adata2.shape}")
        print(f"   🧮 embedding: {embedding.shape}")
        
        # 2. 使用adata1作为基础 (通常是RNA数据，包含更完整的metadata)
        adata_integrated = adata1.copy()
        
        # 3. 添加集成embedding
        adata_integrated.obsm['COSMOS'] = embedding
        
        # 4. 移除可能存在的COSMOS原始聚类结果，避免混淆
        # 我们将使用自己的聚类方法，不保留COSMOS的预测标签
        
        # 5. 计算UMAP坐标
        print("   🗺️  计算UMAP坐标...")
        sc.pp.neighbors(adata_integrated, use_rep='COSMOS', n_neighbors=15)
        sc.tl.umap(adata_integrated)
        
        # 6. 获取聚类数量
        cluster_nums = get_cluster_numbers()
        
        # 确定数据集对应的聚类数量
        dataset_key = f"{dataset_type}_{sample_name}"
        if 'Fusion' in sample_name:  # 水平整合结果
            if 'Mouse_Thymus' in dataset_type:
                n_clusters = cluster_nums['Mouse_Thymus']
            elif 'Mouse_slpeen' in dataset_type:
                n_clusters = cluster_nums['Mouse_Spleen']
            elif 'Mouse_Brain' in dataset_type:
                n_clusters = cluster_nums['Mouse_Brain']
            else:
                # HLN_Fusion, HT_Fusion等，使用平均值
                n_clusters = 8  # 默认值
        else:
            n_clusters = cluster_nums.get(dataset_key, 8)  # 默认8
        
        print(f"   🎯 使用聚类数量: {n_clusters}")
        
        # 7. 读取训练时间
        running_time_path = sample_dir / 'running_time.npy'
        if running_time_path.exists():
            try:
                train_time = np.load(running_time_path)
                adata_integrated.uns['train_time'] = float(train_time)
                print(f"   ⏱️  训练时间: {train_time:.3f} 秒")
            except Exception as e:
                print(f"   ❌ 读取训练时间失败: {e}")
        else:
            print("   ⚠️  未找到running_time.npy文件")
        
        # 8. 进行4种聚类 (不使用COSMOS前缀)
        clustering_methods = ['mclust', 'leiden', 'louvain', 'kmeans']
        
        for method in clustering_methods:
            print(f"   🔬 运行{method}聚类...")
            try:
                adata_integrated = universal_clustering(
                    adata_integrated,
                    n_clusters=n_clusters,
                    used_obsm='COSMOS',
                    method=method,
                    key=method,  # 直接使用方法名，不加COSMOS前缀
                    use_pca=False  # 已经是embedding了
                )
            except Exception as e:
                print(f"      ❌ {method}聚类失败: {e}")
                # 创建占位符
                adata_integrated.obs[method] = '0'
        
        # 9. 方法信息已保存在train_time中，不需要额外元数据
        
        # 10. 尝试加载原始空间坐标
        original_adata = load_original_data(dataset_type, sample_name)
        if original_adata is not None and 'spatial' in original_adata.obsm:
            # 确保索引匹配
            if len(original_adata.obs) == len(adata_integrated.obs):
                adata_integrated.obsm['spatial'] = original_adata.obsm['spatial']
                print("   📍 添加空间坐标")
        
        # 11. 生成输出文件名
        # 生成与我们格式一致的命名
        if 'Fusion' in sample_name:
            # 水平整合结果
            if dataset_type == 'HLN':
                output_name = 't_COSMOS_HLN_Fusion.h5ad'
            elif dataset_type == 'HT':
                output_name = 't_COSMOS_HT_Fusion.h5ad'
            elif dataset_type == 'Mouse_Thymus':
                output_name = 't_COSMOS_MT_Fusion.h5ad'
            elif dataset_type == 'Mouse_slpeen':
                output_name = 't_COSMOS_MS_Fusion.h5ad'
            elif dataset_type == 'Mouse_Brain':
                output_name = 't_COSMOS_MB_Fusion.h5ad'
            else:
                output_name = f't_COSMOS_{dataset_type}_{sample_name}.h5ad'
        else:
            # 垂直整合结果
            if dataset_type == 'HLN':
                output_name = f't_COSMOS_HLN_{sample_name}.h5ad'
            elif dataset_type == 'HT':
                # slice1 -> S1, slice2 -> S2, slice3 -> S3
                slice_map = {'slice1': 'S1', 'slice2': 'S2', 'slice3': 'S3'}
                slice_id = slice_map.get(sample_name, sample_name)
                output_name = f't_COSMOS_HT_{slice_id}.h5ad'
            elif dataset_type.startswith('MISAR'):
                # MISAR_S1_E11 -> t_COSMOS_MISAR_S1_E11.h5ad
                output_name = f't_COSMOS_{dataset_type}_{sample_name}.h5ad'
            elif dataset_type == 'Mouse_Thymus':
                # Dataset3_Mouse_Thymus1 -> Thymus1
                thymus_id = sample_name.split('_')[-1]  # Mouse_Thymus1 -> Thymus1
                output_name = f't_COSMOS_MT_{thymus_id}.h5ad'
            elif dataset_type == 'Mouse_slpeen':
                # Dataset1_Mouse_Spleen1 -> Spleen1
                spleen_id = sample_name.split('_')[-1]  # Mouse_Spleen1 -> Spleen1
                output_name = f't_COSMOS_MS_{spleen_id}.h5ad'
            elif dataset_type == 'Mouse_Brain':
                # Dataset7_Mouse_Brain_ATAC -> MB_ATAC
                brain_type = sample_name.split('_')[-1]  # Mouse_Brain_ATAC -> ATAC
                output_name = f't_COSMOS_MB_{brain_type}.h5ad'
            else:
                output_name = f't_COSMOS_{dataset_type}_{sample_name}.h5ad'
        
        # 12. 保存结果
        output_path = output_dir / output_name
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        adata_integrated.write(output_path)
        print(f"   ✅ 保存到: {output_path}")
        
        return output_path
        
    except Exception as e:
        print(f"   ❌ 转换失败: {e}")
        return None

def main():
    """主函数"""
    print("🔄 开始转换COSMOS结果到SMOBench格式...")
    
    # Set R environment variables for mclust clustering
    os.environ['R_HOME'] = '/home/zhenghong/miniconda3/envs/smobench/lib/R'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    print("✅ R environment configured")
    
    # 输入和输出路径
    cosmos_results_dir = Path("/home/zhenghong/SMOBench-CLEAN/yjxiao/COSMOS/results")
    output_dir = Path("/home/zhenghong/SMOBench-CLEAN/Results/adata/t_COSMOS")
    
    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)
    
    converted_count = 0
    failed_count = 0
    
    # 遍历所有COSMOS结果
    for dataset_dir in cosmos_results_dir.iterdir():
        if dataset_dir.is_dir() and not dataset_dir.name.startswith('.'):
            print(f"\n📂 处理数据集类型: {dataset_dir.name}")
            
            for sample_dir in dataset_dir.iterdir():
                if sample_dir.is_dir():
                    result = convert_cosmos_sample(sample_dir, output_dir)
                    if result:
                        converted_count += 1
                    else:
                        failed_count += 1
    
    print(f"\n✅ 转换完成!")
    print(f"📊 统计: 成功 {converted_count} 个, 失败 {failed_count} 个")
    print(f"📁 结果保存在: {output_dir}")

if __name__ == "__main__":
    main()