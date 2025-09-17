#!/usr/bin/env python3
"""
Convert yrliu (PRESENT) results to SMOBench standard format
将PRESENT方法的结果转换为SMOBench标准格式
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
import sys
warnings.filterwarnings('ignore')

# Add project root to path
sys.path.append('/home/zhenghong/SMOBench-CLEAN')
from Utils.SMOBench_clustering import universal_clustering

def get_dataset_mapping():
    """获取yrliu目录名到SMOBench数据集的映射"""
    return {
        # withGT datasets - 直接映射
        'A1_output': ('HLN', 'A1'),
        'D1_output': ('HLN', 'D1'), 
        'S1_output': ('HT', 'S1'),
        'S2_output': ('HT', 'S2'),
        'S3_output': ('HT', 'S3'),
        'S1_E11_output': ('MISAR_S1', 'E11'),
        'S1_E13_output': ('MISAR_S1', 'E13'),
        'S1_E15_output': ('MISAR_S1', 'E15'),
        'S1_E18_output': ('MISAR_S1', 'E18'),
        'S2_E11_output': ('MISAR_S2', 'E11'),
        'S2_E13_output': ('MISAR_S2', 'E13'),
        'S2_E15_output': ('MISAR_S2', 'E15'),
        'S2_E18_output': ('MISAR_S2', 'E18'),
        
        # woGT datasets - 基于Dataset编号映射
        'Dataset1_output': ('Mouse_Spleen', 'Spleen1'),
        'Dataset2_output': ('Mouse_Spleen', 'Spleen2'),
        'Dataset3_output': ('Mouse_Thymus', 'Thymus1'),
        'Dataset4_output': ('Mouse_Thymus', 'Thymus2'),
        'Dataset5_output': ('Mouse_Thymus', 'Thymus3'),
        'Dataset6_output': ('Mouse_Thymus', 'Thymus4'),
        'Dataset7_output': ('Mouse_Brain', 'ATAC'),
        'Dataset8_output': ('Mouse_Brain', 'H3K4me3'),
        'Dataset9_output': ('Mouse_Brain', 'H3K27ac'),
        'Dataset10_output': ('Mouse_Brain', 'H3K27me3'),
        
        # Fusion datasets (水平整合)
        'HLN_output': ('HLN', 'Fusion'),
        'HT_output': ('HT', 'Fusion'),
        'MISAR_S1_output': ('MISAR_S1', 'Fusion'),
        'MISAR_S2_output': ('MISAR_S2', 'Fusion'),
        'Mouse_Thymus_output': ('Mouse_Thymus', 'Fusion'),
        'Mouse_slpeen_output': ('Mouse_Spleen', 'Fusion'),
        'Mouse_Brain_multi_output': ('Mouse_Brain', 'Fusion'),
    }

def get_cluster_numbers():
    """获取SMOBench标准聚类数"""
    return {
        # withGT datasets
        'HLN_A1': 10, 'HLN_D1': 11,
        'HT_S1': 4, 'HT_S2': 5, 'HT_S3': 5,
        'MISAR_S1_E11': 8, 'MISAR_S1_E13': 12, 'MISAR_S1_E15': 12, 'MISAR_S1_E18': 14,
        'MISAR_S2_E11': 13, 'MISAR_S2_E13': 14, 'MISAR_S2_E15': 15, 'MISAR_S2_E18': 16,
        
        # woGT datasets  
        'Mouse_Thymus': 8, 'Mouse_Spleen': 5, 'Mouse_Brain': 18,
        
        # Fusion datasets (水平整合，使用平均值或合理值)
        'HLN_Fusion': 10, 'HT_Fusion': 5, 
        'MISAR_S1_Fusion': 12, 'MISAR_S2_Fusion': 14,
        'Mouse_Thymus_Fusion': 8, 'Mouse_Spleen_Fusion': 5, 'Mouse_Brain_Fusion': 18
    }

def load_original_adata(dataset_category, subset_name):
    """加载原始AnnData获取metadata和空间坐标"""
    base_path = "/home/zhenghong/SMOBench-CLEAN/Dataset"
    
    # 确定数据类型和路径
    if dataset_category in ['HLN', 'HT']:
        if subset_name == 'Fusion':
            return None  # Fusion数据不加载原始数据
        modality = 'RNA_ADT'
        data_type = 'withGT'
        if dataset_category == 'HLN':
            original_dataset = 'Human_Lymph_Nodes'
        else:
            original_dataset = 'Human_Tonsils'
    elif dataset_category in ['MISAR_S1', 'MISAR_S2']:
        if subset_name == 'Fusion':
            return None
        modality = 'RNA_ATAC'
        data_type = 'withGT'
        if dataset_category == 'MISAR_S1':
            original_dataset = 'Mouse_Embryos_S1'
        else:
            original_dataset = 'Mouse_Embryos_S2'
    elif dataset_category in ['Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain']:
        if subset_name == 'Fusion':
            return None
        data_type = 'woGT'
        if dataset_category == 'Mouse_Brain':
            modality = 'RNA_ATAC'
            original_dataset = f'Mouse_Brain'
        else:
            modality = 'RNA_ADT'
            original_dataset = dataset_category
    else:
        return None
    
    # 构建路径
    rna_path = f"{base_path}/{data_type}/{modality}/{original_dataset}/{subset_name}/adata_RNA.h5ad"
    
    if os.path.exists(rna_path):
        try:
            return sc.read_h5ad(rna_path)
        except Exception as e:
            print(f"Warning: Error loading original data {rna_path}: {e}")
    
    return None

def convert_yrliu_sample(sample_dir, output_dir, method_name="PRESENT"):
    """转换单个yrliu样本到SMOBench格式"""
    
    sample_path = Path(sample_dir)
    sample_name = sample_path.name
    
    print(f"\n🔄 转换: {sample_name}")
    
    # 获取数据集映射
    dataset_mapping = get_dataset_mapping()
    
    if sample_name not in dataset_mapping:
        print(f"⏭️  跳过未知样本: {sample_name}")
        return None
    
    dataset_category, subset_name = dataset_mapping[sample_name]
    print(f"   映射到: {dataset_category}/{subset_name}")
    
    # 检查必需文件
    adata_file = sample_path / 'adata_output.h5ad'
    embeddings_file = sample_path / 'embeddings_output.csv'
    domains_file = sample_path / 'domains_output.csv'
    
    required_files = [adata_file, embeddings_file, domains_file]
    missing_files = [f for f in required_files if not f.exists()]
    
    if missing_files:
        print(f"❌ 缺少必需文件: {[f.name for f in missing_files]}")
        return None
    
    try:
        # 1. 读取PRESENT结果
        adata = sc.read_h5ad(adata_file)
        embeddings_df = pd.read_csv(embeddings_file, index_col=0)
        domains_df = pd.read_csv(domains_file, index_col=0)
        
        print(f"   📊 adata: {adata.shape}")
        print(f"   📊 embeddings: {embeddings_df.shape}")
        print(f"   📊 domains: {domains_df.shape}")
        
        # 2. 添加PRESENT嵌入
        # 假设embeddings_df就是PRESENT的主要嵌入
        embeddings_array = embeddings_df.values
        if embeddings_array.shape[0] == adata.n_obs:
            adata.obsm[method_name] = embeddings_array
            print(f"   ✅ 添加{method_name}嵌入: {embeddings_array.shape}")
        else:
            print(f"   ❌ 嵌入维度不匹配: adata {adata.n_obs}, embeddings {embeddings_array.shape[0]}")
            return None
        
        # 3. 添加PRESENT预测的domains作为初始聚类
        if 'domains' in domains_df.columns and len(domains_df) == adata.n_obs:
            adata.obs['PRESENT_domains'] = pd.Categorical(domains_df['domains'].astype(str))
            print(f"   ✅ 添加PRESENT domains")
        
        # 4. 计算UMAP坐标
        print("   🗺️  计算UMAP坐标...")
        sc.pp.neighbors(adata, use_rep=method_name, n_neighbors=30)
        sc.tl.umap(adata)
        
        # 5. 获取聚类数量
        cluster_nums = get_cluster_numbers()
        dataset_key = f"{dataset_category}_{subset_name}"
        if subset_name == 'Fusion':
            dataset_key = f"{dataset_category}_Fusion"
        elif dataset_category in ['Mouse_Thymus', 'Mouse_Spleen', 'Mouse_Brain']:
            dataset_key = dataset_category
        
        n_clusters = cluster_nums.get(dataset_key, 8)
        print(f"   🎯 使用聚类数量: {n_clusters}")
        
        # 6. 进行4种聚类
        clustering_methods = ['mclust', 'leiden', 'louvain', 'kmeans']
        
        for method in clustering_methods:
            print(f"   🔬 运行{method}聚类...")
            try:
                adata = universal_clustering(
                    adata,
                    n_clusters=n_clusters,
                    used_obsm=method_name,
                    method=method,
                    key=method,
                    use_pca=False
                )
                print(f"      ✅ {method}聚类完成")
            except Exception as e:
                print(f"      ❌ {method}聚类失败: {e}")
                # 创建占位符
                adata.obs[method] = pd.Categorical(['0'] * adata.n_obs)
        
        # 7. 尝试加载原始空间坐标
        original_adata = load_original_adata(dataset_category, subset_name)
        if original_adata is not None and 'spatial' in original_adata.obsm:
            if len(original_adata.obs) == len(adata.obs):
                adata.obsm['spatial'] = original_adata.obsm['spatial']
                print("   📍 添加空间坐标")
        
        # 8. 添加训练时间 (如果有的话，这里用默认值)
        adata.uns['train_time'] = 0.0  # PRESENT结果没有训练时间信息
        
        # 9. 生成输出文件名和路径
        if subset_name == 'Fusion':
            if dataset_category == 'HLN':
                output_filename = f"t_{method_name}_HLN_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "HLN" / "Fusion"
            elif dataset_category == 'HT':
                output_filename = f"t_{method_name}_HT_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "HT" / "Fusion"
            elif dataset_category == 'MISAR_S1':
                output_filename = f"t_{method_name}_MISAR_S1_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "MISAR_S1" / "Fusion"
            elif dataset_category == 'MISAR_S2':
                output_filename = f"t_{method_name}_MISAR_S2_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "MISAR_S2" / "Fusion"
            elif dataset_category == 'Mouse_Thymus':
                output_filename = f"t_{method_name}_MT_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "Mouse_Thymus" / "Fusion"
            elif dataset_category == 'Mouse_Spleen':
                output_filename = f"t_{method_name}_MS_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "Mouse_Spleen" / "Fusion"
            elif dataset_category == 'Mouse_Brain':
                output_filename = f"t_{method_name}_MB_Fusion.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "Mouse_Brain" / "Fusion"
        else:
            # 垂直整合结果
            if dataset_category == 'HLN':
                output_filename = f"t_{method_name}_HLN_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "HLN" / subset_name
            elif dataset_category == 'HT':
                output_filename = f"t_{method_name}_HT_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "HT" / subset_name
            elif dataset_category == 'MISAR_S1':
                output_filename = f"t_{method_name}_MISAR_S1_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "MISAR_S1" / subset_name
            elif dataset_category == 'MISAR_S2':
                output_filename = f"t_{method_name}_MISAR_S2_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "MISAR_S2" / subset_name
            elif dataset_category == 'Mouse_Thymus':
                output_filename = f"t_{method_name}_MT_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "Mouse_Thymus" / subset_name
            elif dataset_category == 'Mouse_Spleen':
                output_filename = f"t_{method_name}_MS_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "Mouse_Spleen" / subset_name
            elif dataset_category == 'Mouse_Brain':
                output_filename = f"t_{method_name}_MB_{subset_name}.h5ad"
                output_dataset_dir = Path(output_dir) / f"t_{method_name}" / "Mouse_Brain" / subset_name
        
        # 10. 保存结果
        output_dataset_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dataset_dir / output_filename
        
        adata.write(output_path)
        print(f"   ✅ 保存到: {output_path}")
        
        return output_path
        
    except Exception as e:
        print(f"   ❌ 转换失败: {e}")
        return None

def main():
    """主函数"""
    print("🔄 开始转换yrliu (PRESENT) 结果到SMOBench格式...")
    
    # Set R environment variables for mclust clustering
    os.environ['R_HOME'] = '/home/zhenghong/miniconda3/envs/smobench/lib/R'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    print("✅ R environment configured")
    
    # 输入和输出路径
    yrliu_dir = Path("/home/zhenghong/SMOBench-CLEAN/yrliu")
    output_dir = Path("/home/zhenghong/SMOBench-CLEAN/Results/adata")
    
    converted_count = 0
    failed_count = 0
    
    # 获取所有输出目录
    output_dirs = [d for d in yrliu_dir.iterdir() 
                   if d.is_dir() and d.name.endswith('_output') 
                   and not d.name.startswith('.')]
    
    output_dirs.sort()
    
    print(f"找到 {len(output_dirs)} 个yrliu输出目录")
    
    # 转换每个样本
    for sample_dir in output_dirs:
        result = convert_yrliu_sample(sample_dir, output_dir)
        if result:
            converted_count += 1
        else:
            failed_count += 1
    
    print(f"\n✅ 转换完成!")
    print(f"📊 统计: 成功 {converted_count} 个, 失败 {failed_count} 个")
    print(f"📁 结果保存在: {output_dir}/t_PRESENT/")

if __name__ == "__main__":
    main()