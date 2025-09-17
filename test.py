#!/usr/bin/env python3
"""
SMOBench 聚类数验证脚本
检查SpaMV和CANDIES所有结果的聚类数是否符合标准
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import glob

def get_standard_cluster_numbers():
    """获取SMOBench标准聚类数"""
    return {
        # withGT datasets
        'HLN_A1': 10, 'HLN_D1': 11,
        'HT_S1': 4, 'HT_S2': 5, 'HT_S3': 5,
        'MISAR_S1_E11': 8, 'MISAR_S1_E13': 12, 'MISAR_S1_E15': 12, 'MISAR_S1_E18': 14,
        'MISAR_S2_E11': 13, 'MISAR_S2_E13': 14, 'MISAR_S2_E15': 15, 'MISAR_S2_E18': 16,
        
        # woGT datasets  
        'Mouse_Thymus1': 8, 'Mouse_Thymus2': 8, 'Mouse_Thymus3': 8, 'Mouse_Thymus4': 8,
        'Mouse_Spleen1': 5, 'Mouse_Spleen2': 5,
        'Mouse_Brain_ATAC': 18, 'Mouse_Brain_H3K4me3': 18, 'Mouse_Brain_H3K27ac': 18, 'Mouse_Brain_H3K27me3': 18
    }

def extract_dataset_key(file_path):
    """从文件路径提取数据集键名"""
    path_parts = file_path.split('/')
    
    if 'SpaMV' in path_parts:
        if 'HLN' in path_parts:
            return f"HLN_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'HT' in path_parts:
            return f"HT_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'MISAR_S1' in path_parts:
            return f"MISAR_S1_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'MISAR_S2' in path_parts:
            return f"MISAR_S2_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'Mouse_Thymus' in path_parts:
            return f"Mouse_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'Mouse_Spleen' in path_parts:
            return f"Mouse_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'Mouse_Brain' in path_parts:
            return f"Mouse_Brain_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
    
    elif 'CANDIES' in path_parts:
        if 'HLN' in path_parts:
            return f"HLN_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'HT' in path_parts:
            return f"HT_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'MISAR_S1' in path_parts:
            return f"MISAR_S1_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'MISAR_S2' in path_parts:
            return f"MISAR_S2_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'Mouse_Thymus' in path_parts:
            return f"Mouse_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'Mouse_Spleen' in path_parts:
            return f"Mouse_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
        elif 'Mouse_Brain' in path_parts:
            return f"Mouse_Brain_{path_parts[-1].split('_')[-1].replace('.h5ad', '')}"
    
    return "Unknown"

def check_clustering_numbers(method_name, file_pattern):
    """检查指定方法的所有聚类数"""
    print(f"\n{'='*60}")
    print(f"🔍 检查 {method_name} 聚类数")
    print(f"{'='*60}")
    
    standard_clusters = get_standard_cluster_numbers()
    clustering_methods = ['mclust', 'leiden', 'louvain', 'kmeans']
    
    files = glob.glob(file_pattern, recursive=True)
    files.sort()
    
    total_files = len(files)
    correct_files = 0
    issues_found = []
    
    print(f"找到 {total_files} 个 {method_name} 结果文件\n")
    
    for file_path in files:
        dataset_key = extract_dataset_key(file_path)
        expected_clusters = standard_clusters.get(dataset_key, "Unknown")
        
        try:
            adata = sc.read_h5ad(file_path)
            
            # 获取简短的文件名用于显示
            short_name = '/'.join(file_path.split('/')[-4:])
            
            print(f"📁 {short_name}")
            print(f"   数据集: {dataset_key} | 期望聚类数: {expected_clusters}")
            
            file_correct = True
            method_results = []
            
            for cluster_method in clustering_methods:
                if cluster_method in adata.obs.columns:
                    actual_clusters = adata.obs[cluster_method].nunique()
                    is_correct = actual_clusters == expected_clusters
                    status = '✅' if is_correct else '❌'
                    
                    method_results.append(f"{cluster_method}: {actual_clusters} {status}")
                    
                    if not is_correct:
                        file_correct = False
                        issues_found.append({
                            'file': short_name,
                            'dataset': dataset_key,
                            'method': cluster_method,
                            'expected': expected_clusters,
                            'actual': actual_clusters
                        })
                else:
                    method_results.append(f"{cluster_method}: 缺失 ❌")
                    file_correct = False
                    issues_found.append({
                        'file': short_name,
                        'dataset': dataset_key,
                        'method': cluster_method,
                        'expected': expected_clusters,
                        'actual': 'Missing'
                    })
            
            print(f"   聚类结果: {' | '.join(method_results)}")
            
            if file_correct:
                correct_files += 1
                print(f"   状态: ✅ 全部正确")
            else:
                print(f"   状态: ❌ 发现问题")
            
            print()
            
        except Exception as e:
            print(f"   ❌ 读取失败: {e}")
            issues_found.append({
                'file': short_name,
                'dataset': dataset_key,
                'method': 'ALL',
                'expected': expected_clusters,
                'actual': f'Error: {e}'
            })
            print()
    
    # 总结报告
    print(f"📊 {method_name} 总结:")
    print(f"   总文件数: {total_files}")
    print(f"   完全正确: {correct_files}")
    print(f"   有问题的: {total_files - correct_files}")
    print(f"   准确率: {correct_files/total_files*100:.1f}%")
    
    if issues_found:
        print(f"\n⚠️  发现的问题:")
        for issue in issues_found:
            print(f"   - {issue['file']} | {issue['dataset']} | {issue['method']}: "
                  f"期望{issue['expected']}, 实际{issue['actual']}")
    else:
        print(f"\n🎉 {method_name} 所有聚类数完全正确!")
    
    return correct_files, total_files, issues_found

def main():
    """主函数"""
    print("🧪 SMOBench 聚类数全面验证")
    print("=" * 60)
    
    all_issues = []
    
    # 检查SpaMV
    spamv_correct, spamv_total, spamv_issues = check_clustering_numbers(
        "SpaMV", 
        "/home/zhenghong/SMOBench-CLEAN/Results/adata/SpaMV/**/*.h5ad"
    )
    all_issues.extend(spamv_issues)
    
    # 检查CANDIES
    candies_correct, candies_total, candies_issues = check_clustering_numbers(
        "CANDIES",
        "/home/zhenghong/SMOBench-CLEAN/Results/adata/CANDIES/**/*.h5ad"
    )
    all_issues.extend(candies_issues)
    
    # 全局总结
    print(f"\n{'='*60}")
    print("🎯 全局总结")
    print(f"{'='*60}")
    
    total_correct = spamv_correct + candies_correct
    total_files = spamv_total + candies_total
    
    print(f"SpaMV:   {spamv_correct}/{spamv_total} 正确 ({spamv_correct/spamv_total*100:.1f}%)")
    print(f"CANDIES: {candies_correct}/{candies_total} 正确 ({candies_correct/candies_total*100:.1f}%)")
    print(f"总计:    {total_correct}/{total_files} 正确 ({total_correct/total_files*100:.1f}%)")
    
    if all_issues:
        print(f"\n❌ 总共发现 {len(all_issues)} 个问题")
        print("建议检查这些文件的聚类参数设置")
    else:
        print(f"\n🎉 所有文件的聚类数都完全正确!")
        print("SpaMV和CANDIES完全符合SMOBench标准!")

if __name__ == "__main__":
    main()