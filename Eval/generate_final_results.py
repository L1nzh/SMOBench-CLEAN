import pandas as pd
import numpy as np
import os
import glob

# 定义指标到类别的映射关系（可配置，需要同时修改/src/demo中配置的指标）
metric_to_category = {
    'Moran Index': 'SC', 'Geary C': 'SC',
    
    'ARI': 'BioC', 'NMI': 'BioC', 'AMI': 'BioC', 'FMI': 'BioC', 
    'Purity': 'BioC', 'Homogeneity': 'BioC', 'Completeness': 'BioC', 
    'V-measure': 'BioC', 'Jaccard Index': 'BioC', 'Dice Index': 'BioC', 
    'F-measure': 'BioC', 'Silhouette Coefficient': 'BioC', 
    'Calinski-Harabaz Index': 'BioC', 'Davies-Bouldin Index': 'BioC', 
    'asw_celltype': 'BioC', 'graph_clisi': 'BioC',
    # 'Isolated Label F1': 'BioC', 'Isolated Silhouette': 'BioC',
    
    'Graph Connectivity': 'BER',
    'asw_batch': 'BER', 'Graph iLISI': 'BER', 'kBET': 'BER'
}

# 定义哪些指标是值越低越好（可配置）
# asw_batch,kBET本来也是越小越好，但是SCALE了
lower_is_better_metrics = ['Davies-Bouldin Index', 'Geary C']



def process_metrics_from_directory_vertical(directory_path):
    """
    对指定目录下的所有指标CSV文件执行完整的处理流程。
    步骤  1: 读取和分类
    步骤 2: 标准化
    步骤 3: 按类别求平均（批次级）
    步骤 4: 合并不同批次，求组织级平均
    """
    
    # --- 步骤 1: 读取文件并进行分类 ---
    print("\n--- 步骤 1: 读取与分类 ---")
    
    # 查找所有匹配的CSV文件
    file_paths_with_gt = glob.glob(os.path.join(directory_path, "*_cluster_metrics_with_GT.csv"))
    file_paths_wo_gt = glob.glob(os.path.join(directory_path, "*_cluster_metrics_wo_GT.csv"))
    file_paths = file_paths_with_gt + file_paths_wo_gt
    
    if not file_paths:
        print(f"错误: 在目录 '{directory_path}' 中未找到任何匹配的CSV文件。")
        return None

    # 读取所有文件到一个DataFrame列表中
    all_data = []
    for path in file_paths:
        filename = os.path.basename(path)
        parts = filename.split('_')
        method, tissue, batch = parts[0], parts[1], parts[2]
        
        df_temp = pd.read_csv(path, header=None, names=['Metric', 'Raw_Value'])
        df_temp['Method'] = method
        df_temp['Tissue'] = tissue
        df_temp['Batch'] = batch
        all_data.append(df_temp)
        
    # 将所有数据合并成一个大的DataFrame
    df_master = pd.concat(all_data, ignore_index=True)
    
    # 添加'Category'列
    df_master['Category'] = df_master['Metric'].map(metric_to_category)
    
    # 检查是否有未分类的指标
    unclassified = df_master[df_master['Category'].isna()]
    if not unclassified.empty:
        print("\n警告: 以下指标未找到分类，将被忽略:")
        print(unclassified['Metric'].unique())
        df_master.dropna(subset=['Category'], inplace=True) # 移除未分类的行
        
    print("数据读取和分类完成。数据预览:")
    print(df_master.head())

    # --- 步骤 2: 对每个数据做标准化 ---
    print("\n--- 步骤 2: 标准化 ---")
    
    def normalize_values(group):
        # 检查指标的方向性
        metric_name = group.name[2] # 获取分组的指标名称
        is_lower_better = metric_name in lower_is_better_metrics
        
        min_val = group.min()
        max_val = group.max()
        
        # 处理所有值都相同的情况，避免除以零
        if max_val == min_val:
            return pd.Series(0.5, index=group.index)
            
        if is_lower_better:
            # 值越低，分数越高
            normalized = (max_val - group) / (max_val - min_val)
        else:
            # 值越高，分数越高
            normalized = (group - min_val) / (max_val - min_val)
            
        return normalized

    # 按 组织-批次-指标 分组，然后对每个组应用标准化
    df_master['Normalized_Value'] = df_master.groupby(
        ['Tissue', 'Batch', 'Metric']
    )['Raw_Value'].transform(normalize_values)
    
    print("数据标准化完成。数据预览 (包含标准化得分):")
    # 显示一个指标在所有方法中的标准化结果作为示例
    print(df_master[df_master['Metric'] == 'ARI'][['Method', 'Raw_Value', 'Normalized_Value']].sort_values('Raw_Value'))
    print(df_master[df_master['Metric'] == 'Davies-Bouldin Index'][['Method', 'Raw_Value', 'Normalized_Value']].sort_values('Raw_Value'))

    # --- 步骤 3: 同类取平均（批次级） ---
    print("\n--- 步骤 3: 同类取平均（批次级） ---")

    # 按 方法-组织-批次-类别 分组，计算标准化得分的平均值
    batch_level_scores = df_master.groupby(
        ['Method', 'Tissue', 'Batch', 'Category']
    )['Normalized_Value'].mean()
    
    # 将结果从长格式转换为更易读的宽格式
    batch_summary = batch_level_scores.unstack(level='Category').reset_index()
    
    # 重命名列以符合期望的输出
    column_mapping = {'SC': 'SC_Score', 'BioC': 'BioC_Score', 'BER': 'BER_Score'}
    # 只重命名存在的列
    existing_columns = {k: v for k, v in column_mapping.items() if k in batch_summary.columns}
    batch_summary.rename(columns=existing_columns, inplace=True)
    
    print("批次级结果（每个批次单独的得分）:")
    print(batch_summary.to_string(index=False))

    # --- 新增步骤 4: 合并不同批次，计算组织级平均 ---
    print("\n--- 步骤 4: 合并批次，计算组织级综合得分 ---")
    
    # 确定实际存在的类别列
    category_columns = [col for col in ['SC_Score', 'BioC_Score', 'BER_Score'] if col in batch_summary.columns]
    
    # 按 方法-组织 分组，对不同批次的得分取平均
    tissue_level_scores = batch_summary.groupby(
        ['Method', 'Tissue']
    )[category_columns].mean().reset_index()
    
    # 保留两位小数（可选，增强可读性）
    tissue_level_scores[category_columns] = tissue_level_scores[category_columns].round(2)
    
    print("组织级综合结果（合并不同批次后的平均得分）:")
    return tissue_level_scores


def process_metrics_from_directory_horizontal(directory_path):
    """
    处理horizontal任务的输出文件（以_BC.csv为后缀）
    步骤 1: 读取和分类
    步骤 2: 标准化（按方法-组织-指标）
    步骤 3: 按类别求平均
    """
    
    # --- 步骤 1: 读取文件并进行分类 ---
    print("\n--- 步骤 1: 读取与分类 (Horizontal) ---")
    
    # 查找所有匹配的BC CSV文件
    file_paths = glob.glob(os.path.join(directory_path, "*_BC.csv"))
    
    if not file_paths:
        print(f"错误: 在目录 '{directory_path}' 中未找到任何BC匹配的CSV文件。")
        return None

    # 读取所有文件到一个DataFrame列表中
    all_data = []
    for path in file_paths:
        filename = os.path.basename(path)
        # 文件名格式如: SpatialGlue_HT_BC.csv
        parts = filename.split('_')
        method, tissue = parts[0], parts[1]
        
        df_temp = pd.read_csv(path, header=None, names=['Metric', 'Raw_Value'])
        df_temp['Method'] = method
        df_temp['Tissue'] = tissue
        all_data.append(df_temp)
        
    # 将所有数据合并成一个大的DataFrame
    df_master = pd.concat(all_data, ignore_index=True)
    
    # 添加'Category'列
    df_master['Category'] = df_master['Metric'].map(metric_to_category)
    
    # 检查是否有未分类的指标
    unclassified = df_master[df_master['Category'].isna()]
    if not unclassified.empty:
        print("\n警告: 以下指标未找到分类，将被忽略:")
        print(unclassified['Metric'].unique())
        df_master.dropna(subset=['Category'], inplace=True) # 移除未分类的行
        
    print("数据读取和分类完成。数据预览:")
    print(df_master.head())

    # --- 步骤 2: 对每个数据做标准化 ---
    print("\n--- 步骤 2: 标准化 (Horizontal) ---")
    
    def normalize_values(group):
        # 检查指标的方向性
        metric_name = group.name[1] # 获取分组的指标名称
        is_lower_better = metric_name in lower_is_better_metrics
        
        min_val = group.min()
        max_val = group.max()
        
        # 处理所有值都相同的情况，避免除以零
        if max_val == min_val:
            return pd.Series(0.5, index=group.index)
            
        if is_lower_better:
            # 值越低，分数越高
            normalized = (max_val - group) / (max_val - min_val)
        else:
            # 值越高，分数越高
            normalized = (group - min_val) / (max_val - min_val)
            
        return normalized

    # 按 组织-指标 分组，然后对每个组应用标准化
    df_master['Normalized_Value'] = df_master.groupby(
        ['Tissue', 'Metric']
    )['Raw_Value'].transform(normalize_values)
    
    print("数据标准化完成。数据预览ARI DBI (包含标准化得分):")
    # 显示一个指标在所有方法中的标准化结果作为示例
    print(df_master[df_master['Metric'] == 'Graph Connectivity'][['Method', 'Raw_Value', 'Normalized_Value']].sort_values('Raw_Value'))

    # --- 步骤 3: 同类取平均 ---
    print("\n--- 步骤 3: 同类取平均 (Horizontal) ---")

    # 按 方法-组织-类别 分组，计算标准化得分的平均值
    tissue_level_scores = df_master.groupby(
        ['Method', 'Tissue', 'Category']
    )['Normalized_Value'].mean()
    
    # 将结果从长格式转换为更易读的宽格式
    final_summary = tissue_level_scores.unstack(level='Category').reset_index()
    
    # 重命名列以符合期望的输出
    column_mapping = {'SC': 'SC_Score', 'BioC': 'BioC_Score', 'BER': 'BER_Score'}
    # 只重命名存在的列
    existing_columns = {k: v for k, v in column_mapping.items() if k in final_summary.columns}
    final_summary.rename(columns=existing_columns, inplace=True)
    
    # 确定实际存在的类别列
    category_columns = [col for col in ['SC_Score', 'BioC_Score', 'BER_Score'] if col in final_summary.columns]
    
    # 保留两位小数（可选，增强可读性）
    final_summary[category_columns] = final_summary[category_columns].round(2)
    
    print("Horizontal任务综合结果:")
    return final_summary

def generate_final_SMOBench(results_dir):
    """
    读取final_vertical.csv和final_horizontal.csv，将同一方法同一组织的得分取平均，
    生成最终的SMOBench综合评估结果
    """
    print("\n" + "="*60)
    print("生成 SMOBench 综合评估结果")
    print("="*60)
    
    # 读取vertical和horizontal结果
    vertical_path = os.path.join(results_dir, 'final_vertical.csv')
    horizontal_path = os.path.join(results_dir, 'final_horizontal.csv')
    
    dataframes = []
    
    # 读取vertical结果
    if os.path.exists(vertical_path):
        df_vertical = pd.read_csv(vertical_path)
        dataframes.append(df_vertical)
        print(f"已加载 Vertical 结果: {vertical_path}")
    else:
        print(f"警告: 未找到 Vertical 结果文件 {vertical_path}")
    
    # 读取horizontal结果
    if os.path.exists(horizontal_path):
        df_horizontal = pd.read_csv(horizontal_path)
        dataframes.append(df_horizontal)
        print(f"已加载 Horizontal 结果: {horizontal_path}")
    else:
        print(f"警告: 未找到 Horizontal 结果文件 {horizontal_path}")
    
    # 如果没有数据文件，返回None
    if not dataframes:
        print("错误: 没有找到任何输入文件")
        return None
    
    # 合并所有数据
    df_combined = pd.concat(dataframes, ignore_index=True)
    
    # 确定实际存在的类别列
    category_columns = [col for col in ['SC_Score', 'BioC_Score', 'BER_Score'] if col in df_combined.columns]
    
    # 按 方法-组织 分组，计算每个类别得分的平均值
    grouped = df_combined.groupby(['Method', 'Tissue'])[category_columns].mean().reset_index()
    
    # 计算SMOBench综合得分（各维度得分的平均值）
    if len(category_columns) > 0:
        grouped['SMOBench'] = grouped[category_columns].mean(axis=1).round(2)
    
    print("\nSMOBench综合评估结果:")
    return grouped


# --- 执行代码 ---
if __name__ == "__main__":
    data_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results/test')
    
    # 创建results目录用于保存结果
    results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
    os.makedirs(results_dir, exist_ok=True)
    
    # 运行vertical处理流程
    print("="*60)
    print("处理 Vertical 任务结果")
    print("="*60)
    vertical_summary = process_metrics_from_directory_vertical(data_directory)
    
    # 保存vertical结果
    if vertical_summary is not None:
        vertical_output_path = os.path.join(results_dir, 'final_vertical.csv')
        vertical_summary.to_csv(vertical_output_path, index=False)
        print(f"\nVertical结果已保存至: {vertical_output_path}")
        
        print("\n" + "="*50)
        print("          Vertical任务组织级综合性能评估得分")
        print("="*50)
        print(vertical_summary.to_string(index=False))
        print("="*50)
    
    # 运行horizontal处理流程
    print("\n" + "="*60)
    print("处理 Horizontal 任务结果")
    print("="*60)
    horizontal_summary = process_metrics_from_directory_horizontal(data_directory)
    
    # 保存horizontal结果
    if horizontal_summary is not None:
        horizontal_output_path = os.path.join(results_dir, 'final_horizontal.csv')
        horizontal_summary.to_csv(horizontal_output_path, index=False)
        print(f"\nHorizontal结果已保存至: {horizontal_output_path}")
        
        print("\n" + "="*50)
        print("          Horizontal任务组织级综合性能评估得分")
        print("="*50)
        print(horizontal_summary.to_string(index=False))
        print("="*50)
    
    # 生成最终的SMOBench综合评估结果
    final_SMOBench_result = generate_final_SMOBench(results_dir)
    
    # 保存并显示最终结果
    if final_SMOBench_result is not None:
        final_SMOBench_path = os.path.join(results_dir, 'final_SMOBench.csv')
        final_SMOBench_result.to_csv(final_SMOBench_path, index=False)
        print(f"\nSMOBench综合结果已保存至: {final_SMOBench_path}")
        
        print("\n" + "="*50)
        print("          SMOBench综合性能评估得分")
        print("="*50)
        print(final_SMOBench_result.to_string(index=False))
        print("="*50)