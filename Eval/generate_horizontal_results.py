#!/usr/bin/env python3
"""
Generate final results for horizontal integration evaluation
"""

import os
import glob
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

# Import functions from vertical integration
from generate_final_results import (
    normalize_metric_value, 
    aggregate_by_dataset,
    save_results
)

# Horizontal integration specific constants
HORIZONTAL_SC_METRICS = ['Moran Index']
HORIZONTAL_BER_METRICS = ['kBET', 'KNN_connectivity', 'bASW', 'iLISI', 'PCR']
HORIZONTAL_BIOC_METRICS_WITHGT = ['ARI', 'NMI', 'asw_celltype', 'graph_clisi']
HORIZONTAL_BIOC_METRICS_WOGT = ['Davies-Bouldin Index', 'Silhouette Coefficient', 'Calinski-Harabaz Index']

LOWER_IS_BETTER = ['Davies-Bouldin Index', 'Geary C']

def calculate_horizontal_normalized_scores(df):
    """
    Calculate normalized scores for horizontal integration with three dimensions
    
    Parameters:
    -----------
    df : pd.DataFrame
        Aggregated results DataFrame
    
    Returns:
    --------
    pd.DataFrame : DataFrame with normalized scores and dimension scores
    """
    df = df.copy()
    
    # Group by clustering and GT availability for normalization
    for (clustering, has_gt), subset_df in df.groupby(['Clustering', 'GT_Available']):
        mask = (df['Clustering'] == clustering) & (df['GT_Available'] == has_gt)
        
        if len(subset_df) < 2:
            continue  # Skip normalization if only one data point
        
        # Process SC metrics (no normalization needed - already in good range)
        for metric in HORIZONTAL_SC_METRICS:
            if metric in subset_df.columns:
                df.loc[mask, 'SC_Score'] = df.loc[mask, metric]
        
        # Process BioC metrics
        if has_gt:
            bioc_metrics = HORIZONTAL_BIOC_METRICS_WITHGT
        else:
            bioc_metrics = HORIZONTAL_BIOC_METRICS_WOGT
        
        bioc_scores = []
        bioc_normalized_scores = []
        for metric in bioc_metrics:
            if metric in subset_df.columns:
                # Only normalize DBI and CHI
                if metric in ['Davies-Bouldin Index', 'Calinski-Harabaz Index']:
                    all_values = subset_df[metric].dropna().tolist()
                    if all_values:
                        normalized_values = [normalize_metric_value(val, metric, all_values) 
                                           for val in subset_df[metric]]
                        df.loc[mask, f'{metric}_normalized'] = normalized_values
                        bioc_normalized_scores.append(f'{metric}_normalized')
                else:
                    # Other BioC metrics use original values
                    bioc_scores.append(metric)
        
        # Calculate BioC score using mix of normalized (DBI/CHI) and original values
        all_bioc_scores = bioc_scores + bioc_normalized_scores
        if all_bioc_scores:
            df.loc[mask, 'BVC_Score'] = df.loc[mask, all_bioc_scores].mean(axis=1)
        
        # Process BER metrics (no normalization needed for most)
        ber_values = []
        for metric in HORIZONTAL_BER_METRICS:
            if metric in subset_df.columns:
                ber_values.append(metric)
        
        if ber_values:
            df.loc[mask, 'BER_Score'] = df.loc[mask, ber_values].mean(axis=1)
        
        # Calculate Final Score (three dimensions)
        score_columns = []
        if 'SC_Score' in df.columns:
            score_columns.append('SC_Score')
        if 'BVC_Score' in df.columns:
            score_columns.append('BVC_Score')
        if 'BER_Score' in df.columns:
            score_columns.append('BER_Score')
        
        if score_columns:
            df.loc[mask, 'Final_Score'] = df.loc[mask, score_columns].mean(axis=1)
    
    return df

def create_horizontal_summary_tables(df, clustering_method):
    """
    Create summary tables for horizontal integration
    
    Parameters:
    -----------
    df : pd.DataFrame
        Processed results DataFrame
    clustering_method : str
        Clustering method name
    
    Returns:
    --------
    tuple : (rna_adt_table, rna_atac_table, comprehensive_table)
    """
    
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    # Split by data type
    rna_adt_df = df[df['Dataset_Type'] == 'RNA_ADT'].copy() if 'Dataset_Type' in df.columns else pd.DataFrame()
    rna_atac_df = df[df['Dataset_Type'] == 'RNA_ATAC'].copy() if 'Dataset_Type' in df.columns else pd.DataFrame()
    
    def create_summary_table(subset_df, table_name):
        if subset_df.empty:
            return pd.DataFrame()
        
        # Pivot table with Method as rows, Dataset as columns, Final_Score as values
        pivot_df = subset_df.pivot_table(
            index='Method',
            columns='Dataset', 
            values='Final_Score',
            aggfunc='first'  # Take first value if duplicates
        )
        
        # Calculate overall score (mean across datasets)
        if not pivot_df.empty:
            pivot_df['Overall_Score'] = pivot_df.mean(axis=1)
            # Sort by overall score (descending)
            pivot_df = pivot_df.sort_values('Overall_Score', ascending=False)
            # Round to 3 decimal places
            pivot_df = pivot_df.round(3)
        
        return pivot_df
    
    # Create summary tables
    rna_adt_table = create_summary_table(rna_adt_df, "RNA_ADT")
    rna_atac_table = create_summary_table(rna_atac_df, "RNA_ATAC") 
    comprehensive_table = create_summary_table(df, "Comprehensive")
    
    return rna_adt_table, rna_atac_table, comprehensive_table

def process_horizontal_results(results_dir):
    """
    Process horizontal integration results
    
    Parameters:
    -----------
    results_dir : str
        Directory containing horizontal integration results
    
    Returns:
    --------
    dict : Results organized by clustering method
    """
    
    print(f"Processing horizontal integration results from: {results_dir}")
    
    # Find all result files
    withgt_files = glob.glob(os.path.join(results_dir, "**/*_withGT.csv"), recursive=True)
    wogt_files = glob.glob(os.path.join(results_dir, "**/*_woGT.csv"), recursive=True)
    
    all_files = withgt_files + wogt_files
    
    if not all_files:
        print(f"No evaluation result files found in {results_dir}")
        return {}
    
    print(f"Found {len(all_files)} evaluation result files")
    
    # Group results by clustering method
    results_by_clustering = {
        'leiden': [],
        'louvain': [],
        'kmeans': [],
        'mclust': []
    }
    
    # Process each file
    for file_path in all_files:
        try:
            # Parse filename to extract metadata
            filename = os.path.basename(file_path)
            parts = filename.replace('.csv', '').split('_')
            
            if len(parts) < 4:
                print(f"Warning: Cannot parse filename {filename}")
                continue
            
            # Extract metadata from filename
            method_name = parts[0]
            dataset_name = parts[1]
            
            # For horizontal integration, format is: METHOD_DATASET_horizontal_CLUSTERING_GT.csv
            if 'horizontal' in parts:
                horizontal_idx = parts.index('horizontal')
                if len(parts) > horizontal_idx + 1:
                    clustering_method = parts[horizontal_idx + 1]
                    gt_status = parts[-1]  # withGT or woGT
                else:
                    continue
            else:
                continue
            
            if clustering_method not in results_by_clustering:
                continue
            
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            if df.empty:
                continue
            
            # Horizontal integration CSV files have Metric,Value format - need to transpose
            if 'Metric' in df.columns and 'Value' in df.columns:
                # Convert key-value pairs to single row dictionary
                metrics_dict = dict(zip(df['Metric'], df['Value']))
                
                # Add metadata
                metrics_dict.update({
                    'Method': method_name,
                    'Dataset': dataset_name,
                    'Clustering': clustering_method,
                    'GT_Available': gt_status == 'withGT',
                    'Slice': 'horizontal',  # For horizontal integration
                    'Dataset_Type': 'RNA_ADT' if dataset_name in ['HLN', 'HT', 'Mouse_Thymus', 'Mouse_Spleen'] else 'RNA_ATAC'
                })
                
                results_by_clustering[clustering_method].append(metrics_dict)
            else:
                # Fallback for other formats
                for _, row in df.iterrows():
                    result_dict = row.to_dict()
                    result_dict.update({
                        'Method': method_name,
                        'Dataset': dataset_name,
                        'Clustering': clustering_method,
                        'GT_Available': gt_status == 'withGT',
                        'Slice': 'horizontal',  # For horizontal integration
                        'Dataset_Type': 'RNA_ADT' if dataset_name in ['HLN', 'HT', 'Mouse_Thymus', 'Mouse_Spleen'] else 'RNA_ATAC'
                    })
                    
                    results_by_clustering[clustering_method].append(result_dict)
        
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    # Convert lists to DataFrames
    for clustering_method in results_by_clustering:
        if results_by_clustering[clustering_method]:
            results_by_clustering[clustering_method] = pd.DataFrame(results_by_clustering[clustering_method])
        else:
            results_by_clustering[clustering_method] = pd.DataFrame()
        
        print(f"  {clustering_method}: {len(results_by_clustering[clustering_method])} results")
    
    return results_by_clustering

def save_horizontal_results(final_results, output_dir):
    """
    Save horizontal integration results to files
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    for clustering_method, results in final_results.items():
        if not results:
            continue
        
        clustering_dir = os.path.join(output_dir, clustering_method)
        os.makedirs(clustering_dir, exist_ok=True)
        
        # Save detailed results
        detailed_file = os.path.join(clustering_dir, f"detailed_results_{clustering_method}.csv")
        results['detailed_results'].to_csv(detailed_file, index=False)
        print(f"Saved detailed results: {detailed_file}")
        
        # Save summary tables
        if not results['rna_adt_summary'].empty:
            rna_adt_file = os.path.join(clustering_dir, f"rna_adt_summary_{clustering_method}.csv") 
            results['rna_adt_summary'].to_csv(rna_adt_file)
            print(f"Saved RNA+ADT summary: {rna_adt_file}")
        
        if not results['rna_atac_summary'].empty:
            rna_atac_file = os.path.join(clustering_dir, f"rna_atac_summary_{clustering_method}.csv")
            results['rna_atac_summary'].to_csv(rna_atac_file)
            print(f"Saved RNA+ATAC summary: {rna_atac_file}")
        
        if not results['comprehensive_summary'].empty:
            comprehensive_file = os.path.join(clustering_dir, f"comprehensive_summary_{clustering_method}.csv")
            results['comprehensive_summary'].to_csv(comprehensive_file)
            print(f"Saved comprehensive summary: {comprehensive_file}")

def generate_horizontal_evaluation():
    """
    Generate comprehensive evaluation results for horizontal integration
    """
    
    print("="*80)
    print("SMOBench Horizontal Integration - Final Results Generation")
    print("="*80)
    
    # Setup directories
    results_dir = "/home/zhenghong/SMOBench-CLEAN/Results/evaluation/horizontal_integration"
    output_dir = "/home/zhenghong/SMOBench-CLEAN/Results/evaluation/horizontal_integration/final_results"
    
    if not os.path.exists(results_dir):
        print(f"Error: Results directory not found: {results_dir}")
        print("Please run eval_horizontal_integration.py first to generate evaluation results.")
        return
    
    # Process individual results
    results_by_clustering = process_horizontal_results(results_dir)
    
    if not any(not df.empty for df in results_by_clustering.values()):
        print("Error: No evaluation results found.")
        return
    
    # Process each clustering method
    final_results = {}
    
    for clustering_method in ['leiden', 'louvain', 'kmeans', 'mclust']:
        if clustering_method not in results_by_clustering or results_by_clustering[clustering_method].empty:
            print(f"Warning: No results found for {clustering_method}")
            continue
        
        print(f"\nProcessing {clustering_method} clustering results...")
        
        # Aggregate by dataset (horizontal integration typically has one result per dataset)
        # Convert DataFrame to list of dictionaries for aggregation
        results_list = results_by_clustering[clustering_method].to_dict('records')
        aggregated_df = aggregate_by_dataset(results_list)
        
        if aggregated_df.empty:
            print(f"Warning: No aggregated results for {clustering_method}")
            continue
        
        # Calculate normalized scores for horizontal integration
        scored_df = calculate_horizontal_normalized_scores(aggregated_df)
        
        # Create summary tables
        rna_adt_table, rna_atac_table, comprehensive_table = create_horizontal_summary_tables(scored_df, clustering_method)
        
        # Store results
        final_results[clustering_method] = {
            'detailed_results': scored_df,
            'rna_adt_summary': rna_adt_table,
            'rna_atac_summary': rna_atac_table,
            'comprehensive_summary': comprehensive_table
        }
        
        print(f"  Processed {len(scored_df)} dataset-level results")
        if not rna_adt_table.empty:
            print(f"  RNA_ADT summary: {rna_adt_table.shape[0]} methods × {rna_adt_table.shape[1]-1} datasets")
        if not rna_atac_table.empty:
            print(f"  RNA_ATAC summary: {rna_atac_table.shape[0]} methods × {rna_atac_table.shape[1]-1} datasets")
        if not comprehensive_table.empty:
            print(f"  Comprehensive summary: {comprehensive_table.shape[0]} methods × {comprehensive_table.shape[1]-1} datasets")
    
    # Save all results
    save_horizontal_results(final_results, output_dir)
    
    # Display summary tables
    print("\n" + "="*80)
    print("HORIZONTAL INTEGRATION - FINAL SUMMARY TABLES")
    print("="*80)
    
    for clustering_method in ['leiden', 'louvain', 'kmeans', 'mclust']:
        if clustering_method not in final_results:
            continue
        
        print(f"\n{clustering_method.upper()} CLUSTERING RESULTS:")
        print("-" * 50)
        
        if not final_results[clustering_method]['comprehensive_summary'].empty:
            print(f"\nComprehensive Summary (All Supported Datasets):")
            print(final_results[clustering_method]['comprehensive_summary'].to_string())
        
        if not final_results[clustering_method]['rna_adt_summary'].empty:
            print(f"\nRNA + ADT Integration:")
            print(final_results[clustering_method]['rna_adt_summary'].to_string())
        
        if not final_results[clustering_method]['rna_atac_summary'].empty:
            print(f"\nRNA + ATAC Integration:")
            print(final_results[clustering_method]['rna_atac_summary'].to_string())
    
    print(f"\n" + "="*80)
    print("Horizontal Integration Results Generation Complete!")
    print(f"All results saved to: {output_dir}")
    print("="*80)

if __name__ == "__main__":
    generate_horizontal_evaluation()