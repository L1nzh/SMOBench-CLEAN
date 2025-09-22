#!/usr/bin/env python3
"""
Fix Data Alignment Issues in Vertical Integration Results

This script identifies and fixes datasets with mismatched array lengths
between embeddings, ground truth labels, predictions, and spatial coordinates.

Usage:
    python fix_data_alignment.py [--dry-run] [--method METHOD_NAME]
    
    --dry-run : Show what would be fixed without making changes
    --method  : Fix only specific method (e.g., PRAGA, SpaMV)
"""

import os
import sys
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def check_data_alignment(adata_path, method_name):
    """
    Check if data arrays are properly aligned
    
    Returns:
    --------
    dict : Status and details about alignment issues
    """
    try:
        adata = sc.read_h5ad(adata_path)
        
        # Get array lengths
        n_cells = adata.n_obs
        
        # Find embeddings
        embedding_key = None
        for key in [method_name, 'X_integrated', 'X_emb', 'embeddings']:
            if key in adata.obsm:
                embedding_key = key
                break
        
        if embedding_key is None:
            return {'status': 'no_embeddings', 'details': 'No embeddings found'}
        
        n_embeddings = adata.obsm[embedding_key].shape[0]
        
        # Check spatial coordinates
        if 'spatial' not in adata.obsm:
            return {'status': 'no_spatial', 'details': 'No spatial coordinates found'}
        
        n_spatial = adata.obsm['spatial'].shape[0]
        
        # Check ground truth (if available)
        n_gt = None
        if 'Spatial_Label' in adata.obs:
            n_gt = len(adata.obs['Spatial_Label'])
        
        # Check clustering results
        clustering_methods = ['leiden', 'louvain', 'kmeans', 'mclust']
        clustering_lengths = {}
        for method in clustering_methods:
            if method in adata.obs:
                clustering_lengths[method] = len(adata.obs[method])
        
        # Analyze alignment
        lengths = {
            'cells': n_cells,
            'embeddings': n_embeddings,
            'spatial': n_spatial
        }
        
        if n_gt is not None:
            lengths['ground_truth'] = n_gt
        
        lengths.update({f'clustering_{k}': v for k, v in clustering_lengths.items()})
        
        # Check if all are equal
        unique_lengths = set(lengths.values())
        
        if len(unique_lengths) == 1:
            return {'status': 'aligned', 'details': lengths}
        else:
            return {
                'status': 'misaligned',
                'details': lengths,
                'embedding_key': embedding_key
            }
            
    except Exception as e:
        return {'status': 'error', 'details': str(e)}

def fix_data_alignment(adata_path, method_name, dry_run=False):
    """
    Fix data alignment by truncating to minimum common length
    
    Returns:
    --------
    dict : Fix results
    """
    try:
        print(f"Processing: {adata_path}")
        
        # Check current alignment
        alignment_check = check_data_alignment(adata_path, method_name)
        
        if alignment_check['status'] != 'misaligned':
            return {'status': 'no_fix_needed', 'details': alignment_check}
        
        # Load data
        adata = sc.read_h5ad(adata_path)
        lengths = alignment_check['details']
        embedding_key = alignment_check['embedding_key']
        
        print(f"  Current lengths: {lengths}")
        
        # Find minimum length (excluding 'cells' which is just adata.n_obs)
        relevant_lengths = {k: v for k, v in lengths.items() if k != 'cells'}
        min_length = min(relevant_lengths.values())
        
        print(f"  Truncating to minimum length: {min_length}")
        
        if dry_run:
            return {
                'status': 'would_fix',
                'original_lengths': lengths,
                'target_length': min_length
            }
        
        # Create backup
        backup_path = adata_path + '.backup'
        if not os.path.exists(backup_path):
            print(f"  Creating backup: {backup_path}")
            adata.write(backup_path)
        
        # Truncate all arrays to minimum length
        adata_fixed = adata[:min_length].copy()
        
        # Ensure all obsm arrays are truncated
        for key in adata_fixed.obsm.keys():
            if adata_fixed.obsm[key].shape[0] > min_length:
                adata_fixed.obsm[key] = adata_fixed.obsm[key][:min_length]
        
        # Verify fix
        new_check = check_data_alignment_direct(adata_fixed, method_name)
        
        if new_check['status'] == 'aligned':
            # Save fixed data
            print(f"  Saving fixed data to: {adata_path}")
            adata_fixed.write(adata_path)
            
            return {
                'status': 'fixed',
                'original_lengths': lengths,
                'new_length': min_length,
                'backup_created': backup_path
            }
        else:
            return {
                'status': 'fix_failed',
                'details': new_check
            }
            
    except Exception as e:
        return {'status': 'error', 'details': str(e)}

def check_data_alignment_direct(adata, method_name):
    """Check alignment directly on loaded AnnData object"""
    try:
        n_cells = adata.n_obs
        
        embedding_key = None
        for key in [method_name, 'X_integrated', 'X_emb', 'embeddings']:
            if key in adata.obsm:
                embedding_key = key
                break
        
        if embedding_key is None:
            return {'status': 'no_embeddings'}
        
        n_embeddings = adata.obsm[embedding_key].shape[0]
        n_spatial = adata.obsm['spatial'].shape[0] if 'spatial' in adata.obsm else 0
        
        lengths = {
            'cells': n_cells,
            'embeddings': n_embeddings,
            'spatial': n_spatial
        }
        
        if 'Spatial_Label' in adata.obs:
            lengths['ground_truth'] = len(adata.obs['Spatial_Label'])
        
        unique_lengths = set(lengths.values())
        
        if len(unique_lengths) == 1:
            return {'status': 'aligned', 'details': lengths}
        else:
            return {'status': 'misaligned', 'details': lengths}
            
    except Exception as e:
        return {'status': 'error', 'details': str(e)}

def scan_and_fix_datasets(results_dir, target_method=None, dry_run=False):
    """
    Scan and fix all datasets with alignment issues
    """
    print("="*80)
    print("SMOBench Data Alignment Fixer")
    print("="*80)
    
    if dry_run:
        print("DRY RUN MODE - No changes will be made")
    
    # Methods to check
    methods = ['CANDIES', 'COSMOS', 'PRAGA', 'PRESENT', 'SpaMV', 'SpaMosaic', 'SpatialGlue']
    if target_method:
        methods = [target_method]
    
    # Track results
    total_files = 0
    misaligned_files = 0
    fixed_files = 0
    failed_fixes = 0
    
    fix_log = []
    
    for method_name in methods:
        method_dir = os.path.join(results_dir, method_name)
        if not os.path.exists(method_dir):
            print(f"Warning: Method directory not found: {method_dir}")
            continue
        
        print(f"\n--- Checking Method: {method_name} ---")
        
        # Find all .h5ad files for this method
        h5ad_files = list(Path(method_dir).rglob("*.h5ad"))
        
        for h5ad_file in h5ad_files:
            # Skip backup files
            if h5ad_file.name.endswith('.backup'):
                continue
                
            total_files += 1
            file_path = str(h5ad_file)
            
            # Check alignment
            alignment_check = check_data_alignment(file_path, method_name)
            
            if alignment_check['status'] == 'misaligned':
                misaligned_files += 1
                
                print(f"\nMisaligned: {file_path}")
                print(f"  Lengths: {alignment_check['details']}")
                
                # Attempt fix
                fix_result = fix_data_alignment(file_path, method_name, dry_run)
                
                if fix_result['status'] == 'fixed':
                    fixed_files += 1
                    print(f"  FIXED: Truncated to {fix_result['new_length']} cells")
                elif fix_result['status'] == 'would_fix':
                    print(f"  → WOULD FIX: Truncate to {fix_result['target_length']} cells")
                else:
                    failed_fixes += 1
                    print(f"  FAILED: {fix_result.get('details', 'Unknown error')}")
                
                fix_log.append({
                    'file': file_path,
                    'method': method_name,
                    'status': fix_result['status'],
                    'details': fix_result
                })
            
            elif alignment_check['status'] == 'aligned':
                print(f"Aligned: {h5ad_file.name}")
            else:
                print(f"? Issue: {h5ad_file.name} - {alignment_check['status']}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total files checked: {total_files}")
    print(f"Misaligned files found: {misaligned_files}")
    
    if dry_run:
        print(f"Files that would be fixed: {len([r for r in fix_log if r['status'] == 'would_fix'])}")
    else:
        print(f"Files successfully fixed: {fixed_files}")
        print(f"Files failed to fix: {failed_fixes}")
    
    if fix_log:
        print(f"\nDetailed log saved to: alignment_fix_log.txt")
        
        # Save detailed log
        with open('alignment_fix_log.txt', 'w') as f:
            f.write("SMOBench Data Alignment Fix Log\n")
            f.write("="*50 + "\n\n")
            
            for entry in fix_log:
                f.write(f"File: {entry['file']}\n")
                f.write(f"Method: {entry['method']}\n")
                f.write(f"Status: {entry['status']}\n")
                f.write(f"Details: {entry['details']}\n")
                f.write("-" * 30 + "\n")
    
    return {
        'total_files': total_files,
        'misaligned_files': misaligned_files,
        'fixed_files': fixed_files,
        'failed_fixes': failed_fixes,
        'fix_log': fix_log
    }

def main():
    parser = argparse.ArgumentParser(description='Fix data alignment issues in SMOBench results')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be fixed without making changes')
    parser.add_argument('--method', type=str, help='Fix only specific method (e.g., PRAGA, SpaMV)')
    parser.add_argument('--results-dir', type=str, default='/home/zhenghong/SMOBench-CLEAN/Results/adata/vertical_integration',
                       help='Path to results directory')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results_dir):
        print(f"Error: Results directory not found: {args.results_dir}")
        sys.exit(1)
    
    # Run the scan and fix
    results = scan_and_fix_datasets(
        results_dir=args.results_dir,
        target_method=args.method,
        dry_run=args.dry_run
    )
    
    print("\nData alignment fix completed!")

if __name__ == "__main__":
    main()