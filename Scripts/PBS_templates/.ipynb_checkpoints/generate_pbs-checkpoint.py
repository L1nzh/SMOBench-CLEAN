#!/usr/bin/env python3

import csv
import os
from string import Template

def generate_pbs_scripts(template_path, tasks_csv, output_dir="./"):
    """
    Generate PBS scripts from template and task configuration CSV
    
    Parameters:
        template_path: Path to PBS template file
        tasks_csv: Path to CSV file with task parameters
        output_dir: Directory to save generated PBS scripts
    """
    
    # Read PBS template
    with open(template_path, "r") as f:
        template = Template(f.read())
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read tasks from CSV and generate scripts
    generated_count = 0
    with open(tasks_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Handle empty values (set to empty string for optional parameters)
            for key, value in row.items():
                if not value or value.lower() == 'none':
                    row[key] = ""
            
            # Generate PBS script content
            pbs_content = template.safe_substitute(row)
            
            # Create PBS filename
            pbs_filename = f"run_{row['TASK_NAME']}.pbs"
            pbs_path = os.path.join(output_dir, pbs_filename)
            
            # Write PBS script
            with open(pbs_path, "w") as pbs_file:
                pbs_file.write(pbs_content)
                
            print(f"Generated: {pbs_filename}")
            generated_count += 1
    
    print(f"\nGenerated {generated_count} PBS scripts in {output_dir}")
    print("To submit all jobs, run:")
    print(f"cd {output_dir} && for pbs in run_*.pbs; do qsub $pbs; done")


def create_sample_tasks_csv(output_path="tasks_sample.csv"):
    """Create a sample tasks CSV file for reference"""
    
    sample_tasks = [
        {
            'TASK_NAME': 'SpatialGlue_HT_S1',
            'SCRIPT_PATH': 'Scripts/integration/SpatialGlue/run_SpatialGlue.py',
            'DATA_TYPE': '10x',
            'RNA_PATH': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad',
            'ADT_PATH': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad',
            'ATAC_PATH': '',
            'SAVE_PATH': 'Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1.h5ad',
            'LOG_PATH': 'Results/logs/SpatialGlue_HT_S1.log',
            'DEVICE': 'cuda:0',
            'SEED': '2024',
            'N_CLUSTERS': '6'
        },
        {
            'TASK_NAME': 'SpatialGlue_MB_ATAC',
            'SCRIPT_PATH': 'Scripts/integration/SpatialGlue/run_SpatialGlue.py',
            'DATA_TYPE': 'Spatial-epigenome-transcriptome',
            'RNA_PATH': 'Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_RNA.h5ad',
            'ADT_PATH': '',
            'ATAC_PATH': 'Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_peaks_normalized.h5ad',
            'SAVE_PATH': 'Results/adata/SpatialGlue/MB/ATAC/SpatialGlue_MB_ATAC.h5ad',
            'LOG_PATH': 'Results/logs/SpatialGlue_MB_ATAC.log',
            'DEVICE': 'cuda:0',
            'SEED': '2024',
            'N_CLUSTERS': '8'
        }
    ]
    
    # Write sample CSV
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = sample_tasks[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for task in sample_tasks:
            writer.writerow(task)
    
    print(f"Created sample tasks CSV: {output_path}")


def generate_full_tasks_csv(output_path="tasks_full.csv"):
    """Generate complete tasks CSV for all SMOBench datasets and methods"""
    
    tasks = []
    
    # Dataset definitions with their parameters
    datasets = {
        # With Ground Truth - RNA+ADT
        'HT_S1': {
            'data_type': '10x',
            'rna_path': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad',
            'adt_path': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad',
            'atac_path': '',
            'n_clusters': '6'
        },
        'HT_S2': {
            'data_type': '10x',
            'rna_path': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_RNA.h5ad',
            'adt_path': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_ADT.h5ad',
            'atac_path': '',
            'n_clusters': '6'
        },
        'HT_S3': {
            'data_type': '10x', 
            'rna_path': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_RNA.h5ad',
            'adt_path': 'Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_ADT.h5ad',
            'atac_path': '',
            'n_clusters': '6'
        },
        'HLN_A1': {
            'data_type': '10x',
            'rna_path': 'Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_RNA.h5ad',
            'adt_path': 'Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_ADT.h5ad',
            'atac_path': '',
            'n_clusters': '7'
        },
        'HLN_D1': {
            'data_type': '10x',
            'rna_path': 'Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_RNA.h5ad',
            'adt_path': 'Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_ADT.h5ad',
            'atac_path': '',
            'n_clusters': '7'
        },
        # With Ground Truth - RNA+ATAC
        'ME_S1_E11': {
            'data_type': 'Spatial-epigenome-transcriptome',
            'rna_path': 'Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E11/adata_RNA.h5ad',
            'adt_path': '',
            'atac_path': 'Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E11/adata_ATAC.h5ad',
            'n_clusters': '8'
        },
        # Without Ground Truth - RNA+ADT
        'MS_Spleen1': {
            'data_type': '10x',
            'rna_path': 'Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen1/adata_RNA.h5ad',
            'adt_path': 'Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen1/adata_ADT.h5ad',
            'atac_path': '',
            'n_clusters': '8'
        },
        # Without Ground Truth - RNA+ATAC  
        'MB_ATAC': {
            'data_type': 'Spatial-epigenome-transcriptome',
            'rna_path': 'Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_RNA.h5ad',
            'adt_path': '',
            'atac_path': 'Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_peaks_normalized.h5ad',
            'n_clusters': '10'
        }
    }
    
    # Method definitions
    methods = {
        'SpatialGlue': {
            'script_path': 'Scripts/integration/SpatialGlue/run_SpatialGlue.py'
        },
        'SpaMosaic': {
            'script_path': 'Scripts/integration/SpaMosaic/run_SpaMosaic.py'
        },
        'PRAGA': {
            'script_path': 'Scripts/integration/PRAGA/run_PRAGA.py'
        }
    }
    
    # Generate tasks for each method-dataset combination
    for method_name, method_info in methods.items():
        for dataset_name, dataset_info in datasets.items():
            task = {
                'TASK_NAME': f"{method_name}_{dataset_name}",
                'SCRIPT_PATH': method_info['script_path'],
                'DATA_TYPE': dataset_info['data_type'],
                'RNA_PATH': dataset_info['rna_path'],
                'ADT_PATH': dataset_info['adt_path'],
                'ATAC_PATH': dataset_info['atac_path'],
                'SAVE_PATH': f"Results/adata/{method_name}/{dataset_name.split('_')[0]}/{dataset_name.split('_')[-1]}/{method_name}_{dataset_name}.h5ad",
                'LOG_PATH': f"Results/logs/{method_name}_{dataset_name}.log",
                'DEVICE': 'cuda:0',
                'SEED': '2024',
                'N_CLUSTERS': dataset_info['n_clusters']
            }
            tasks.append(task)
    
    # Write full CSV
    with open(output_path, 'w', newline='') as csvfile:
        if tasks:
            fieldnames = tasks[0].keys()
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for task in tasks:
                writer.writerow(task)
    
    print(f"Generated full tasks CSV with {len(tasks)} tasks: {output_path}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate PBS scripts from template and CSV")
    parser.add_argument("--template", default="Scripts/PBS_templates/template.pbs", help="PBS template file")
    parser.add_argument("--tasks", help="Tasks CSV file")
    parser.add_argument("--output-dir", default="Scripts/PBS_templates/generated/", help="Output directory for PBS scripts")
    parser.add_argument("--create-sample", action="store_true", help="Create sample tasks CSV")
    parser.add_argument("--create-full", action="store_true", help="Create full tasks CSV")
    
    args = parser.parse_args()
    
    if args.create_sample:
        create_sample_tasks_csv()
    elif args.create_full:
        generate_full_tasks_csv()
    elif args.tasks:
        generate_pbs_scripts(args.template, args.tasks, args.output_dir)
    else:
        print("Please specify --tasks CSV file, --create-sample, or --create-full")
        print("\nExamples:")
        print("python generate_pbs.py --create-sample")
        print("python generate_pbs.py --create-full") 
        print("python generate_pbs.py --tasks tasks_full.csv")