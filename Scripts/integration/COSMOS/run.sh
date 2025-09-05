#!/bin/bash

# Comprehensive COSMOS run script for all datasets
# Use relative paths for portability - script should be run from SMOBench root directory
# Cluster numbers: HLN:6 HT:6 mouse_thymus:8 mouse_spleen:5 mouse_brain:18 3M:5 misar:16

echo "Starting comprehensive COSMOS processing..."
echo "Start time: $(date)"

# Create base results directory
mkdir -p Results/adata/COSMOS Results/plot/COSMOS

# === withGT RNA+ADT Datasets ===

echo "Processing Human_Lymph_Nodes datasets..."

# Human_Lymph_Nodes A1
echo "Processing Human_Lymph_Nodes A1..."
python Scripts/integration/COSMOS/run_cosmos.py \
  --data_type 10x \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_ADT.h5ad \
  --save_path Results/adata/COSMOS/HLN/A1/COSMOS_HLN_A1.h5ad \
  --method COSMOS \
  --dataset Human_Lymph_Nodes/A1 \
  --cluster_nums 6

# Human_Lymph_Nodes D1
echo "Processing Human_Lymph_Nodes D1..."
python Scripts/integration/COSMOS/run_cosmos.py \
  --data_type 10x \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_ADT.h5ad \
  --save_path Results/adata/COSMOS/HLN/D1/COSMOS_HLN_D1.h5ad \
  --method COSMOS \
  --dataset Human_Lymph_Nodes/D1 \
  --cluster_nums 6

echo "Processing Human_Tonsils datasets..."

# Human_Tonsils S1
echo "Processing Human_Tonsils S1..."
python Scripts/integration/COSMOS/run_cosmos.py \
  --data_type 10x \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad \
  --save_path Results/adata/COSMOS/HT/S1/COSMOS_HT_S1.h5ad \
  --method COSMOS \
  --dataset Human_Tonsils/S1 \
  --cluster_nums 6

# Human_Tonsils S2
echo "Processing Human_Tonsils S2..."
python Scripts/integration/COSMOS/run_cosmos.py \
  --data_type 10x \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_ADT.h5ad \
  --save_path Results/adata/COSMOS/HT/S2/COSMOS_HT_S2.h5ad \
  --method COSMOS \
  --dataset Human_Tonsils/S2 \
  --cluster_nums 6

# Human_Tonsils S3
echo "Processing Human_Tonsils S3..."
python Scripts/integration/COSMOS/run_cosmos.py \
  --data_type 10x \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_ADT.h5ad \
  --save_path Results/adata/COSMOS/HT/S3/COSMOS_HT_S3.h5ad \
  --method COSMOS \
  --dataset Human_Tonsils/S3 \
  --cluster_nums 6

# === woGT RNA+ADT Datasets ===

echo "Processing Mouse_Thymus datasets..."

# Mouse_Thymus datasets
for thymus_id in 1 2 3 4; do
    echo "Processing Mouse_Thymus${thymus_id}..."
    python Scripts/integration/COSMOS/run_cosmos.py \
      --data_type Stereo-CITE-seq \
      --RNA_path Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus${thymus_id}/adata_RNA.h5ad \
      --ADT_path Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus${thymus_id}/adata_ADT.h5ad \
      --save_path Results/adata/COSMOS/Mouse_Thymus/Thymus${thymus_id}/COSMOS_MT_Thymus${thymus_id}.h5ad \
      --method COSMOS \
      --dataset Mouse_Thymus/Mouse_Thymus${thymus_id} \
      --cluster_nums 8
done

echo "Processing Mouse_Spleen datasets..."

# Mouse_Spleen datasets
for spleen_id in 1 2; do
    echo "Processing Mouse_Spleen${spleen_id}..."
    python Scripts/integration/COSMOS/run_cosmos.py \
      --data_type SPOTS \
      --RNA_path Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen${spleen_id}/adata_RNA.h5ad \
      --ADT_path Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen${spleen_id}/adata_ADT.h5ad \
      --save_path Results/adata/COSMOS/Mouse_Spleen/Spleen${spleen_id}/COSMOS_MS_Spleen${spleen_id}.h5ad \
      --method COSMOS \
      --dataset Mouse_Spleen/Mouse_Spleen${spleen_id} \
      --cluster_nums 5
done

# === RNA+ATAC Datasets ===

echo "Processing Mouse_Embryos RNA+ATAC datasets..."

# Check if Mouse_Embryos datasets exist
if [ -d "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1" ]; then
    for embryo_stage in E11 E13 E15 E18; do
        if [ -f "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/${embryo_stage}/adata_RNA.h5ad" ]; then
            echo "Processing Mouse_Embryos_S1 ${embryo_stage}..."
            python Scripts/integration/COSMOS/run_cosmos.py \
              --data_type MISAR \
              --RNA_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/${embryo_stage}/adata_RNA.h5ad \
              --ATAC_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/${embryo_stage}/adata_ATAC.h5ad \
              --save_path Results/adata/COSMOS/MISAR_S1/${embryo_stage}/COSMOS_MISAR_S1_${embryo_stage}.h5ad \
              --method COSMOS \
              --dataset Mouse_Embryos_S1/${embryo_stage} \
              --cluster_nums 16
        fi
    done
fi

if [ -d "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2" ]; then
    for embryo_stage in E11 E13 E15 E18; do
        if [ -f "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/${embryo_stage}/adata_RNA.h5ad" ]; then
            echo "Processing Mouse_Embryos_S2 ${embryo_stage}..."
            python Scripts/integration/COSMOS/run_cosmos.py \
              --data_type MISAR \
              --RNA_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/${embryo_stage}/adata_RNA.h5ad \
              --ATAC_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/${embryo_stage}/adata_ATAC.h5ad \
              --save_path Results/adata/COSMOS/MISAR_S2/${embryo_stage}/COSMOS_MISAR_S2_${embryo_stage}.h5ad \
              --method COSMOS \
              --dataset Mouse_Embryos_S2/${embryo_stage} \
              --cluster_nums 16
        fi
    done
fi

echo "Processing Mouse_Brain RNA+ATAC datasets..."

# Mouse_Brain datasets (woGT RNA+ATAC)
brain_types=("ATAC" "H3K4me3" "H3K27ac" "H3K27me3")
for brain_type in "${brain_types[@]}"; do
    if [ -f "Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_${brain_type}/adata_RNA.h5ad" ]; then
        echo "Processing Mouse_Brain ${brain_type}..."
        # Check if the ATAC file exists with different naming
        atac_file="Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_${brain_type}/adata_ATAC.h5ad"
        if [ ! -f "$atac_file" ]; then
            atac_file="Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_${brain_type}/adata_peaks_normalized.h5ad"
        fi
        
        if [ -f "$atac_file" ]; then
            python Scripts/integration/COSMOS/run_cosmos.py \
              --data_type Spatial-epigenome-transcriptome \
              --RNA_path Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_${brain_type}/adata_RNA.h5ad \
              --ATAC_path "$atac_file" \
              --save_path Results/adata/COSMOS/Mouse_Brain/${brain_type}/COSMOS_MB_${brain_type}.h5ad \
              --method COSMOS \
              --dataset Mouse_Brain/Mouse_Brain_${brain_type} \
              --cluster_nums 18
        fi
    fi
done

# # === 3M Simulation Dataset ===

# echo "Processing 3M Simulation dataset..."

# if [ -f "Dataset/withGT/3M_Simulation/adata_RNA.h5ad" ] && [ -f "Dataset/withGT/3M_Simulation/adata_ADT.h5ad" ]; then
#     echo "Processing 3M Simulation (RNA+ADT)..."
#     python Scripts/integration/COSMOS/run_cosmos.py \
#       --data_type simulation \
#       --RNA_path Dataset/withGT/3M_Simulation/adata_RNA.h5ad \
#       --ADT_path Dataset/withGT/3M_Simulation/adata_ADT.h5ad \
#       --save_path Results/adata/COSMOS/3M_Simulation/COSMOS_3M_Sim.h5ad \
#       --method COSMOS \
#       --dataset 3M_Simulation/Simulation \
#       --cluster_nums 5
# fi

echo "COSMOS processing completed!"
echo "End time: $(date)"

# Generate summary report
echo "=== PROCESSING SUMMARY ==="
echo "Results saved to Results/adata/COSMOS/"
echo "Plots saved to Results/plot/COSMOS/"
echo ""
echo "Processed datasets:"
find Results/adata/COSMOS -name "*.h5ad" | sort | while read file; do
    echo "  - $file"
done

echo ""
echo "Total results: $(find Results/adata/COSMOS -name "*.h5ad" | wc -l) datasets processed"