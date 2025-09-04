#!/bin/bash

# Comprehensive SpatialGlue run script for all datasets
# Use relative paths for portability - script should be run from SMOBench root directory
# Cluster numbers: HLN:6 HT:6 mouse_thymus:8 mouse_spleen:5 mouse_brain:18 3M:5 misar:16

echo "Starting comprehensive SpatialGlue processing..."
echo "Start time: $(date)"

# Create base results directory
mkdir -p Results/adata/SpatialGlue Results/plot/SpatialGlue

# # === withGT RNA+ADT Datasets ===

# echo "Processing Human_Lymph_Nodes datasets..."

# # Human_Lymph_Nodes A1
# echo "Processing Human_Lymph_Nodes A1..."
# python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_RNA.h5ad \
#   --ADT_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_ADT.h5ad \
#   --save_path Results/adata/SpatialGlue/HLN/A1/SpatialGlue_HLN_A1.h5ad \
#   --method SpatialGlue \
#   --dataset Human_Lymph_Nodes/A1 \
#   --cluster_nums 6

# # Human_Lymph_Nodes D1
# echo "Processing Human_Lymph_Nodes D1..."
# python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_RNA.h5ad \
#   --ADT_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_ADT.h5ad \
#   --save_path Results/adata/SpatialGlue/HLN/D1/SpatialGlue_HLN_D1.h5ad \
#   --method SpatialGlue \
#   --dataset Human_Lymph_Nodes/D1 \
#   --cluster_nums 6

# echo "Processing Human_Tonsils datasets..."

# # Human_Tonsils S1
# echo "Processing Human_Tonsils S1..."
# python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad \
#   --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad \
#   --save_path Results/adata/SpatialGlue/HT/S1/SpatialGlue_HT_S1.h5ad \
#   --method SpatialGlue \
#   --dataset Human_Tonsils/S1 \
#   --cluster_nums 6

# # Human_Tonsils S2
# echo "Processing Human_Tonsils S2..."
# python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_RNA.h5ad \
#   --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_ADT.h5ad \
#   --save_path Results/adata/SpatialGlue/HT/S2/SpatialGlue_HT_S2.h5ad \
#   --method SpatialGlue \
#   --dataset Human_Tonsils/S2 \
#   --cluster_nums 6

# # Human_Tonsils S3
# echo "Processing Human_Tonsils S3..."
# python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_RNA.h5ad \
#   --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_ADT.h5ad \
#   --save_path Results/adata/SpatialGlue/HT/S3/SpatialGlue_HT_S3.h5ad \
#   --method SpatialGlue \
#   --dataset Human_Tonsils/S3 \
#   --cluster_nums 6

# # === woGT RNA+ADT Datasets ===

# echo "Processing Mouse_Thymus datasets..."

# # Mouse_Thymus datasets
# for thymus_id in 1 2 3 4; do
#     echo "Processing Mouse_Thymus${thymus_id}..."
#     python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#       --data_type Stereo-CITE-seq \
#       --RNA_path Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus${thymus_id}/adata_RNA.h5ad \
#       --ADT_path Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus${thymus_id}/adata_ADT.h5ad \
#       --save_path Results/adata/SpatialGlue/Mouse_Thymus/Thymus${thymus_id}/SpatialGlue_MT_Thymus${thymus_id}.h5ad \
#       --method SpatialGlue \
#       --dataset Mouse_Thymus/Mouse_Thymus${thymus_id} \
#       --cluster_nums 8
# done

# echo "Processing Mouse_Spleen datasets..."

# # Mouse_Spleen datasets
# for spleen_id in 1 2; do
#     echo "Processing Mouse_Spleen${spleen_id}..."
#     python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
#       --data_type SPOTS \
#       --RNA_path Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen${spleen_id}/adata_RNA.h5ad \
#       --ADT_path Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen${spleen_id}/adata_ADT.h5ad \
#       --save_path Results/adata/SpatialGlue/Mouse_Spleen/Spleen${spleen_id}/SpatialGlue_MS_Spleen${spleen_id}.h5ad \
#       --method SpatialGlue \
#       --dataset Mouse_Spleen/Mouse_Spleen${spleen_id} \
#       --cluster_nums 5
# done

# === RNA+ATAC Datasets ===

echo "Processing Mouse_Embryos RNA+ATAC datasets..."

# Check if Mouse_Embryos datasets exist
if [ -d "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1" ]; then
    for embryo_stage in E11 E13 E15 E18; do
        if [ -f "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/${embryo_stage}/adata_RNA.h5ad" ]; then
            echo "Processing Mouse_Embryos_S1 ${embryo_stage}..."
            python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
              --data_type MISAR \
              --RNA_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/${embryo_stage}/adata_RNA.h5ad \
              --ATAC_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/${embryo_stage}/adata_ATAC.h5ad \
              --save_path Results/adata/SpatialGlue/MISAR_S1/${embryo_stage}/SpatialGlue_MISAR_S1_${embryo_stage}.h5ad \
              --method SpatialGlue \
              --dataset Mouse_Embryos_S1/${embryo_stage} \
              --cluster_nums 16
        fi
    done
fi

if [ -d "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2" ]; then
    for embryo_stage in E11 E13 E15 E18; do
        if [ -f "Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/${embryo_stage}/adata_RNA.h5ad" ]; then
            echo "Processing Mouse_Embryos_S2 ${embryo_stage}..."
            python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
              --data_type MISAR \
              --RNA_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/${embryo_stage}/adata_RNA.h5ad \
              --ATAC_path Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/${embryo_stage}/adata_ATAC.h5ad \
              --save_path Results/adata/SpatialGlue/MISAR_S2/${embryo_stage}/SpatialGlue_MISAR_S2_${embryo_stage}.h5ad \
              --method SpatialGlue \
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
            python Scripts/integration/SpatialGlue/run_SpatialGlue.py \
              --data_type Spatial-epigenome-transcriptome \
              --RNA_path Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_${brain_type}/adata_RNA.h5ad \
              --ATAC_path "$atac_file" \
              --save_path Results/adata/SpatialGlue/Mouse_Brain/${brain_type}/SpatialGlue_MB_${brain_type}.h5ad \
              --method SpatialGlue \
              --dataset Mouse_Brain/Mouse_Brain_${brain_type} \
              --cluster_nums 18
        fi
    fi
done

echo "SpatialGlue processing completed!"
echo "End time: $(date)"

# Generate summary report
echo "=== PROCESSING SUMMARY ==="
echo "Results saved to Results/adata/SpatialGlue/"
echo "Plots saved to Results/plot/SpatialGlue/"
echo ""
echo "Processed datasets:"
find Results/adata/SpatialGlue -name "*.h5ad" | sort | while read file; do
    echo "  - $file"
done

echo ""
echo "Total results: $(find Results/adata/SpatialGlue -name "*.h5ad" | wc -l) datasets processed"