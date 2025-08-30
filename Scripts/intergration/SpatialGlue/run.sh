#!/bin/bash

# 聚类数HLN:6 HT:6 mouse_thymus:8 mouse_spleen:5 mouse_brain:18 3M:5 misar？
# 定义结果保存根目录的绝对路径
RESULT_ROOT="/home/users/nus/e1503317/projects/dmeng/zhlin/SMOBench/Results/adata/SpatialGlue"

## Human_Lymph_Nodes (withGT RNA+ADT)
# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path ../../../Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_RNA.h5ad \
#   --ADT_path ../../../Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_ADT.h5ad \
#   --save_path "${RESULT_ROOT}/HLN/A1/SpatialGlue_HLN_A1.h5ad" \
#   --method SpatialGlue \
#   --dataset Human_Lymph_Nodes/A1\
#   --cluster_nums 6

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path ../../../Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_RNA.h5ad \
#   --ADT_path ../../../Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_ADT.h5ad \
#   --save_path "${RESULT_ROOT}/HLN/D1/SpatialGlue_HLN_D1.h5ad" \
#   --method SpatialGlue \
#   --dataset Human_Lymph_Nodes/D1\
#   --cluster_nums 6

## Human_Tonsils (withGT RNA+ADT)
# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path ../../../Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad \
#   --ADT_path ../../../Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad \
#   --save_path "${RESULT_ROOT}/HT/S1/SpatialGlue_HT_S1.h5ad" \
#   --method SpatialGlue \
#   --dataset Human_Tonsils/S1 \
#   --cluster_nums 6

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path ../../../Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_RNA.h5ad \
#   --ADT_path ../../../Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_ADT.h5ad \
#   --save_path "${RESULT_ROOT}/HT/S2/SpatialGlue_HT_S2.h5ad" \
#   --method SpatialGlue \
#   --dataset Human_Tonsils/S2 \
#   --cluster_nums 6

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type 10x \
#   --RNA_path ../../../Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_RNA.h5ad \
#   --ADT_path ../../../Dataset/withGT/RNA_ADT/Human_Tonsils/S3/adata_ADT.h5ad \
#   --save_path "${RESULT_ROOT}/HT/S3/SpatialGlue_HT_S3.h5ad" \
#   --method SpatialGlue \
#   --dataset Human_Tonsils/S3 \
#   --cluster_nums 6

# # Mouse_Embryos_S1 (withGT RNA+ATAC)
# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E11/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E11/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S1/E11/SpatialGlue_MISAR_S1_E11.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S1/E11

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E13/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E13/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S1/E13/SpatialGlue_MISAR_S1_E13.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S1/E13

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E15/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E15/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S1/E15/SpatialGlue_MISAR_S1_E15.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S1/E15

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E18/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S1/E18/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S1/E18/SpatialGlue_MISAR_S1_E18.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S1/E18

# # Mouse_Embryos_S2 (withGT RNA+ATAC)
# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E11/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E11/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S2/E11/SpatialGlue_MISAR_S2_E11.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S2/E11

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E13/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E13/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S2/E13/SpatialGlue_MISAR_S2_E13.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S2/E13

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E15/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E15/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S2/E15/SpatialGlue_MISAR_S2_E15.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S2/E15

# python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
#   --data_type MISAR \
#   --RNA_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E18/adata_RNA.h5ad \
#   --ATAC_path ../../../Dataset/withGT/RNA_ATAC/Mouse_Embryos_S2/E18/adata_ATAC.h5ad \
#   --save_path "${RESULT_ROOT}/MISAR_S2/E18/SpatialGlue_MISAR_S2_E18.h5ad" \
#   --method SpatialGlue \
#   --dataset Mouse_Embryos_S2/E18

# Mouse_Spleen (woGT RNA+ADT)
python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type SPOTS \
  --RNA_path ../../../Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen1/adata_RNA.h5ad \
  --ADT_path ../../../Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Spleen/Spleen1/SpatialGlue_MS_Spleen1.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Spleen/Mouse_Spleen1 \
  --cluster_nums 5

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type SPOTS \
  --RNA_path ../../../Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen2/adata_RNA.h5ad \
  --ADT_path ../../../Dataset/woGT/RNA_ADT/Mouse_Spleen/Mouse_Spleen2/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Spleen/Spleen2/SpatialGlue_MS_Spleen2.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Spleen/Mouse_Spleen2 \
  --cluster_nums 5


# Mouse_Thymus (woGT RNA+ADT)
python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Stereo-CITE-seq \
  --RNA_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus1/adata_RNA.h5ad \
  --ADT_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Thymus/Thymus1/SpatialGlue_MT_Thymus1.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Thymus/Mouse_Thymus1 \
  --cluster_nums 8

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Stereo-CITE-seq \
  --RNA_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus2/adata_RNA.h5ad \
  --ADT_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus2/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Thymus/Thymus2/SpatialGlue_MT_Thymus2.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Thymus/Mouse_Thymus2 \
  --cluster_nums 8

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Stereo-CITE-seq \
  --RNA_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus3/adata_RNA.h5ad \
  --ADT_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus3/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Thymus/Thymus3/SpatialGlue_MT_Thymus3.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Thymus/Mouse_Thymus3 \
  --cluster_nums 8

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Stereo-CITE-seq \
  --RNA_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus4/adata_RNA.h5ad \
  --ADT_path ../../../Dataset/woGT/RNA_ADT/Mouse_Thymus/Mouse_Thymus4/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Thymus/Thymus4/SpatialGlue_MT_Thymus4.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Thymus/Mouse_Thymus4 \
  --cluster_nums 8

# Mouse_Brain (woGT RNA+ATAC)
python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Spatial-epigenome-transcriptome \
  --RNA_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_RNA.h5ad \
  --ATAC_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_ATAC/adata_peaks_normalized.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Brain/ATAC/SpatialGlue_MB_ATAC.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Brain/Mouse_Brain_ATAC \
  --cluster_nums 18

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Spatial-epigenome-transcriptome \
  --RNA_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K4me3/adata_RNA.h5ad \
  --ATAC_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K4me3/adata_peaks_normalized.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Brain/H3K4me3/SpatialGlue_MB_H3K4me3.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Brain/Mouse_Brain_H3K4me3 \
  --cluster_nums 18

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Spatial-epigenome-transcriptome \
  --RNA_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K27ac/adata_RNA.h5ad \
  --ATAC_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K27ac/adata_peaks_normalized.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Brain/H3K27ac/SpatialGlue_MB_H3K27ac.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Brain/Mouse_Brain_H3K27ac \
  --cluster_nums 18

python ../../../Scripts/intergration/SpatialGlue/run_SpatialGlue.py \
  --data_type Spatial-epigenome-transcriptome \
  --RNA_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K27me3/adata_RNA.h5ad \
  --ATAC_path ../../../Dataset/woGT/RNA_ATAC/Mouse_Brain/Mouse_Brain_H3K27me3/adata_peaks_normalized.h5ad \
  --save_path "${RESULT_ROOT}/Mouse_Brain/H3K27me3/SpatialGlue_MB_H3K27me3.h5ad" \
  --method SpatialGlue \
  --dataset Mouse_Brain/Mouse_Brain_H3K27me3 \
  --cluster_nums 18