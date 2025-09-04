#!/bin/bash

# 聚类数HLN:6 HT:6 mouse_thymus:8 mouse_spleen:5 mouse_brain:18 3M:5 misar:16
# Use relative paths for portability - script should be run from SMOBench root directory
RESULT_ROOT="Results/adata/SpaMV"

## Human_Lymph_Nodes (withGT RNA+ADT)
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/A1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/HLN/A1/SpaMV_HLN_A1.h5ad" \
  --method SpaMV \
  --dataset Human_Lymph_Nodes/A1 \
  --cluster_nums 6

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Lymph_Nodes/D1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/HLN/D1/SpaMV_HLN_D1.h5ad" \
  --method SpaMV \
  --dataset Human_Lymph_Nodes/D1 \
  --cluster_nums 6

## Human_Tonsils (withGT RNA+ADT)
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/HT/S1/SpaMV_HT_S1.h5ad" \
  --method SpaMV \
  --dataset Human_Tonsils/S1 \
  --cluster_nums 6

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/Human_Tonsils/S2/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/HT/S2/SpaMV_HT_S2.h5ad" \
  --method SpaMV \
  --dataset Human_Tonsils/S2 \
  --cluster_nums 6

## mouse_thymus (withGT RNA+ADT)
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/mouse_thymus/thymus1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/mouse_thymus/thymus1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/mouse_thymus/thymus1/SpaMV_mouse_thymus_thymus1.h5ad" \
  --method SpaMV \
  --dataset mouse_thymus/thymus1 \
  --cluster_nums 8

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/mouse_thymus/thymus2/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/mouse_thymus/thymus2/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/mouse_thymus/thymus2/SpaMV_mouse_thymus_thymus2.h5ad" \
  --method SpaMV \
  --dataset mouse_thymus/thymus2 \
  --cluster_nums 8

## mouse_spleen (withGT RNA+ADT)
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/mouse_spleen/spleen1/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/mouse_spleen/spleen1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/mouse_spleen/spleen1/SpaMV_mouse_spleen_spleen1.h5ad" \
  --method SpaMV \
  --dataset mouse_spleen/spleen1 \
  --cluster_nums 5

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ADT/mouse_spleen/spleen2/adata_RNA.h5ad \
  --ADT_path Dataset/withGT/RNA_ADT/mouse_spleen/spleen2/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/mouse_spleen/spleen2/SpaMV_mouse_spleen_spleen2.h5ad" \
  --method SpaMV \
  --dataset mouse_spleen/spleen2 \
  --cluster_nums 5

## mouse_brain (withGT RNA+ATAC)  
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ATAC/mouse_brain/E18_Rep1/adata_RNA.h5ad \
  --ATAC_path Dataset/withGT/RNA_ATAC/mouse_brain/E18_Rep1/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/mouse_brain/E18_Rep1/SpaMV_mouse_brain_E18_Rep1.h5ad" \
  --method SpaMV \
  --dataset mouse_brain/E18_Rep1 \
  --cluster_nums 18

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/withGT/RNA_ATAC/mouse_brain/E18_Rep2/adata_RNA.h5ad \
  --ATAC_path Dataset/withGT/RNA_ATAC/mouse_brain/E18_Rep2/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/mouse_brain/E18_Rep2/SpaMV_mouse_brain_E18_Rep2.h5ad" \
  --method SpaMV \
  --dataset mouse_brain/E18_Rep2 \
  --cluster_nums 18

## 3M datasets (woGT)
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/3M_1/adata_RNA.h5ad \
  --ADT_path Dataset/woGT/3M_1/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/3M_1/SpaMV_3M_1.h5ad" \
  --method SpaMV \
  --dataset 3M/1 \
  --cluster_nums 5

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/3M_2/adata_RNA.h5ad \
  --ADT_path Dataset/woGT/3M_2/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/3M_2/SpaMV_3M_2.h5ad" \
  --method SpaMV \
  --dataset 3M/2 \
  --cluster_nums 5

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/3M_3/adata_RNA.h5ad \
  --ADT_path Dataset/woGT/3M_3/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/3M_3/SpaMV_3M_3.h5ad" \
  --method SpaMV \
  --dataset 3M/3 \
  --cluster_nums 5

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/3M_4/adata_RNA.h5ad \
  --ADT_path Dataset/woGT/3M_4/adata_ADT.h5ad \
  --save_path "${RESULT_ROOT}/3M_4/SpaMV_3M_4.h5ad" \
  --method SpaMV \
  --dataset 3M/4 \
  --cluster_nums 5

## MISAR datasets (woGT RNA+ATAC)
python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_1/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_1/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_1/SpaMV_MISAR_misar_1.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_1 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_2/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_2/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_2/SpaMV_MISAR_misar_2.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_2 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_3/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_3/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_3/SpaMV_MISAR_misar_3.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_3 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_4/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_4/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_4/SpaMV_MISAR_misar_4.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_4 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_5/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_5/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_5/SpaMV_MISAR_misar_5.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_5 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_6/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_6/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_6/SpaMV_MISAR_misar_6.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_6 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_7/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_7/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_7/SpaMV_MISAR_misar_7.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_7 \
  --cluster_nums 16

python Scripts/integration/SpaMV/run_SpaMV.py \
  --RNA_path Dataset/woGT/MISAR/misar_8/adata_RNA.h5ad \
  --ATAC_path Dataset/woGT/MISAR/misar_8/adata_ATAC.h5ad \
  --save_path "${RESULT_ROOT}/MISAR/misar_8/SpaMV_MISAR_misar_8.h5ad" \
  --method SpaMV \
  --dataset MISAR/misar_8 \
  --cluster_nums 16