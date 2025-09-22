# SMOBench Horizontal Integration 架构设计

## 📁 新的项目结构

```
SMOBench-CLEAN/
├── Dataset/                           # 数据集目录
│   ├── withGT/                       # 带Ground Truth的数据集
│   │   ├── RNA_ADT/                  # 垂直整合原始数据
│   │   │   ├── Human_Lymph_Nodes/
│   │   │   ├── Human_Tonsils/
│   │   │   └── ...
│   │   ├── RNA_ATAC/
│   │   │   ├── MISAR/
│   │   │   └── ...
│   │   └── Fusion/                   # 🆕 水平整合输入数据
│   │       ├── RNA_ADT/
│   │       │   ├── HLN_Fusion/       # HLN A1+D1合并的RNA和ADT
│   │       │   ├── HT_Fusion/        # HT S1+S2+S3合并的RNA和ADT
│   │       │   └── ...
│   │       └── RNA_ATAC/
│   │           ├── MISAR_S1_Fusion/  # MISAR S1 E11+E13+E15+E18合并
│   │           ├── MISAR_S2_Fusion/
│   │           └── ...
│   ├── woGT/                         # 无Ground Truth的数据集
│   │   ├── RNA_ADT/                  # 垂直整合原始数据
│   │   ├── RNA_ATAC/
│   │   └── Fusion/                   # 🆕 水平整合输入数据
│   │       ├── RNA_ADT/
│   │       │   ├── MT_Fusion/        # Mouse Thymus合并
│   │       │   └── MS_Fusion/        # Mouse Spleen合并  
│   │       └── RNA_ATAC/
│   │           └── MB_Fusion/        # Mouse Brain合并
│   └── data_info/                    # 数据集信息文件
├── Scripts/                          # 执行脚本目录
│   ├── vertical_integration/         # 🔄 垂直整合脚本 (原integration)
│   │   ├── SpatialGlue/
│   │   ├── SpaMosaic/
│   │   ├── PRESENT/
│   │   ├── COSMOS/
│   │   └── ...
│   ├── horizontal_integration/       # 🆕 水平整合脚本
│   │   ├── SpatialGlue/
│   │   ├── SpaMosaic/
│   │   ├── PRESENT/
│   │   ├── COSMOS/
│   │   └── ...
│   ├── data_preparation/             # 🆕 数据准备脚本
│   │   ├── generate_fusion_data.py   # 生成Fusion数据
│   │   └── validate_fusion_data.py   # 验证Fusion数据
│   ├── evaluation/                   # 评估脚本
│   └── visualization/                # 可视化脚本
├── Results/                          # 结果目录
│   ├── adata/                        # AnnData结果
│   │   ├── vertical_integration/     # 🔄 垂直整合结果 (原有结果)
│   │   │   ├── SpatialGlue/
│   │   │   ├── SpaMosaic/
│   │   │   ├── PRESENT/
│   │   │   ├── COSMOS/
│   │   │   └── ...
│   │   └── horizontal_integration/   # 🆕 水平整合结果
│   │       ├── SpatialGlue/
│   │       ├── SpaMosaic/
│   │       ├── PRESENT/
│   │       ├── COSMOS/
│   │       └── ...
│   └── plot/                        # 可视化结果
│       ├── vertical_integration/
│       └── horizontal_integration/   # 🆕
├── Methods/                          # 集成方法 (保持不变)
├── Utils/                            # 工具函数 (保持不变)
└── docs/                             # 文档 (保持不变)
```

## 🔄 数据流程

### 1. 垂直整合 (Vertical Integration)
```
原始数据 → 方法整合 → 聚类 → 评估
Dataset/withGT/RNA_ADT/HLN/A1/ → Scripts/vertical_integration/SpatialGlue/ → Results/adata/vertical_integration/SpatialGlue/
```

### 2. 水平整合 (Horizontal Integration)  
```
Step 1: 数据合并
Dataset/withGT/RNA_ADT/HLN/A1/ + Dataset/withGT/RNA_ADT/HLN/D1/ → Dataset/withGT/Fusion/RNA_ADT/HLN_Fusion/

Step 2: 方法整合
Dataset/withGT/Fusion/RNA_ADT/HLN_Fusion/ → Scripts/horizontal_integration/SpatialGlue/ → Results/adata/horizontal_integration/SpatialGlue/
```

## 📋 Fusion数据生成规则

### withGT数据集
- **HLN_Fusion**: A1 + D1 → 12个聚类
- **HT_Fusion**: S1 + S2 + S3 → 5个聚类  
- **MISAR_S1_Fusion**: E11 + E13 + E15 + E18 → 14个聚类
- **MISAR_S2_Fusion**: E11 + E13 + E15 + E18 → 16个聚类

### woGT数据集
- **MT_Fusion**: Thymus1 + Thymus2 + Thymus3 + Thymus4 → 8个聚类
- **MS_Fusion**: Spleen1 + Spleen2 → 5个聚类
- **MB_Fusion**: Brain_ATAC + Brain_H3K4me3 + Brain_H3K27ac + Brain_H3K27me3 → 18个聚类

## 🔧 实施步骤

1. **重组现有结构**: 移动integration → vertical_integration
2. **创建Fusion数据**: 合并RNA和ADT/ATAC数据
3. **开发horizontal integration脚本**: 基于现有vertical脚本修改
4. **更新评估框架**: 支持horizontal结果评估
5. **保持向后兼容**: 确保现有功能正常运行