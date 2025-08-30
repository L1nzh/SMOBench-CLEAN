import os
import sys

# 添加 PRAGA 方法的根目录（包含 PRAGA/PRAGA/ 的目录）
praga_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../Methods/PRAGA/PRAGA"))
sys.path.append(praga_root)
present_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../Methods/PRESENT"))
sys.path.append(present_root)

# 添加 SMOBench 根目录（用于 PRESENT、SpatialGlue 等）
smobench_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../.."))
sys.path.append(smobench_root)

# 现在可以安全导入
from PRAGA.preprocess import clr_normalize_each_cell, pca, lsi
from PRAGA.Train_model import Train
from PRESENT.Utils import combine_BC
from SpatialGlue.preprocess import fix_seed

print("All imports successful!")