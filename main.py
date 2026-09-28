"""总入口：配置在project_config.py中修改。"""

import traceback

import project_config as cfg
from common import discover_files, save_config_snapshot
from ecg_preprocess import preprocess_all_ecg
from ecg_qc import run_ecg_qc
from ecg_finalize import finalize_all_ecg
from eda_preprocess import preprocess_all_eda
from eda_qc import run_eda_qc
from eda_finalize import finalize_all_eda


def main():
    files = discover_files()
    snapshot = save_config_snapshot()
    print("实验：{}，文件：{}，配置：{}".format(cfg.EXPERIMENT_NAME, len(files), snapshot))

    if cfg.RUN_ECG:
        print("\n[1/6] ECG初始预处理与可视化")
        preprocess_all_ecg(files)
        print("\n[2/6] ECG自动QC与可视化")
        run_ecg_qc()
        print("\n[3/6] ECG按QC区间重新处理并计算最终特征")
        finalize_all_ecg(files)

    if cfg.RUN_EDA:
        print("\n[4/6] EDA初始预处理与可视化")
        preprocess_all_eda(files)
        print("\n[5/6] EDA自动QC与可视化")
        run_eda_qc()
        print("\n[6/6] EDA按QC区间重新处理并计算最终特征")
        finalize_all_eda(files)

    print("\n完成。结果位于：{}".format(cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
