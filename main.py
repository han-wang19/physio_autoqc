"""项目入口：运行全部流程，或单独运行预处理、QC、最终计算。"""

import argparse
import traceback

import project_config as cfg
from processing.common import discover_files, save_config_snapshot
from processing.ecg_preprocess import preprocess_all_ecg
from processing.ecg_qc import run_ecg_qc
from processing.ecg_finalize import finalize_all_ecg
from processing.eda_preprocess import preprocess_all_eda
from processing.eda_qc import run_eda_qc
from processing.eda_finalize import finalize_all_eda


def _parse_args():
    """读取命令行阶段选择。"""
    parser = argparse.ArgumentParser(description="ECG/EDA preprocessing and quality control")
    parser.add_argument('--experiment', choices=tuple(cfg.EXPERIMENTS), default=cfg.EXPERIMENT_NAME)
    parser.add_argument(
        "--stage", choices=("all", "preprocess", "qc", "finalize"), default="all",
        help="all=完整流程；preprocess=仅生成NPZ；qc=读取NPZ重新质检绘图；finalize=最终计算",
    )
    return parser.parse_args()


def main(stage="all"):
    """按所选阶段处理配置中启用的 ECG 和 EDA。"""
    needs_source_files = stage in ("all", "preprocess", "finalize")
    files = discover_files() if needs_source_files else []
    snapshot = save_config_snapshot()
    print("实验：{}，阶段：{}，配置：{}".format(cfg.EXPERIMENT_NAME, stage, snapshot))
    if needs_source_files:
        print("源文件：{}".format(len(files)))

    if cfg.RUN_ECG:
        if stage in ("all", "preprocess"):
            print("\nECG | 初始预处理")
            _check_preprocessing(preprocess_all_ecg(files), 'ECG')
        if stage in ("all", "qc"):
            print("\nECG | 从预处理NPZ执行QC与绘图")
            run_ecg_qc()
        if stage in ("all", "finalize"):
            print("\nECG | 按QC区间重新处理并计算最终特征")
            finalize_all_ecg(files)

    if cfg.RUN_EDA:
        if stage in ("all", "preprocess"):
            print("\nEDA | 初始预处理")
            _check_preprocessing(preprocess_all_eda(files), 'EDA')
        if stage in ("all", "qc"):
            print("\nEDA | 从预处理NPZ执行QC与绘图")
            run_eda_qc()
        if stage in ("all", "finalize"):
            print("\nEDA | 按QC区间重新处理并计算最终特征")
            finalize_all_eda(files)

    print("\n完成。结果位于：{}".format(cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME))


def _check_preprocessing(rows, signal):
    failed = [str(row['record_id']) for row in rows if row.get('status') != 'OK']
    if failed:
        raise RuntimeError('{} preprocessing failed for {}; see preprocessing report'.format(
            signal, ', '.join(failed)))


if __name__ == "__main__":
    try:
        args = _parse_args()
        cfg.select_experiment(args.experiment)
        main(args.stage)
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
