"""ECG和EDA共用的文件、输出和绘图辅助函数。"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import project_config as cfg


def experiment_dir():
    """返回当前实验的输出根目录。"""
    return cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME


def signal_dir(signal):
    """返回并创建 ECG 或 EDA 的输出目录。"""
    path = experiment_dir() / signal.lower()
    path.mkdir(parents=True, exist_ok=True)
    return path


def output_path(signal, name, root=None):
    """Resolve ECG/EDA deliverables and internal storage consistently."""
    root = Path(root) if root is not None else signal_dir(signal)
    if signal in ('ecg', 'eda'):
        if name in ('preprocessed', 'final', 'manual_rpeaks', 'excluded_artifacts'):
            (root / 'internal').mkdir(parents=True, exist_ok=True)
            return root / 'internal' / name
        if name in ('reviewed_fif', 'final_ecg_hrv_features.csv',
                    'ecg_valid_segments.csv', 'ecg_rr_intervals.csv', 'final_eda_features.csv'):
            (root / 'results').mkdir(parents=True, exist_ok=True)
            return root / 'results' / name
    return root / name


def plot_dir(signal, stage):
    """返回并创建指定信号、指定阶段的图片目录。"""
    path = signal_dir(signal) / "plots" / stage
    path.mkdir(parents=True, exist_ok=True)
    return path


ECG_AUDIT_CSV = {
    'ecg_preprocessing_report.csv', 'ecg_window_qc.csv', 'ecg_rr_bad_runs.csv',
    'ecg_qc_report.csv', 'ecg_qc_effective_summary.csv',
    'ecg_qc_effective_intervals.csv', 'final_ecg_processing_report.csv',
}

EDA_AUDIT_CSV = {'eda_preprocessing_report.csv', 'eda_window_qc.csv',
                 'eda_qc_report.csv', 'eda_qc_effective_summary.csv',
                 'eda_qc_effective_intervals.csv', 'final_eda_processing_report.csv'}


def csv_path(signal, filename, experiment_root=None):
    """Keep process reports in audit/, with legacy CSV read compatibility."""
    root = (Path(experiment_root) / signal if experiment_root is not None else signal_dir(signal))
    if (signal == 'ecg' and filename in ECG_AUDIT_CSV) or (signal == 'eda' and filename in EDA_AUDIT_CSV):
        target = root / 'audit' / filename
        legacy = root / filename
        if legacy.exists() and not target.exists():
            return legacy
        target.parent.mkdir(parents=True, exist_ok=True)
        return target
    return root / filename


def discover_files():
    """按配置查找全部源 FIF 文件。"""
    if not cfg.DATA_DIR.is_dir():
        raise FileNotFoundError("数据目录不存在：{}".format(cfg.DATA_DIR))
    files = sorted(cfg.DATA_DIR.glob(cfg.FILE_PATTERN))
    if not files:
        raise FileNotFoundError("没有找到文件：{} / {}".format(cfg.DATA_DIR, cfg.FILE_PATTERN))
    return files


def record_id(path):
    """从源文件名提取记录编号，如 A1_music1。"""
    suffix = "_physio_raw.fif"
    return path.name[:-len(suffix)] if path.name.endswith(suffix) else path.stem


def save_figure(fig, path):
    """按配置保存、显示并关闭 Matplotlib 图。"""
    if cfg.SAVE_PLOTS:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=cfg.PLOT_DPI, bbox_inches="tight", facecolor="white")
    if cfg.SHOW_PLOTS:
        plt.show()
    plt.close(fig)


def plot_indices(size, max_points=30000):
    """为长信号生成等距抽样切片，仅用于降低绘图负担。"""
    step = max(1, int(np.ceil(size / float(max_points))))
    return slice(None, None, step)


def write_csv(rows, path):
    """将字典列表写为 Excel 兼容的 UTF-8 CSV。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def save_config_snapshot():
    """保存本次运行的全部大写配置项。"""
    keys = [name for name in dir(cfg) if name.isupper()]
    values = {}
    for name in keys:
        value = getattr(cfg, name)
        if isinstance(value, Path):
            value = str(value)
        values[name] = value
    path = experiment_dir() / "config_snapshot.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(values, handle, ensure_ascii=False, indent=2, default=str)
    return path


def longest_true_interval(mask, window_starts, window_seconds):
    """返回连续True窗口覆盖的最长时间区间。"""
    best = None
    run_start = None
    for index, value in enumerate(mask):
        if value and run_start is None:
            run_start = index
        at_end = index == len(mask) - 1
        if run_start is not None and ((not value) or at_end):
            run_end = index if value and at_end else index - 1
            start_s = float(window_starts[run_start])
            end_s = float(window_starts[run_end] + window_seconds)
            candidate = (start_s, end_s)
            if best is None or candidate[1] - candidate[0] > best[1] - best[0]:
                best = candidate
            run_start = None
    return best


def window_starts(duration, window_seconds, step_seconds):
    """生成覆盖记录末尾的窗口起点。"""
    if duration <= window_seconds:
        return np.array([0.0])
    starts = list(np.arange(0.0, duration - window_seconds + 1e-9, step_seconds))
    final_start = duration - window_seconds
    if not starts or final_start - starts[-1] > 1e-6:
        starts.append(final_start)
    return np.asarray(starts, dtype=float)
