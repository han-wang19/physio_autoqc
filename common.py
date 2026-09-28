"""ECG和EDA共用的文件、输出和绘图辅助函数。"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import project_config as cfg


def experiment_dir():
    return cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME


def signal_dir(signal):
    path = experiment_dir() / signal.lower()
    path.mkdir(parents=True, exist_ok=True)
    return path


def plot_dir(signal, stage):
    path = signal_dir(signal) / "plots" / stage
    path.mkdir(parents=True, exist_ok=True)
    return path


def discover_files():
    if not cfg.DATA_DIR.is_dir():
        raise FileNotFoundError("数据目录不存在：{}".format(cfg.DATA_DIR))
    files = sorted(cfg.DATA_DIR.glob(cfg.FILE_PATTERN))
    if not files:
        raise FileNotFoundError("没有找到文件：{} / {}".format(cfg.DATA_DIR, cfg.FILE_PATTERN))
    return files


def record_id(path):
    suffix = "_physio_raw.fif"
    return path.name[:-len(suffix)] if path.name.endswith(suffix) else path.stem


def save_figure(fig, path):
    if cfg.SAVE_PLOTS:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=cfg.PLOT_DPI, bbox_inches="tight", facecolor="white")
    if cfg.SHOW_PLOTS:
        plt.show()
    plt.close(fig)


def plot_indices(size, max_points=30000):
    step = max(1, int(np.ceil(size / float(max_points))))
    return slice(None, None, step)


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def save_config_snapshot():
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
        json.dump(values, handle, ensure_ascii=False, indent=2)
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
