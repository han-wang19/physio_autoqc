"""EDA 初始预处理：修复缺失、重采样、分解信号并保存 SCR。"""
from processing.common import output_path, csv_path

import mne
import neurokit2 as nk
import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import record_id, signal_dir, write_csv


def _interpolate_nonfinite(values):
    """线性插值 NaN/Inf，同时返回插值前的有限值掩码。"""
    finite = np.isfinite(values)
    if finite.all():
        return values.copy(), finite
    repaired = pd.Series(values).replace([np.inf, -np.inf], np.nan).interpolate(limit_direction="both").to_numpy()
    return repaired, finite


def _prepare_resampled_eda(raw):
    """换算为 µS、修复缺失、重采样，并映射原始缺失掩码。"""
    if cfg.EDA_SCALE_TO_US <= 0:
        raise ValueError("EDA_SCALE_TO_US必须为正数")
    source_values = raw.get_data()[0] * cfg.EDA_SCALE_TO_US
    if not np.isfinite(source_values).any():
        raise ValueError("EDA全部为非有限值")
    repaired_source, source_finite = _interpolate_nonfinite(source_values)
    source_time = np.arange(source_values.size) / float(raw.info["sfreq"])
    info = mne.create_info([cfg.EDA_CHANNEL], raw.info["sfreq"], ch_types=["misc"])
    working = mne.io.RawArray(repaired_source[np.newaxis, :], info, verbose=False)
    working.resample(cfg.EDA_TARGET_SFREQ, npad="auto", verbose=False)
    resampled = working.get_data()[0]
    target_time = np.arange(resampled.size) / cfg.EDA_TARGET_SFREQ
    finite_score = np.interp(target_time, source_time, source_finite.astype(float))
    target_finite = finite_score >= 0.999
    return resampled, target_finite, float(np.mean(source_finite))


def _filter_scr_info(info):
    """仅保留达到绝对幅度阈值的 SCR 及其配套字段。"""
    peaks = np.asarray(info.get("SCR_Peaks", []), dtype=int)
    amplitudes = np.asarray(info.get("SCR_Amplitude", []), dtype=float)
    n = min(peaks.size, amplitudes.size)
    keep = np.isfinite(amplitudes[:n]) & (amplitudes[:n] >= cfg.EDA_SCR_ABSOLUTE_MIN_US)
    filtered = {}
    for key in ("SCR_Onsets", "SCR_Peaks", "SCR_Recovery", "SCR_Amplitude", "SCR_RiseTime", "SCR_RecoveryTime"):
        values = np.asarray(info.get(key, []))
        filtered[key] = values[:n][keep] if values.size >= n else np.array([])
    return filtered


def _interpolate_masked(values, valid_mask):
    """为连续滤波/分解插值无效样本，同时由调用方保留原 validity mask。"""
    values = np.asarray(values, dtype=float)
    valid = np.asarray(valid_mask, dtype=bool) & np.isfinite(values)
    if np.sum(valid) < 2:
        raise ValueError("有效EDA样本不足，无法局部修复")
    samples = np.arange(values.size)
    return np.interp(samples, samples[valid], values[valid])


def _decompose_eda(values):
    """使用项目统一参数执行 clean、tonic/phasic 分解与 SCR 检测。"""
    clean = np.asarray(
        nk.eda_clean(values, sampling_rate=cfg.EDA_TARGET_SFREQ, method="neurokit")
    )
    components = nk.eda_phasic(
        clean, sampling_rate=cfg.EDA_TARGET_SFREQ, method="highpass"
    )
    tonic = components["EDA_Tonic"].to_numpy()
    phasic = components["EDA_Phasic"].to_numpy()
    from scipy.signal import find_peaks
    if find_peaks(phasic)[0].size == 0:
        # NeuroKit's relative threshold divides by max(prominence), which
        # raises for a flat/monotonic signal with no candidate peaks.
        empty = {key: np.array([]) for key in (
            'SCR_Onsets', 'SCR_Peaks', 'SCR_Recovery', 'SCR_Amplitude',
            'SCR_RiseTime', 'SCR_RecoveryTime')}
        return clean, tonic, phasic, empty
    _, info = nk.eda_peaks(
        phasic, sampling_rate=cfg.EDA_TARGET_SFREQ,
        method="neurokit", amplitude_min=cfg.EDA_SCR_DETECTOR_MIN,
    )
    return clean, tonic, phasic, _filter_scr_info(info)


def preprocess_one(path):
    """预处理一个 FIF 文件中的 EDA 通道。"""
    rid = record_id(path)

    # 读取信号，修复缺失并统一到目标采样率。
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.EDA_CHANNEL not in raw.ch_names:
        raise ValueError("缺少EDA通道：{}".format(cfg.EDA_CHANNEL))
    raw.pick([cfg.EDA_CHANNEL])
    source_sfreq = float(raw.info["sfreq"])
    raw_values, finite_mask, source_finite_fraction = _prepare_resampled_eda(raw)

    # NeuroKit 清理后分解 tonic/phasic，再从 phasic 中检测 SCR。
    clean, tonic, phasic, scr = _decompose_eda(raw_values)

    # 保存 QC 和最终处理所需的全部中间量。
    out_dir = output_path("eda", "preprocessed")
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / "{}.npz".format(rid), raw=raw_values, finite_mask=finite_mask,
        clean=clean, tonic=tonic, phasic=phasic,
        scr_peaks=scr["SCR_Peaks"].astype(int), scr_amplitudes=scr["SCR_Amplitude"],
        scr_rise=scr["SCR_RiseTime"], scr_recovery=scr["SCR_RecoveryTime"],
        sfreq=cfg.EDA_TARGET_SFREQ, unit=cfg.EDA_OUTPUT_UNIT, input_unit=cfg.EDA_INPUT_UNIT,
        scale_to_us=cfg.EDA_SCALE_TO_US,
        calibration_source=cfg.EDA_CALIBRATION_SOURCE,
    )
    return {
        "record_id": rid, "source_file": str(path), "status": "OK",
        "source_sfreq_hz": source_sfreq, "target_sfreq_hz": cfg.EDA_TARGET_SFREQ,
        "unit": cfg.EDA_OUTPUT_UNIT, "input_unit": cfg.EDA_INPUT_UNIT,
        "scale_to_us": cfg.EDA_SCALE_TO_US,
        "calibration_source": cfg.EDA_CALIBRATION_SOURCE,
        "duration_s": raw_values.size / cfg.EDA_TARGET_SFREQ,
        "source_finite_fraction": source_finite_fraction,
        "finite_fraction_after_mapping": float(np.mean(finite_mask)),
        "scr_count": int(len(scr["SCR_Peaks"])),
        "raw_range_us": float(np.nanmax(raw_values) - np.nanmin(raw_values)), "reason": "",
    }


def preprocess_all_eda(files):
    """依次预处理全部 EDA 文件并写出汇总报告。"""
    rows = []
    for index, path in enumerate(files, 1):
        rid = record_id(path)
        print("  EDA {}/{}: {}".format(index, len(files), rid), flush=True)
        try:
            rows.append(preprocess_one(path))
        except Exception as exc:
            rows.append({"record_id": rid, "source_file": str(path), "status": "ERROR", "reason": str(exc)})
            print("    ERROR: {}".format(exc), flush=True)
    write_csv(rows, csv_path("eda", "eda_preprocessing_report.csv"))
    return rows
