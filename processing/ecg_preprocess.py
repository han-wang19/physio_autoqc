"""ECG 初始预处理：读取、滤波、R 峰检测并保存中间结果。"""
from processing.common import output_path

from processing.common import csv_path

from contextlib import redirect_stdout
from io import StringIO
import json

import mne
import neurokit2 as nk
import numpy as np

import project_config as cfg
from processing.common import record_id, signal_dir, write_csv


def _detect_peaks(ecg, sfreq, method):
    """检测候选 R 峰，并隐藏 ProMAC 内部子检测器的已捕获报错。"""
    if method.lower() == "promac":
        # ProMAC 会自行跳过失败的成员检测器，却把错误列表打印到标准输出。
        # 这里仅隐藏这段内部日志；若 ProMAC 整体失败，异常仍会向外抛出。
        with redirect_stdout(StringIO()):
            _, info = nk.ecg_peaks(ecg, sampling_rate=sfreq, method=method)
    else:
        _, info = nk.ecg_peaks(ecg, sampling_rate=sfreq, method=method)
    return np.asarray(info["ECG_R_Peaks"], dtype=int)


def _refine_peaks(ecg, peaks, sfreq):
    """修正 RR 异常，按统一的双相形态/主偏转规则对齐峰位。"""
    _, fixed = nk.signal_fixpeaks(peaks, sampling_rate=sfreq, iterative=True)
    fixed = np.asarray(fixed, dtype=int)
    half_window = max(1, int(round(cfg.ECG_LOCAL_REFINEMENT_MS / 1000.0 * sfreq)))

    # 先从整条记录估计主 QRS 极性，再对所有心搏使用相同的 fiducial
    # point；逐搏选择绝对极值会在双相 QRS 中交替吸附到正、负峰。
    positive_excursions, negative_excursions, rs_candidates = [], [], []
    for peak in fixed:
        left = max(0, peak - half_window)
        right = min(ecg.size, peak + half_window + 1)
        local = ecg[left:right]
        if local.size:
            baseline = np.median(local)
            high, low = int(np.argmax(local)), int(np.argmin(local))
            positive, negative = float(local[high]-baseline), float(baseline-local[low])
            positive_excursions.append(positive)
            negative_excursions.append(negative)
            separation = (low-high) / sfreq * 1000.
            rs_candidates.append(
                0 < high < local.size-1 and 0 < low < local.size-1
                and negative > 0 and positive >= cfg.ECG_RS_MIN_POSITIVE_RATIO * negative
                and cfg.ECG_RS_SEPARATION_MS[0] <= separation <= cfg.ECG_RS_SEPARATION_MS[1])
    positive_qrs = (
        not negative_excursions
        or np.median(positive_excursions) >= np.median(negative_excursions)
    )
    # A consistent, substantial leading positive lobe is the R candidate of
    # an RS morphology even when its following negative lobe is deeper.
    # This is a shared heuristic, not a record-specific polarity exception.
    if rs_candidates and np.mean(rs_candidates) >= cfg.ECG_RS_MIN_CONSISTENCY:
        positive_qrs = True

    refined = []
    for peak in fixed:
        left = max(0, peak - half_window)
        right = min(ecg.size, peak + half_window + 1)
        local = ecg[left:right]
        if local.size:
            extremum = np.argmax(local) if positive_qrs else np.argmin(local)
            refined.append(left + int(extremum))
    return np.unique(np.asarray(refined, dtype=int)), fixed


def _load_filtered_ecg(path):
    """读取 ECG 通道并完成带通和工频陷波。"""
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.ECG_CHANNEL not in raw.ch_names:
        raise ValueError("缺少ECG通道：{}".format(cfg.ECG_CHANNEL))

    raw.pick([cfg.ECG_CHANNEL])
    from processing.units import calibrate_ecg_raw
    calibration = calibrate_ecg_raw(raw, path)
    sfreq = float(raw.info["sfreq"])
    duration = raw.n_times / sfreq
    raw.filter(*cfg.ECG_BANDPASS_HZ, picks=[cfg.ECG_CHANNEL], method="iir", verbose=False)
    raw.notch_filter(cfg.ECG_NOTCH_HZ, picks=[cfg.ECG_CHANNEL], method="iir", verbose=False)
    return raw.get_data()[0], sfreq, duration, calibration


def preprocess_one(path):
    """预处理一个 FIF 文件。"""
    rid = record_id(path)
    filtered, sfreq, duration, calibration = _load_filtered_ecg(path)
    initial_peaks = _detect_peaks(filtered, sfreq, cfg.ECG_PRIMARY_DETECTOR)
    refined_peaks, fixed_peaks = _refine_peaks(filtered, initial_peaks, sfreq)

    out_dir = output_path("ecg", "preprocessed")
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / "{}.npz".format(rid), filtered=filtered,
        primary_peaks=refined_peaks, initial_peaks=initial_peaks,
        sfreq=sfreq,
        peak_alignment_rule='consistent_rs_v2',
        unit='V', voltage_calibrated=True, calibration_json=json.dumps(calibration),
    )
    return {
        "record_id": rid, "source_file": str(path), "status": "OK",
        "sfreq_hz": sfreq, "duration_s": duration,
        "initial_peak_count": initial_peaks.size,
        "fixed_peak_count": fixed_peaks.size, "final_peak_count": refined_peaks.size,
        "ecg_sd": float(np.nanstd(filtered)), "ecg_range": float(np.nanmax(filtered) - np.nanmin(filtered)),
        "reason": "",
        **calibration,
    }


def preprocess_all_ecg(files):
    """预处理全部文件并写出汇总报告。"""
    rows = []
    for index, path in enumerate(files, 1):
        rid = record_id(path)
        print("  ECG {}/{}: {}".format(index, len(files), rid), flush=True)
        try:
            rows.append(preprocess_one(path))
        except Exception as exc:
            rows.append({"record_id": rid, "source_file": str(path), "status": "ERROR", "reason": str(exc)})
            print("    ERROR: {}".format(exc), flush=True)
    write_csv(rows, csv_path("ecg", "ecg_preprocessing_report.csv"))
    return rows
