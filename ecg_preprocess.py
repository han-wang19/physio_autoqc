"""ECG读取、滤波、R峰检测和阶段可视化。"""

import matplotlib.pyplot as plt
import mne
import neurokit2 as nk
import numpy as np

import project_config as cfg
from common import plot_dir, plot_indices, record_id, save_figure, signal_dir, write_csv


def _detect_peaks(ecg, sfreq, method):
    _, info = nk.ecg_peaks(ecg, sampling_rate=sfreq, method=method)
    return np.asarray(info["ECG_R_Peaks"], dtype=int)


def _refine_peaks(ecg, peaks, sfreq):
    """先修正异常峰，再在局部搜索绝对幅度最大的极值，兼容反相ECG。"""
    _, fixed = nk.signal_fixpeaks(peaks, sampling_rate=sfreq, iterative=True)
    fixed = np.asarray(fixed, dtype=int)
    half_window = max(1, int(round(cfg.ECG_LOCAL_REFINEMENT_MS / 1000.0 * sfreq)))
    refined = []
    for peak in fixed:
        left = max(0, peak - half_window)
        right = min(ecg.size, peak + half_window + 1)
        local = ecg[left:right]
        if local.size:
            refined.append(left + int(np.argmax(np.abs(local - np.median(local)))))
    return np.unique(np.asarray(refined, dtype=int)), fixed


def _rpeak_plot(ecg, peaks, sfreq, rid):
    time = np.arange(ecg.size) / sfreq
    view = plot_indices(ecg.size)
    heart_rate = nk.ecg_rate(peaks, sampling_rate=sfreq) if peaks.size > 1 else np.array([])
    peak_time = peaks / sfreq
    hr_time = peak_time[1:] if heart_rate.size == peaks.size - 1 else peak_time[:heart_rate.size]
    fig, axes = plt.subplots(2, 1, figsize=(18, 8), sharex=True)
    axes[0].plot(time[view], ecg[view], linewidth=0.6, color="#2563eb", label="Filtered ECG")
    axes[0].plot(peak_time, ecg[peaks], "o", color="#dc2626", markersize=2, label="R peaks")
    axes[0].set(title="{} | Refined R peaks".format(rid), ylabel="Amplitude")
    axes[0].legend()
    if heart_rate.size:
        axes[1].plot(hr_time, heart_rate, color="#16a34a", linewidth=1.0)
    axes[1].set(xlabel="Time (s)", ylabel="Heart rate (bpm)")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("ecg", "01_preprocess") / "{}_ecg_peaks_hr.png".format(rid))


def preprocess_one(path):
    rid = record_id(path)
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.ECG_CHANNEL not in raw.ch_names:
        raise ValueError("缺少ECG通道：{}".format(cfg.ECG_CHANNEL))
    raw.pick([cfg.ECG_CHANNEL])
    sfreq = float(raw.info["sfreq"])
    raw_ecg = raw.get_data()[0]
    filtered_raw = raw.copy().filter(*cfg.ECG_BANDPASS_HZ, method="iir", verbose=False)
    filtered_raw.notch_filter(cfg.ECG_NOTCH_HZ, method="iir", verbose=False)
    filtered = filtered_raw.get_data()[0]

    primary_initial = _detect_peaks(filtered, sfreq, cfg.ECG_PRIMARY_DETECTOR)
    primary, fixed = _refine_peaks(filtered, primary_initial, sfreq)
    out_dir = signal_dir("ecg") / "preprocessed"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / "{}.npz".format(rid), filtered=filtered,
        primary_peaks=primary, initial_peaks=primary_initial,
        fixed_peaks=fixed, sfreq=sfreq,
    )
    _rpeak_plot(filtered, primary, sfreq, rid)
    return {
        "record_id": rid, "source_file": str(path), "status": "OK",
        "sfreq_hz": sfreq, "duration_s": raw_ecg.size / sfreq,
        "initial_peak_count": primary_initial.size,
        "fixed_peak_count": fixed.size, "final_peak_count": primary.size,
        "ecg_sd": float(np.nanstd(filtered)), "ecg_range": float(np.nanmax(filtered) - np.nanmin(filtered)),
        "reason": "",
    }


def preprocess_all_ecg(files):
    rows = []
    for index, path in enumerate(files, 1):
        rid = record_id(path)
        print("  ECG {}/{}: {}".format(index, len(files), rid), flush=True)
        try:
            rows.append(preprocess_one(path))
        except Exception as exc:
            rows.append({"record_id": rid, "source_file": str(path), "status": "ERROR", "reason": str(exc)})
            print("    ERROR: {}".format(exc), flush=True)
    write_csv(rows, signal_dir("ecg") / "ecg_preprocessing_report.csv")
    return rows
