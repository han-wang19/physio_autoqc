"""EDA读取、重采样、清理、分解、SCR检测和阶段可视化。"""

import matplotlib.pyplot as plt
import mne
import neurokit2 as nk
import numpy as np
import pandas as pd

import project_config as cfg
from common import plot_dir, plot_indices, record_id, save_figure, signal_dir, write_csv


def _interpolate_nonfinite(values):
    finite = np.isfinite(values)
    if finite.all():
        return values.copy(), finite
    repaired = pd.Series(values).replace([np.inf, -np.inf], np.nan).interpolate(limit_direction="both").to_numpy()
    return repaired, finite


def _prepare_resampled_eda(raw):
    """先在源采样率修复非有限值，再重采样，并保留目标时间轴缺失掩码。"""
    source_values = raw.get_data()[0]
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
    peaks = np.asarray(info.get("SCR_Peaks", []), dtype=int)
    amplitudes = np.asarray(info.get("SCR_Amplitude", []), dtype=float)
    n = min(peaks.size, amplitudes.size)
    keep = np.isfinite(amplitudes[:n]) & (amplitudes[:n] >= cfg.EDA_SCR_ABSOLUTE_MIN)
    filtered = {}
    for key in ("SCR_Onsets", "SCR_Peaks", "SCR_Recovery", "SCR_Amplitude", "SCR_RiseTime", "SCR_RecoveryTime"):
        values = np.asarray(info.get(key, []))
        filtered[key] = values[:n][keep] if values.size >= n else np.array([])
    return filtered


def _plot_preprocess(raw, clean, tonic, phasic, info, sfreq, rid):
    """一张图展示EDA预处理所需的全部信息。"""
    time = np.arange(clean.size) / sfreq
    peaks = np.asarray(info["SCR_Peaks"], dtype=int)
    amplitudes = np.asarray(info["SCR_Amplitude"], dtype=float)
    view = plot_indices(clean.size)
    fig, axes = plt.subplots(3, 1, figsize=(18, 10), sharex=True)
    axes[0].plot(time[view], raw[view], color="#9ca3af", linewidth=0.5, label="Raw/resampled")
    axes[0].plot(time[view], clean[view], color="#2563eb", linewidth=0.8, label="Clean")
    axes[0].plot(time[view], tonic[view], color="#f97316", linewidth=1.0, label="Tonic")
    axes[0].set(title="{} | EDA preprocessing".format(rid), ylabel="Source units")
    axes[0].legend(ncol=3)
    axes[1].plot(time[view], phasic[view], color="#16a34a", linewidth=0.8)
    if peaks.size:
        axes[1].plot(time[peaks], phasic[peaks], "o", color="#dc2626", markersize=3)
        axes[2].stem(time[peaks], amplitudes, linefmt="#7c3aed", markerfmt="o", basefmt=" ")
    axes[1].set(ylabel="Phasic + SCR")
    axes[2].set(
        title="SCR >= {:.3g} source units".format(cfg.EDA_SCR_ABSOLUTE_MIN),
        xlabel="Time (s)", ylabel="SCR amplitude",
    )
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("eda", "01_preprocess") / "{}_eda_preprocess.png".format(rid))


def preprocess_one(path):
    rid = record_id(path)
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.EDA_CHANNEL not in raw.ch_names:
        raise ValueError("缺少EDA通道：{}".format(cfg.EDA_CHANNEL))
    raw.pick([cfg.EDA_CHANNEL])
    source_sfreq = float(raw.info["sfreq"])
    raw_values, finite_mask, source_finite_fraction = _prepare_resampled_eda(raw)
    clean = np.asarray(nk.eda_clean(raw_values, sampling_rate=cfg.EDA_TARGET_SFREQ, method="neurokit"))
    components = nk.eda_phasic(clean, sampling_rate=cfg.EDA_TARGET_SFREQ, method="highpass")
    tonic = components["EDA_Tonic"].to_numpy()
    phasic = components["EDA_Phasic"].to_numpy()
    _, info = nk.eda_peaks(
        phasic, sampling_rate=cfg.EDA_TARGET_SFREQ,
        method="neurokit", amplitude_min=cfg.EDA_SCR_DETECTOR_MIN,
    )
    scr = _filter_scr_info(info)

    out_dir = signal_dir("eda") / "preprocessed"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / "{}.npz".format(rid), raw=raw_values, finite_mask=finite_mask,
        clean=clean, tonic=tonic, phasic=phasic,
        scr_peaks=scr["SCR_Peaks"].astype(int), scr_amplitudes=scr["SCR_Amplitude"],
        scr_rise=scr["SCR_RiseTime"], scr_recovery=scr["SCR_RecoveryTime"],
        sfreq=cfg.EDA_TARGET_SFREQ,
    )
    _plot_preprocess(raw_values, clean, tonic, phasic, scr, cfg.EDA_TARGET_SFREQ, rid)
    return {
        "record_id": rid, "source_file": str(path), "status": "OK",
        "source_sfreq_hz": source_sfreq, "target_sfreq_hz": cfg.EDA_TARGET_SFREQ,
        "duration_s": raw_values.size / cfg.EDA_TARGET_SFREQ,
        "source_finite_fraction": source_finite_fraction,
        "finite_fraction_after_mapping": float(np.mean(finite_mask)),
        "scr_count": int(len(scr["SCR_Peaks"])),
        "raw_range": float(np.nanmax(raw_values) - np.nanmin(raw_values)), "reason": "",
    }


def preprocess_all_eda(files):
    rows = []
    for index, path in enumerate(files, 1):
        rid = record_id(path)
        print("  EDA {}/{}: {}".format(index, len(files), rid), flush=True)
        try:
            rows.append(preprocess_one(path))
        except Exception as exc:
            rows.append({"record_id": rid, "source_file": str(path), "status": "ERROR", "reason": str(exc)})
            print("    ERROR: {}".format(exc), flush=True)
    write_csv(rows, signal_dir("eda") / "eda_preprocessing_report.csv")
    return rows
