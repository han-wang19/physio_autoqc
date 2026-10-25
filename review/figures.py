"""Build review figures without starting a UI or saving research outputs."""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import project_config as cfg
from review.manual_qc import effective_manual_intervals


def _axis_range(values, scale):
    """以整条记录的纵轴范围为 1 倍，返回独立纵轴缩放范围。"""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return [-1.0, 1.0]
    low, high = float(np.min(values)), float(np.max(values))
    span = max(high - low, np.finfo(float).eps)
    center = (low + high) / 2.0
    half = 0.55 * span * float(scale)
    return [center - half, center + half]


def _subtract_intervals(start, end, exclusions):
    """Subtract non-overlapping display intervals from one interval."""
    pieces = [(float(start), float(end))]
    for exclusion in exclusions.itertuples(index=False):
        cut_start, cut_end = float(exclusion.start_s), float(exclusion.end_s)
        updated = []
        for left, right in pieces:
            if cut_end <= left or cut_start >= right:
                updated.append((left, right))
            else:
                if left < cut_start:
                    updated.append((left, cut_start))
                if cut_end < right:
                    updated.append((cut_end, right))
        pieces = updated
    return pieces


def _add_manual_intervals(fig, intervals, rows):
    """Draw only the left/right boundaries of each manual annotation."""
    colors = {"PASS": "#16a34a", "REVIEW": "#f59e0b", "FAIL": "#dc2626"}
    intervals = effective_manual_intervals(intervals)
    shapes = list(fig.layout.shapes or ())
    for item in intervals.itertuples(index=False):
        color = colors.get(str(item.label), "#64748b")
        for row_index in range(1, rows + 1):
            for boundary in (float(item.start_s), float(item.end_s)):
                suffix = str(row_index) if row_index > 1 else ""
                shapes.append(dict(
                    type="line", x0=boundary, x1=boundary, y0=0, y1=1,
                    xref="x" + suffix, yref="y" + suffix + " domain",
                    line={"color": color, "width": 2, "dash": "dash"},
                    layer="above",
                ))
    fig.update_layout(shapes=shapes)


def _add_effective_qc_spans(fig, bad_windows, auto_pass_intervals, manual, rows):
    """Paint each time interval once, using the final manual override where present."""
    manual = effective_manual_intervals(manual)
    colors = {"PASS": "#22c55e", "REVIEW": "#f59e0b", "FAIL": "#ef4444"}
    spans = []
    for item in bad_windows.itertuples(index=False):
        for left, right in _subtract_intervals(item.start_s, item.end_s, auto_pass_intervals):
            for start, end in _subtract_intervals(left, right, manual):
                spans.append((start, end, item.status))
    for item in auto_pass_intervals.itertuples(index=False):
        for start, end in _subtract_intervals(item.start_s, item.end_s, manual):
            spans.append((start, end, "PASS"))
    spans.extend((item.start_s, item.end_s, item.label) for item in manual.itertuples(index=False))
    # Different failure reasons share one background color. Merge them before
    # creating shapes; adding vrects individually repeatedly copies the layout.
    merged = []
    for start, end, status in sorted(spans):
        if end <= start:
            continue
        if merged and merged[-1][2] == status and start <= merged[-1][1] + 1e-9:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]), status)
        else:
            merged.append((start, end, status))
    shapes = list(fig.layout.shapes or ())
    for start, end, status in merged:
        for row in range(1, rows + 1):
            suffix = str(row) if row > 1 else ""
            shapes.append(dict(type="rect", x0=start, x1=end, y0=0, y1=1,
                               xref="x" + suffix, yref="y" + suffix + " domain",
                               fillcolor=colors[status], opacity=0.09,
                               line_width=0, layer="below"))
    fig.update_layout(shapes=shapes)
    for status in ("PASS", "REVIEW", "FAIL"):
        if any(item[2] == status for item in spans):
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=status,
                                    marker={"color": colors[status], "symbol": "square"},
                                    hoverinfo="skip"), row=1, col=1)


def plot_indices(signal, peaks, max_points=250_000):
    """长记录使用 min/max 包络降采样，同时保留 R 峰附近的形态。"""
    size = signal.size
    if size <= max_points:
        return np.arange(size)
    stride = max(2, int(np.ceil(size / (max_points / 2.0))))
    complete = size // stride * stride
    blocks = signal[:complete].reshape(-1, stride)
    offsets = np.arange(blocks.shape[0]) * stride
    envelope = np.concatenate((offsets + np.argmin(blocks, axis=1), offsets + np.argmax(blocks, axis=1)))
    tail = np.arange(complete, size)
    peak_neighborhood = (peaks[:, None] + np.arange(-2, 3)).ravel() if peaks.size else np.array([], dtype=int)
    peak_neighborhood = peak_neighborhood[(peak_neighborhood >= 0) & (peak_neighborhood < size)]
    return np.unique(np.concatenate((envelope, tail, peak_neighborhood)))


def rr_diagnostics(peaks, sfreq):
    """返回瞬时心率及与自动 QC 规则一致的异常 RR 掩码。"""
    peak_times = peaks / sfreq
    rr_seconds = np.diff(peak_times)
    heart_rate = np.divide(
        60.0, rr_seconds,
        out=np.full(rr_seconds.shape, np.nan, dtype=float),
        where=rr_seconds > 0,
    )
    local_medians = np.empty(rr_seconds.size, dtype=float)
    for index in range(rr_seconds.size):
        left, right = max(0, index - 2), min(rr_seconds.size, index + 3)
        local_medians[index] = np.median(rr_seconds[left:right])
    deviation = np.abs(rr_seconds - local_medians) / np.maximum(local_medians, np.finfo(float).eps)
    low, high = cfg.ECG_PLAUSIBLE_HR_BPM
    abnormal = (~np.isfinite(heart_rate)) | (heart_rate < low) | (heart_rate > high)
    abnormal |= deviation > cfg.ECG_RR_LOCAL_DEVIATION
    return peak_times, heart_rate, abnormal


def ecg_figure(data, peaks, initial_start, initial_end, record_id, bad_windows,
               auto_pass_intervals, manual_intervals,
               waveform_y_scale, heart_rate_y_scale, drag_mode):
    ecg, sfreq = data["filtered"], float(data["sfreq"])
    # The browser draws a pixel envelope from the full waveform for the current
    # viewport. Keep the initial overview small to avoid serializing it twice.
    samples = plot_indices(ecg, peaks, max_points=30_000)
    peak_times, heart_rate, abnormal = rr_diagnostics(peaks, sfreq)

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        row_heights=[0.68, 0.32],
    )
    fig.add_trace(go.Scattergl(
        x=samples / sfreq, y=ecg[samples], mode="lines", name="Filtered ECG",
        line={"width": 1},
    ), row=1, col=1)
    fig.add_trace(go.Scattergl(
        x=peak_times, y=ecg[peaks], mode="markers", name="R peaks",
        marker={"color": "red", "size": 7},
    ), row=1, col=1)
    fig.add_trace(go.Scattergl(
        x=peak_times[1:], y=heart_rate, mode="lines+markers", name="Heart rate",
        line={"color": "#16a34a", "width": 1.5},
        marker={"size": 4}, connectgaps=False,
        hovertemplate="%{x:.3f} s<br>%{y:.1f} bpm<extra></extra>",
    ), row=2, col=1)
    fig.add_trace(go.Scattergl(
        x=peak_times[1:][abnormal], y=heart_rate[abnormal], mode="markers",
        name="Abnormal RR", marker={"color": "#dc2626", "size": 9, "symbol": "x"},
        hovertemplate="异常 RR<br>%{x:.3f} s<br>%{y:.1f} bpm<extra></extra>",
    ), row=2, col=1)
    low, high = cfg.ECG_PLAUSIBLE_HR_BPM
    for limit in (low, high):
        fig.add_hline(y=limit, line_color="#94a3b8", line_dash="dot", row=2, col=1)
    _add_effective_qc_spans(fig, bad_windows, auto_pass_intervals, manual_intervals, rows=2)
    _add_manual_intervals(fig, manual_intervals, rows=2)
    fig.update_yaxes(
        title_text="Filtered ECG", fixedrange=True,
        range=_axis_range(ecg, waveform_y_scale),
        uirevision="ecg-wave-y-{}-{}".format(record_id, waveform_y_scale), row=1, col=1,
    )
    finite_hr = heart_rate[np.isfinite(heart_rate)]
    fig.update_yaxes(
        title_text="Heart rate (bpm)", fixedrange=True,
        range=_axis_range(
            finite_hr if finite_hr.size else np.asarray(cfg.ECG_PLAUSIBLE_HR_BPM),
            heart_rate_y_scale,
        ),
        uirevision="ecg-hr-y-{}-{}".format(record_id, heart_rate_y_scale), row=2, col=1,
    )
    fig.update_xaxes(title_text="Time (s)", range=[initial_start, initial_end], row=2, col=1)
    fig.update_xaxes(fixedrange=False, uirevision="ecg-time-{}".format(record_id))
    fig.update_layout(
        height=650,
        dragmode=drag_mode,
        clickmode="event+select",
        selectdirection="h",
        hovermode="x unified",
        margin={"t": 20, "b": 45},
        uirevision="ecg-{}".format(record_id),
    )
    return fig, int(np.count_nonzero(abnormal)), int(heart_rate.size), int(samples.size)


def eda_figure(data, initial_start, initial_end, record_id, bad_windows,
               auto_pass_intervals, manual_intervals,
               level_y_scale, phasic_y_scale, drag_mode):
    # Load the heavier processing dependencies only for EDA review.
    from processing.eda_preprocess import _decompose_eda, _interpolate_masked
    from processing.eda_finalize import _qc_valid_mask

    sfreq = float(data["sfreq"])
    valid_mask = _qc_valid_mask(
        record_id, 0.0, data["raw"].size, sfreq, data["finite_mask"], manual_intervals
    )
    if np.sum(valid_mask) >= 2:
        repaired = _interpolate_masked(data["raw"], valid_mask)
        clean, tonic, phasic, scr = _decompose_eda(repaired)
        peaks = np.asarray(scr["SCR_Peaks"], dtype=int)
        peaks = peaks[(peaks >= 0) & (peaks < valid_mask.size)]
        peaks = peaks[valid_mask[peaks]]
    else:
        clean, tonic, phasic = data["clean"], data["tonic"], data["phasic"]
        peaks = data["scr_peaks"].astype(int)
    samples = plot_indices(data["clean"], np.array([], dtype=int))
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
        row_heights=[0.62, 0.38],
    )
    for key, label, color in (
        ("raw", "Raw/resampled", "#9ca3af"),
    ):
        fig.add_trace(go.Scattergl(
            x=samples / sfreq, y=data[key][samples], mode="lines", name=label,
            line={"color": color, "width": 1},
        ), row=1, col=1)
    fig.add_trace(go.Scattergl(
        x=samples / sfreq, y=clean[samples], mode="lines",
        name="Repaired clean (final pipeline)",
        line={"color": "#2563eb", "width": 1},
    ), row=1, col=1)
    fig.add_trace(go.Scattergl(
        x=samples / sfreq, y=tonic[samples], mode="lines", name="Tonic",
        line={"color": "#f97316", "width": 1},
    ), row=1, col=1)
    fig.add_trace(go.Scattergl(
        x=samples / sfreq, y=phasic[samples], mode="lines", name="Phasic",
        line={"color": "#16a34a", "width": 1},
    ), row=2, col=1)
    if peaks.size:
        fig.add_trace(go.Scattergl(
            x=peaks / sfreq, y=phasic[peaks], mode="markers", name="SCR peaks",
            marker={"color": "#dc2626", "size": 7},
        ), row=2, col=1)
    _add_effective_qc_spans(fig, bad_windows, auto_pass_intervals, manual_intervals, rows=2)
    _add_manual_intervals(fig, manual_intervals, rows=2)
    level_values = np.concatenate((data["raw"], clean, tonic))
    fig.update_yaxes(
        title_text="EDA level (µS)", fixedrange=True,
        range=_axis_range(level_values, level_y_scale),
        uirevision="eda-level-y-{}-{}".format(record_id, level_y_scale), row=1, col=1,
    )
    fig.update_yaxes(
        title_text="Phasic (µS)", fixedrange=True,
        range=_axis_range(phasic, phasic_y_scale),
        uirevision="eda-phasic-y-{}-{}".format(record_id, phasic_y_scale), row=2, col=1,
    )
    fig.update_xaxes(fixedrange=False, uirevision="eda-time-{}".format(record_id))
    fig.update_xaxes(title_text="Time (s)", range=[initial_start, initial_end], row=2, col=1)
    fig.update_layout(
        height=650, dragmode=drag_mode, hovermode="x unified",
        selectdirection="h",
        margin={"t": 20, "b": 45}, uirevision="eda-{}".format(record_id),
    )
    return fig
