"""Generate final subject-by-music QC coverage figures for every concert/modality."""

from processing.common import csv_path

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd

import project_config as cfg


STATUS_VALUE = {"EXCLUDED": 0, "TRIMMED": 1, "PASS": 2}
STATUS_LABEL = {0: "Excluded", 1: "Trimmed", 2: "Pass"}
CELL_LABEL = {0: "X", 1: "T", 2: "P"}
COLORS = ["#D3D3D3", "#F4A261", "#2A9D8F"]
MANUAL_STATUS = {
    "INCLUDE_FULL": "PASS",
    "INCLUDE_TRIMMED": "TRIMMED",
    "EXCLUDE": "EXCLUDED",
    "SECOND_REVIEW": "EXCLUDED",
}


def _natural_key(value):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", str(value))]


from reports.record_labels import record_parts as _record_parts, segment_label, segment_key


def _final_status(auto_status):
    status = str(auto_status).upper()
    if status in ("PASS", "MANUAL_PASS"):
        return "PASS"
    if status in ("TRIMMED", "MANUAL_TRIMMED"):
        return "TRIMMED"
    return "EXCLUDED"


def load_final_dispositions(experiment_dir, signal):
    """Merge automatic report and optional manual decisions into final dispositions."""
    report_path = csv_path(signal, "{}_qc_report.csv".format(signal), experiment_dir)
    report = pd.read_csv(report_path, encoding="utf-8-sig")
    status = {
        str(row.record_id): _final_status(row.auto_status)
        for row in report.itertuples(index=False)
    }

    expected = {str(row.record_id): str(row.auto_status) for row in report.itertuples(index=False)}
    reviewed_at = {}
    manual_path = experiment_dir / "manual_qc_decisions.csv"
    if manual_path.exists():
        manual = pd.read_csv(manual_path, encoding="utf-8-sig").fillna("")
        manual = manual[manual["signal"].str.lower() == signal.lower()]
        for row in manual.itertuples(index=False):
            status[str(row.record_id)] = MANUAL_STATUS.get(
                str(row.manual_status).upper(), "EXCLUDED"
            )
            expected[str(row.record_id)] = {
                'INCLUDE_FULL': 'MANUAL_PASS', 'INCLUDE_TRIMMED': 'MANUAL_TRIMMED',
                'EXCLUDE': 'MANUAL_EXCLUDE', 'SECOND_REVIEW': 'SECOND_REVIEW',
            }.get(str(row.manual_status), 'SECOND_REVIEW')
            reviewed_at[str(row.record_id)] = str(row.reviewed_at)
    final_path = csv_path(signal, 'final_{}_processing_report.csv'.format(signal), experiment_dir)
    if final_path.exists():
        final = pd.read_csv(final_path, encoding='utf-8-sig').fillna('')
        for row in final.itertuples(index=False):
            rid = str(row.record_id)
            if str(row.qc_status) != expected.get(rid):
                continue  # An older automatic run cannot override a new manual decision.
            if rid in reviewed_at and str(getattr(row, 'manual_reviewed_at', '')) != reviewed_at[rid]:
                continue
            if str(row.processing_status) != 'INCLUDED':
                status[rid] = 'EXCLUDED'
    return status


def save_summary_figure(dispositions, concert, signal, output_path):
    parsed = [(_record_parts(record_id), status)
              for record_id, status in dispositions.items()]
    subjects = sorted({parts[0] for parts, _ in parsed}, key=_natural_key)
    music = sorted({parts[1] for parts, _ in parsed}, key=segment_key)
    subject_index = {label: index for index, label in enumerate(subjects)}
    music_index = {label: index for index, label in enumerate(music)}

    # Only combinations represented in the QC report enter the denominator.
    matrix = np.full((len(subjects), len(music)), np.nan)
    for (subject, track), status in parsed:
        matrix[subject_index[subject], music_index[track]] = STATUS_VALUE[status]

    cmap = ListedColormap(COLORS)
    cmap.set_bad("white")
    norm = BoundaryNorm(np.arange(-0.5, 3.5, 1.0), cmap.N)
    height = max(7.5, 4.8 + 0.52 * len(subjects))
    fig = plt.figure(figsize=(15, height), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.8])

    ax = fig.add_subplot(grid[0])
    image = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(music)), [segment_label(label) for label in music])
    ax.set_yticks(range(len(subjects)), ["A{}".format(label) for label in subjects])
    ax.set(xlabel="Segment / condition", ylabel="Subject")
    modality = "ECG/HRV" if signal == "ecg" else signal.upper()
    ax.set_title(
        "{} | {} Final Inclusion After QC and Processing".format(concert, modality),
        fontsize=16, weight="bold", pad=14,
    )
    ax.set_xticks(np.arange(-0.5, len(music), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(subjects), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            if np.isnan(matrix[row, column]):
                continue
            value = int(matrix[row, column])
            color = "white" if value == STATUS_VALUE["PASS"] else "#333333"
            ax.text(column, row, CELL_LABEL[value], ha="center", va="center",
                    color=color, weight="bold", fontsize=11)

    colorbar = fig.colorbar(image, ax=ax, ticks=range(3), shrink=0.85, pad=0.02)
    colorbar.ax.set_yticklabels([STATUS_LABEL[value] for value in range(3)])
    colorbar.set_label("QC disposition")

    ax_count = fig.add_subplot(grid[1])
    x = np.arange(len(subjects))
    bottom = np.zeros(len(subjects), dtype=int)
    for value, label in ((2, "Pass"), (1, "Trimmed"), (0, "Excluded")):
        counts = np.sum(matrix == value, axis=1)
        ax_count.bar(x, counts, bottom=bottom, color=COLORS[value], label=label)
        bottom += counts

    recorded = np.sum(np.isfinite(matrix), axis=1)
    included = np.sum(matrix >= STATUS_VALUE["TRIMMED"], axis=1)
    label_y = max(recorded.max(initial=0), 1) + 0.18
    for index, count in enumerate(included):
        ax_count.text(index, label_y,
                      "Included {}/{}".format(int(count), int(recorded[index])),
                      ha="center", va="bottom", fontsize=10, weight="bold")
    ax_count.set_xticks(x, ["A{}".format(subject) for subject in subjects])
    ax_count.set_ylim(0, label_y + 1.2)
    ax_count.set(
        ylabel="Number of segments", xlabel="Subject",
        title="Coverage Count by Subject",
    )
    ax_count.grid(axis="y", alpha=0.25)
    ax_count.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)

    total_included = int(np.sum(matrix >= STATUS_VALUE["TRIMMED"]))
    total_recorded = int(np.sum(np.isfinite(matrix)))
    percent = 100.0 * total_included / total_recorded if total_recorded else 0.0
    fig.suptitle(
        "Included: {}/{} segments ({:.1f}%)  |  P = pass, T = trimmed, X = excluded".format(
            total_included, total_recorded, percent
        ),
        y=1.025, fontsize=11,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=cfg.PLOT_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate_all(output_root=None):
    """Generate one final summary per available concert and signal modality."""
    output_root = Path(output_root or cfg.OUTPUT_ROOT)
    generated = []
    for experiment_dir in sorted(
            (path for path in output_root.iterdir() if path.is_dir()),
            key=lambda path: _natural_key(path.name)):
        for signal in ("ecg", "eda"):
            report = csv_path(signal, "{}_qc_report.csv".format(signal), experiment_dir)
            if not report.exists():
                continue
            dispositions = load_final_dispositions(experiment_dir, signal)
            path = experiment_dir / signal / "plots" / "qc_summary" / "subject_music_qc_coverage.png"
            save_summary_figure(dispositions, experiment_dir.name, signal, path)
            generated.append(path)
            print("Generated: {}".format(path))
    return generated


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate concert-by-modality QC summary figures"
    )
    parser.add_argument("--output-root", type=Path, default=cfg.OUTPUT_ROOT)
    args = parser.parse_args()
    generate_all(args.output_root)
