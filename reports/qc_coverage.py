"""绘制被试 × 曲目的记录级 QC 覆盖汇总图。"""

import re

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

from processing.common import plot_dir, save_figure


STATUS_VALUE = {"MISSING": 0, "FAIL": 1, "REVIEW": 2, "TRIMMED": 3, "PASS": 4}
STATUS_LABEL = {0: "Missing", 1: "Fail", 2: "Review", 3: "Trimmed", 4: "Pass"}
CELL_LABEL = {0: "–", 1: "X", 2: "R", 3: "T", 4: "P"}
COLORS = ["#F2F2F2", "#D3D3D3", "#E9C46A", "#F4A261", "#2A9D8F"]


from reports.record_labels import record_parts as _record_parts, segment_label, segment_key


def _natural_key(value):
    """让字符串中的数字按数值顺序排列。"""
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def save_qc_coverage_plot(reports, signal, stage, title, filename="all_subjects_music_qc_coverage.png"):
    """保存本场实验全部被试 × 曲目的 QC 矩阵和覆盖计数。"""
    parsed = [(_record_parts(row["record_id"]), row["auto_status"]) for row in reports]
    subjects = sorted({parts[0] for parts, _ in parsed}, key=_natural_key)
    music = sorted({parts[1] for parts, _ in parsed}, key=segment_key)

    matrix = np.zeros((len(subjects), len(music)), dtype=int)
    subject_index = {subject: index for index, subject in enumerate(subjects)}
    music_index = {label: index for index, label in enumerate(music)}
    for (subject, label), status in parsed:
        matrix[subject_index[subject], music_index[label]] = STATUS_VALUE[status]

    cmap = ListedColormap(COLORS)
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1.0), cmap.N)
    fig = plt.figure(figsize=(15, 7.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.8])

    ax = fig.add_subplot(grid[0])
    image = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(music)), [segment_label(label) for label in music])
    ax.set_yticks(range(len(subjects)), ["A{}".format(subject) for subject in subjects])
    ax.set(xlabel="Segment / condition", ylabel="Subject")
    ax.set_title(title, fontsize=16, weight="bold", pad=14)
    ax.set_xticks(np.arange(-0.5, len(music), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(subjects), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            text_color = "white" if value == 4 else "#333333"
            ax.text(column, row, CELL_LABEL[value], ha="center", va="center",
                    color=text_color, weight="bold", fontsize=11)

    colorbar = fig.colorbar(image, ax=ax, ticks=range(5), shrink=0.85, pad=0.02)
    colorbar.ax.set_yticklabels([STATUS_LABEL[value] for value in range(5)])
    colorbar.set_label("QC disposition")

    ax_count = fig.add_subplot(grid[1])
    x = np.arange(len(subjects))
    bottom = np.zeros(len(subjects), dtype=int)
    stack_order = [(4, "Pass"), (3, "Trimmed"), (2, "Review"), (1, "Fail"), (0, "Missing")]
    for value, label in stack_order:
        counts = np.sum(matrix == value, axis=1)
        ax_count.bar(x, counts, bottom=bottom, color=COLORS[value], label=label)
        bottom += counts

    included = np.sum(matrix >= STATUS_VALUE["TRIMMED"], axis=1)
    for index, count in enumerate(included):
        ax_count.text(index, len(music) + 0.2, "Included {}/{}".format(count, len(music)),
                      ha="center", va="bottom", fontsize=10, weight="bold")
    ax_count.set_xticks(x, ["A{}".format(subject) for subject in subjects])
    ax_count.set_ylim(0, len(music) + 1.2)
    ax_count.set(
        ylabel="Number of segments", xlabel="Subject",
        title="Coverage Count by Subject",
    )
    ax_count.grid(axis="y", alpha=0.25)
    ax_count.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)

    total_included = int(np.sum(matrix >= STATUS_VALUE["TRIMMED"]))
    total_recorded = int(np.sum(matrix > STATUS_VALUE["MISSING"]))
    fig.suptitle(
        "Included: {}/{} recorded segments | P=pass, T=trimmed, R=review, X=fail, –=missing".format(
            total_included, total_recorded
        ),
        y=1.035, fontsize=11,
    )
    save_figure(
        fig,
        plot_dir(signal, stage) / filename,
    )
