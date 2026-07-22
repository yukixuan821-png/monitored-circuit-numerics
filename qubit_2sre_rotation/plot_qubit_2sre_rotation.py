"""Plot the committed Figure 4 data table."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DATA_FILE = Path(__file__).with_name("qubit_2sre_rotation_data.txt")
OUTPUT_PNG = Path(__file__).with_name("qubit_2sre_rotation.png")
OUTPUT_PDF = Path(__file__).with_name("qubit_2sre_rotation.pdf")
COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]
MARKERS = ["o", "s", "^", "D"]


def read_results(path=DATA_FILE):
    rows = []
    with open(path, encoding="utf-8") as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        for row in reader:
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def plot_results(path=DATA_FILE, output_png=OUTPUT_PNG, output_pdf=OUTPUT_PDF):
    plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix", "font.size": 13, "axes.labelsize": 13, "legend.fontsize": 12, "xtick.labelsize": 13, "ytick.labelsize": 13, "pdf.fonttype": 42, "svg.fonttype": "none"})
    rows = read_results(path)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.17, top=0.96)
    all_means = []
    for color, marker, n in zip(COLORS, MARKERS, (6, 8, 10, 12)):
        group = sorted(
            (row for row in rows if int(row["N"]) == n),
            key=lambda row: row["theta_M"],
        )
        theta = np.array([row["theta_M"] for row in group])
        mean = np.array([row["S2_mean"] for row in group])
        sem = np.array([row["S2_sem"] for row in group])
        if np.any(mean <= 0):
            raise ValueError("log plot requires strictly positive aggregate means")
        lower_endpoint = np.maximum.reduce((mean-sem, mean*1e-3, np.full_like(mean, np.finfo(float).tiny)))
        asymmetric_yerr = np.vstack((mean-lower_endpoint, sem))
        ax.errorbar(theta, mean, yerr=asymmetric_yerr, color=color, marker=marker, linewidth=2, markersize=6, capsize=2, label=rf"$N={n}$")
        all_means.extend(mean)
    theta_grid = np.unique([row["theta_M"] for row in rows])
    anchor_theta = theta_grid[len(theta_grid)//2]
    anchor_y = float(np.median(all_means))
    ax.plot(theta_grid, anchor_y*(theta_grid/anchor_theta)**2, color="#52514e", linestyle="--", linewidth=1.5, label=r"slope $2$")
    ax.set(xscale="log", yscale="log", xlabel=r"$\theta_M$", ylabel=r"$\overline{S}_2(\theta_M)$", title="Steady-state 2-SRE vs rotation angle")
    ax.legend(frameon=False)
    ax.grid(True, which="major", color="#e1e0d9", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    for output, extra_metadata in ((output_png, {"Description": f"plot-only from {Path(path).name}"}), (output_pdf, {"Subject": f"plot-only from {Path(path).name}"})):
        metadata = {"Title": "Figure 4 angle scan", **extra_metadata}
        fig.savefig(output, dpi=300 if Path(output).suffix == ".png" else None, metadata=metadata)
    plt.close(fig)
    return output_png, output_pdf


if __name__ == "__main__":
    print("Saved " + ", ".join(map(str, plot_results())))
