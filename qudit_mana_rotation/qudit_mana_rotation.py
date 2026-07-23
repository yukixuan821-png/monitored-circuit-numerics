"""Odd-prime qudit mana weak-magic-injection angle scan.

Each monitored cycle independently samples the labelled public Clifford frame,
including both spectral offsets.  Gross mana is reported in bits.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

CASES = {(3, 2): (50, 500), (3, 3): (50, 500), (3, 4): (30, 500), (5, 2): (50, 500)}
THETA_VALUES = np.logspace(-4, 0, 10)


def _rng(d, n, theta_index, replica, master_seed=1):
    identity = [3, f"d{d}_n{n}", theta_index, replica]
    words = np.frombuffer(hashlib.sha256(json.dumps(identity, separators=(",", ":")).encode()).digest(), dtype="<u4")
    seed = np.random.SeedSequence([master_seed, *map(int, words)])
    return np.random.Generator(np.random.PCG64DXSM(seed))


def _digits(d, n):
    values = np.arange(d**n, dtype=np.int64)[:, None]
    powers = d ** np.arange(n - 1, -1, -1, dtype=np.int64)
    return (values // powers) % d


def _indices(labels, d):
    powers = d ** np.arange(labels.shape[-1] - 1, -1, -1, dtype=np.int64)
    return np.asarray(labels @ powers, dtype=np.int64)


def sample_public_frame(n, d, rng):
    while True:
        x = rng.integers(0, d, size=2*n, dtype=np.int64)
        if np.any(x):
            break
    coeff = np.concatenate((-x[n:], x[:n])) % d
    pivot = int(np.flatnonzero(coeff)[0])
    z = np.empty_like(x)
    mask = np.arange(2*n) != pivot
    z[mask] = rng.integers(0, d, size=2*n-1, dtype=np.int64)
    remainder = int(coeff[mask] @ z[mask]) % d
    z[pivot] = ((1-remainder) * pow(int(coeff[pivot]), -1, d)) % d
    offsets = rng.integers(0, d, size=2, dtype=np.int64)
    return x, z, int(offsets[0]), int(offsets[1])


def _displacement(state, axis, d, power=1, offset=0):
    n = axis.size // 2
    x, z = axis[:n] % d, axis[n:] % d
    r = power % d
    q = _digits(d, n)
    target = _indices((q + r*x) % d, d)
    omega = np.exp(2j*np.pi/d); tau = -np.exp(1j*np.pi/d)
    phase = omega ** ((r*(q @ z) + r*offset) % d)
    phase *= tau ** (r*r*int(x @ z))
    out = np.empty_like(state)
    out[target] = phase * state
    return out


def _spectral_components(state, axis, d, offset):
    orbit = np.stack([_displacement(state, axis, d, r, offset) for r in range(d)])
    return np.fft.fft(orbit, axis=0) / d


def injection_phases(d, theta):
    if d == 3:
        return np.array([1, np.exp(2j*np.pi*theta/9), np.exp(-2j*np.pi*theta/9)])
    labels = np.arange(d, dtype=np.int64)
    cubic_lift = (labels**3) % d
    return np.exp(-2j*np.pi*theta*cubic_lift/d)


def monitored_cycle(state, d, theta, rng):
    n = round(np.log(state.size)/np.log(d))
    x, z, x_offset, z_offset = sample_public_frame(n, d, rng)
    components = _spectral_components(state, x, d, x_offset)
    state = np.einsum("k,kj->j", injection_phases(d, theta), components, optimize=True)
    state /= np.linalg.norm(state)
    components = _spectral_components(state, z, d, z_offset)
    probabilities = np.maximum(np.sum(np.abs(components)**2, axis=1).real, 0)
    probabilities /= probabilities.sum()
    outcome = int(rng.choice(d, p=probabilities))
    return components[outcome] / np.sqrt(probabilities[outcome])


def gross_mana(state, d, q_block=16):
    """Return log2 ||W||_1 using the pure-state multidimensional FFT kernel."""
    n = round(np.log(state.size)/np.log(d)); dimension = d**n
    labels = _digits(d, n); half = (pow(2, -1, d)*labels) % d
    result = np.empty((dimension, dimension), dtype=float)
    axes = tuple(range(1, n+1)); fft_shape = (d,)*n
    for begin in range(0, dimension, q_block):
        q = labels[begin:begin+q_block, None, :]
        plus = _indices((q + half[None, :, :]) % d, d)
        minus = _indices((q - half[None, :, :]) % d, d)
        correlations = state[plus] * np.conj(state[minus])
        transformed = np.fft.fftn(correlations.reshape((-1,)+fft_shape), axes=axes) / dimension
        result[begin:begin+q.shape[0]] = transformed.reshape(q.shape[0], dimension).real
    return float(np.log2(np.sum(np.abs(result))))


def run_trajectory(d, n, theta, cycles, rng):
    state = np.zeros(d**n, dtype=np.complex128); state[0] = 1
    for _ in range(cycles):
        state = monitored_cycle(state, d, theta, rng)
    return gross_mana(state, d)


def plot_results(data_path="qudit_mana_rotation_data.txt", output_png=None, output_pdf=None):
    data_path = Path(data_path)
    output_png = data_path.with_name("mana_strict_clifford.png") if output_png is None else Path(output_png)
    output_pdf = data_path.with_name("mana_strict_clifford.pdf") if output_pdf is None else Path(output_pdf)
    with open(data_path, encoding="utf-8") as stream:
        reader = csv.DictReader((line for line in stream if not line.startswith("#")), delimiter="\t")
        rows = [{key: float(value) for key, value in row.items()} for row in reader]
    colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]
    markers = ["o", "s", "^", "D"]
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix", "font.size": 13, "axes.labelsize": 13, "legend.fontsize": 12, "xtick.labelsize": 13, "ytick.labelsize": 13, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.17, top=0.90)
    all_means = []
    for color, marker, (d, n) in zip(colors, markers, ((3, 2), (3, 3), (3, 4), (5, 2))):
        group = sorted(
            (row for row in rows if int(row["d"]) == d and int(row["N"]) == n),
            key=lambda row: row["theta_M"],
        )
        theta = np.array([row["theta_M"] for row in group])
        mean = np.array([row["mana_mean"] for row in group])
        sem = np.array([row["mana_sem"] for row in group])
        if np.any(mean <= 0):
            raise ValueError("log plot requires strictly positive aggregate means")
        lower_endpoint = np.maximum.reduce((mean-sem, mean*1e-3, np.full_like(mean, np.finfo(float).tiny)))
        asymmetric_yerr = np.vstack((mean-lower_endpoint, sem))
        ax.errorbar(theta, mean, yerr=asymmetric_yerr, color=color, marker=marker, linewidth=2, markersize=6, capsize=2, label=rf"$d={d},\ N={n}$")
        all_means.extend(mean)
    theta_grid = np.unique([row["theta_M"] for row in rows])
    anchor_theta = theta_grid[len(theta_grid)//2]
    anchor_y = float(np.median(all_means))
    ax.plot(theta_grid, anchor_y*(theta_grid/anchor_theta), color="#52514e", linestyle="--", linewidth=1.5, label=r"slope $1$")
    ax.set(xscale="log", yscale="log", xlabel=r"$\theta_M$", ylabel=r"$\overline{\mathcal{M}}(\theta_M)$", title="Steady-state Mana vs rotation angle")
    ax.legend(frameon=False)
    ax.grid(True, which="major", color="#e1e0d9", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    for path, extra_metadata in ((output_png, {"Description": f"plot-only from {data_path.name}"}), (output_pdf, {"Subject": f"plot-only from {data_path.name}"})):
        metadata = {"Title": "Odd-prime qudit mana angle scan", **extra_metadata}
        fig.savefig(path, dpi=300 if Path(path).suffix == ".png" else None, metadata=metadata)
    plt.close(fig)
    return output_png, output_pdf


def run_scan(output="qudit_mana_rotation_data.txt", cases=CASES, theta_values=THETA_VALUES):
    rows = []
    for (d, n), (replicas, cycles) in cases.items():
        for i, theta in enumerate(theta_values):
            values = [run_trajectory(d, n, float(theta), cycles, _rng(d, n, i, r)) for r in range(replicas)]
            mean = float(np.mean(values)); sem = float(np.std(values, ddof=1)/np.sqrt(replicas))
            rows.append((d, n, cycles, float(theta), replicas, mean, sem))
            print(f"d={d}, N={n}, theta={theta:.6g}, mana={mean:.8g} +/- {sem:.3g}")
    with open(output, "w", encoding="utf-8", newline="") as stream:
        stream.write("# Terminal base-two Gross mana; error is the standard error of the mean.\n")
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["d", "N", "T", "theta_M", "Nr", "mana_mean", "mana_sem"])
        writer.writerows(rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=None)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--plot-only", action="store_true", help="rebuild PNG and vector PDF from the data table")
    args = parser.parse_args(argv)
    folder = Path(__file__).resolve().parent
    if args.plot_only:
        data_path = Path(args.output) if args.output else folder / "qudit_mana_rotation_data.txt"
        print("Saved " + ", ".join(map(str, plot_results(data_path))))
    else:
        output = Path(args.output) if args.output else folder / ("qudit_mana_rotation_smoke.txt" if args.smoke else "qudit_mana_rotation_data.txt")
        run_scan(output, {(3, 2): (2, 3)} if args.smoke else CASES, np.array([1e-3, 1e-2]) if args.smoke else THETA_VALUES)
        if not args.smoke:
            print("Saved " + ", ".join(map(str, plot_results(output))))


if __name__ == "__main__":
    main()
