"""Qubit weak-magic-injection scan used for Figure 4.

The implementation samples a fresh labelled public Clifford frame at every
monitored cycle.  It works directly with the induced ordered hyperbolic Pauli
pair and therefore does not require Qiskit or a brick-wall approximation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

CASES = {6: (50, 200), 8: (50, 500), 10: (30, 1000), 12: (20, 2000)}
THETA_VALUES = np.logspace(-4, 0, 10)


def _rng(n, theta_index, replica, master_seed=1):
    identity = [4, f"n{n}", theta_index, replica]
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


def sample_public_frame(n, rng):
    while True:
        x = rng.integers(0, 2, size=2*n, dtype=np.int64)
        if np.any(x):
            break
    coeff = np.concatenate((-x[n:], x[:n])) % 2
    pivot = int(np.flatnonzero(coeff)[0])
    z = np.empty_like(x)
    mask = np.arange(2*n) != pivot
    z[mask] = rng.integers(0, 2, size=2*n-1, dtype=np.int64)
    z[pivot] = (1 - int(coeff[mask] @ z[mask])) % 2
    offsets = rng.integers(0, 2, size=2, dtype=np.int64)
    return x, z, int(offsets[0]), int(offsets[1])


def _displacement(state, axis, power=1, offset=0):
    n = axis.size // 2
    x, z = axis[:n] % 2, axis[n:] % 2
    q = _digits(2, n)
    target = _indices((q + power*x) % 2, 2)
    phase = ((-1.0) ** ((power*(q @ z) + power*offset) % 2)).astype(np.complex128)
    phase *= (-1j) ** (power*power*int(x @ z))
    out = np.empty_like(state)
    out[target] = phase * state
    return out


def _spectral_components(state, axis, offset):
    return np.fft.fft(np.stack([_displacement(state, axis, r, offset) for r in range(2)]), axis=0) / 2


def monitored_cycle(state, theta, rng):
    n = state.size.bit_length() - 1
    x, z, x_offset, z_offset = sample_public_frame(n, rng)
    components = _spectral_components(state, x, x_offset)
    angle = np.pi * theta / 8
    state = np.einsum("k,kj->j", [np.exp(-1j*angle), np.exp(1j*angle)], components, optimize=True)
    state /= np.linalg.norm(state)
    components = _spectral_components(state, z, z_offset)
    probabilities = np.maximum(np.sum(np.abs(components)**2, axis=1).real, 0)
    probabilities /= probabilities.sum()
    outcome = int(rng.choice(2, p=probabilities))
    return components[outcome] / np.sqrt(probabilities[outcome])


def _fwht(array):
    step = 1
    while step < array.shape[-1]:
        view = array.reshape(*array.shape[:-1], -1, 2*step)
        left, right = view[..., :step].copy(), view[..., step:2*step].copy()
        view[..., :step], view[..., step:2*step] = left + right, left - right
        step *= 2


def total_s2(state, block_size=64):
    state = np.asarray(state, dtype=np.complex128)
    state /= np.linalg.norm(state)
    size = state.size
    n = size.bit_length() - 1
    indices = np.arange(size, dtype=np.int64)
    total = 0.0
    for begin in range(0, size, block_size):
        shifts = np.arange(begin, min(begin+block_size, size), dtype=np.int64)
        correlations = np.conj(state)[None, :] * state[np.bitwise_xor(shifts[:, None], indices[None, :])]
        _fwht(correlations)
        total += float(np.sum(np.abs(correlations)**4))
    return float(-np.log2(min(total / (2**n), 1.0)))


def run_trajectory(n, theta, cycles, rng):
    state = np.zeros(2**n, dtype=np.complex128)
    state[0] = 1
    for _ in range(cycles):
        state = monitored_cycle(state, theta, rng)
    return total_s2(state)


def run_scan(output="qubit_2sre_rotation_data.txt", cases=CASES, theta_values=THETA_VALUES):
    rows = []
    for n, (replicas, cycles) in cases.items():
        for i, theta in enumerate(theta_values):
            values = [run_trajectory(n, float(theta), cycles, _rng(n, i, r)) for r in range(replicas)]
            mean = float(np.mean(values)); sem = float(np.std(values, ddof=1)/np.sqrt(replicas))
            rows.append((n, cycles, float(theta), replicas, mean, sem))
            print(f"N={n}, theta={theta:.6g}, S2={mean:.8g} +/- {sem:.3g}")
    with open(output, "w", encoding="utf-8", newline="") as stream:
        stream.write("# Terminal total base-two S2; error is the standard error of the mean.\n")
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["N", "T", "theta_M", "Nr", "S2_mean", "S2_sem"])
        writer.writerows(rows)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=None)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(argv)
    folder = Path(__file__).resolve().parent
    output = Path(args.output) if args.output else folder / ("qubit_2sre_rotation_smoke.txt" if args.smoke else "qubit_2sre_rotation_data.txt")
    run_scan(output, {2: (2, 3)} if args.smoke else CASES, np.array([1e-3, 1e-2]) if args.smoke else THETA_VALUES)
    if not args.smoke:
        from plot_qubit_2sre_rotation import plot_results
        print("Saved " + ", ".join(map(str, plot_results(output))))


if __name__ == "__main__":
    main()
