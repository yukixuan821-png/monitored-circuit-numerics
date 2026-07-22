# Steady-state 2-SRE versus rotation angle in the qubit monitored circuit

This directory contains Python code for studying how the qubit stabilizer Rényi entropy depends on the local magic-injection angle $\theta_M$. Starting from the computational-basis product state $|0\rangle^{\otimes N}$, the numerical task scans $\theta_M$ in the small-angle regime and compares the terminal response for several system sizes.

Each monitored cycle independently samples a labelled public Clifford frame, applies the local injection projectively equivalent to $\exp(-i\pi\theta_M P/8)$, and measures the conjugate labelled Pauli axis. The exact protocol-level sampler uses a uniformly sampled ordered hyperbolic Pauli pair with independent spectral offsets. It has no Qiskit-dependent or finite-depth brick-wall fallback.

The current script uses the following system sizes:

* $N = 6, 8, 10, 12$.

The rotation angles are taken from a ten-point logarithmically spaced scan over

* $\theta_M = 10^{-4} \sim 1$.

The numbers of trajectories and evolution steps for these system sizes are:

* `N = 6`: `Nr = 50`, `T_max = 200`
* `N = 8`: `Nr = 50`, `T_max = 500`
* `N = 10`: `Nr = 30`, `T_max = 1000`
* `N = 12`: `Nr = 20`, `T_max = 2000`

For a fixed $N$, all values of $\theta_M$ use the same `Nr` and `T_max`, so the angle dependence is compared under the same numerical budget. Each data point is the ensemble mean of the terminal value over `Nr` independent trajectories, and the error bars are standard errors $\mathrm{std}/\sqrt{N_r}$.

The observable is the total base-two stabilizer Rényi entropy

$$
S_2=-\log_2\left(2^{-N}\sum_{P\in\mathcal P_N}|\langle P\rangle|^4\right),
$$

not an entropy density. The code evaluates this quantity using a blocked in-place Walsh–Hadamard transform.

The main files in this directory are:

* `qubit_2sre_rotation.py`: main simulation script;
* `plot_qubit_2sre_rotation.py`: redraws the figure from the released table;
* `qubit_2sre_rotation_data.txt`: final tabulated means and standard errors;
* `qubit_2sre_rotation.png`: final plot preview;
* `qubit_2sre_rotation.pdf`: final vector figure.

## Installation

The code was tested with **Python 3.11**. It requires Python 3.8 or later and NumPy 1.21 or later. Matplotlib 3.4 or later is required for the default full run and for plotting, but not for `--smoke`. SciPy and Qiskit are not required.

## Usage

Run `python qubit_2sre_rotation.py` to perform the full scan, save `qubit_2sre_rotation_data.txt`, and generate both figure formats. Run `python plot_qubit_2sre_rotation.py` to redraw the PNG and PDF directly from the released table without rerunning the simulation.

For a small numerical check, run `python qubit_2sre_rotation.py --smoke`. This writes `qubit_2sre_rotation_smoke.txt` and leaves the released paper table and figures unchanged. An explicit `--output` path may be supplied when a different data filename is desired.
