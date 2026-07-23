# Steady-state 2-SRE versus rotation angle in the qubit monitored circuit

This directory contains Python code for studying how the qubit stabilizer Rényi entropy depends on the local magic-injection angle $\theta_M$. Starting from the computational-basis product state $|0\rangle^{\otimes N}$, the numerical task scans $\theta_M$ in the small-angle regime and compares the terminal response for several system sizes.

Each monitored cycle independently samples a labelled public Clifford frame, applies the local injection projectively equivalent to $\exp(-i\pi\theta_M P/8)$, and measures the conjugate labelled Pauli axis. The exact protocol-level sampler uses a uniformly sampled ordered hyperbolic Pauli pair with independent spectral offsets.

The current script uses the following system sizes:

* $N = 6, 8, 10, 12$.

The rotation angles are taken from a ten-point logarithmically spaced scan over

* $\theta_M = 10^{-4} \sim 1$.

The numbers of trajectories and evolution steps for these system sizes are:

* `N = 6`: `Nr = 50`, `T_max = 200`
* `N = 8`: `Nr = 50`, `T_max = 500`
* `N = 10`: `Nr = 30`, `T_max = 1000`
* `N = 12`: `Nr = 20`, `T_max = 2000`

For a fixed $N$, all values of $\theta_M$ use the same `Nr` and `T_max`, so the angle dependence is compared under the same numerical budget. Here `T_max` is the finite equilibration horizon preceding the observation. Each independent trajectory contributes one $S_2$ value evaluated at the terminal time `T_max`; no post-equilibration time-window average is taken. Each data point is the ensemble mean of these `Nr` terminal values, and the error bars are standard errors $\mathrm{std}(\text{terminal values},\mathrm{ddof}=1)/\sqrt{N_r}$.

The observable is the total base-two stabilizer Rényi entropy

$$
S_2=-\log_2\left(2^{-N}\sum_{P\in\mathcal P_N}|\langle P\rangle|^4\right),
$$

not an entropy density. The code evaluates this quantity using a blocked in-place Walsh–Hadamard transform.

In the figure, the four colored curves correspond to $N=6$, $8$, $10$, and $12$. Every marker is the ensemble mean of one terminal total base-two $S_2$ value from each independent trajectory, and its vertical error bar is the corresponding SEM. The dashed line is a scaling guide proportional to $\theta_M^2$, corresponding to slope two on the log--log axes. Its normalization is chosen only for visual comparison and is not obtained by fitting the numerical data. The weak-injection results are compared with the predicted stationary scaling $S_2(\theta_M)\propto\theta_M^2$ used in the manuscript.

The released `qubit_2sre_rotation_data.txt` table has one row for each parameter point and the columns:

| Column | Meaning |
|---|---|
| `N` | number of qubits |
| `T` | terminal evolution horizon |
| `theta_M` | local magic-injection angle |
| `Nr` | number of independent trajectories |
| `S2_mean` | mean terminal total base-two $S_2$ |
| `S2_sem` | standard error of the terminal mean |

The table contains terminal aggregate statistics rather than complete time series or individual trajectory records.

The main files in this directory are:

* `qubit_2sre_rotation.py`: self-contained simulation and plotting program;
* `qubit_2sre_rotation_data.txt`: final tabulated means and standard errors;
* `qubit_2sre_rotation.png`: final plot preview;
* `qubit_2sre_rotation.pdf`: final vector figure.

## Installation

The code was tested with **Python 3.11**. It requires Python 3.8 or later, NumPy 1.21 or later, and Matplotlib 3.4 or later.

## Usage

Run `python qubit_2sre_rotation.py` to perform the full scan, save `qubit_2sre_rotation_data.txt`, and generate both figure formats. Run `python qubit_2sre_rotation.py --plot-only` to redraw the PNG and PDF directly from the released table without rerunning the simulation.
