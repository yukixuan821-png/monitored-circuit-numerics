# Steady-state mana versus rotation angle in the odd-prime qudit monitored circuit

This directory contains Python code for studying how the Gross mana in an odd-prime qudit random monitored quantum circuit depends on the local magic-injection angle $\theta_M$. Starting from the computational-basis product state $|0\rangle^{\otimes N}$, the numerical task scans $\theta_M$ in the small-angle regime and compares the terminal mana response for different local dimensions and system sizes.

Each monitored cycle independently samples a labelled public Clifford frame: a uniformly sampled ordered hyperbolic Pauli pair with independent spectral offsets. The local non-Clifford injection is followed by measurement of the conjugate labelled Pauli axis. This exact protocol-level implementation replaces both the former finite-depth brick-wall approximation and the former single-frame-per-trajectory implementation.

The current script uses the following parameter blocks:

* `(d, N) = (3, 2)`
* `(d, N) = (3, 3)`
* `(d, N) = (3, 4)`
* `(d, N) = (5, 2)`

The rotation angles are taken from a ten-point logarithmically spaced scan over

* $\theta_M = 10^{-4} \sim 1$.

All four parameter blocks use `T_max = 500`. The corresponding trajectory numbers are:

* `(3, 2)`: `Nr = 50`
* `(3, 3)`: `Nr = 50`
* `(3, 4)`: `Nr = 30`
* `(5, 2)`: `Nr = 50`

For a fixed $(d,N)$, all values of $\theta_M$ use the same `Nr` and `T_max`, so the angle dependence is compared under the same numerical budget. Here `T_max` is the finite equilibration horizon preceding the observation. Each independent trajectory contributes one mana value evaluated at the terminal time `T_max`; no post-equilibration time-window average is taken. Each data point is the ensemble mean of these `Nr` terminal values, and the error bars are standard errors $\mathrm{std}(\text{terminal values},\mathrm{ddof}=1)/\sqrt{N_r}$.

For $d=3$, the injection phases are

$$
\left(1,e^{+2\pi i\theta_M/9},e^{-2\pi i\theta_M/9}\right).
$$

For $d=5$, the code fixes $a=1$ and uses the standard representatives of $k^3\in\mathbb F_5$,

$$
\left(\tau_{1,5}(0),\tau_{1,5}(1),\tau_{1,5}(2),\tau_{1,5}(3),\tau_{1,5}(4)\right)
=(0,1,3,2,4).
$$

Gross mana is computed in bits as

$$
\mathcal M(\psi)=\log_2\lVert W_\psi\rVert_1,
$$

using a pure-state multidimensional FFT kernel.

The dashed line is a scaling guide proportional to $\theta_M$, corresponding to slope one on the log--log axes. Its normalization is chosen only for visual comparison and is not obtained by fitting the numerical data. The weak-injection results are compared with the predicted linear scaling $\overline{\mathcal M}_T(\theta_M)\propto\theta_M$.

The main files in this directory are:

* `qudit_mana_rotation.py`: self-contained simulation and plotting program;
* `qudit_mana_rotation_data.txt`: final tabulated means and standard errors;
* `mana_strict_clifford.png`: final plot preview;
* `mana_strict_clifford.pdf`: final vector figure.

## Installation

The code was tested with **Python 3.11**. It requires Python 3.8 or later, NumPy 1.21 or later, and Matplotlib 3.4 or later.

## Usage

Run `python qudit_mana_rotation.py` to perform the full scan, save `qudit_mana_rotation_data.txt`, and generate both figure formats. Run `python qudit_mana_rotation.py --plot-only` to redraw the PNG and PDF directly from the released table without rerunning the simulation.
