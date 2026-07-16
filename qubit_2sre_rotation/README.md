# Qubit steady-state 2-SRE versus rotation angle

This directory contains the Python code for the qubit monitored-circuit simulation of the steady-state second Rényi entropy (2-SRE) as a function of the rotation angle `theta_M`.

## Files

* `qubit_2sre_rotation.py` contains the main simulation of the monitored qubit protocol and computes the steady-state 2-SRE for different system sizes.
* `plot_qubit_2sre_rotation.py` replots the computed 2-SRE data with connecting lines and a fitted small-angle reference slope.
* `replot_qubit_2sre_rotation.py` produces a relabeled version of the 2-SRE figure with updated notation.
* `qubit_2sre_rotation_data.txt` contains tabulated numerical data associated with the plotted 2-SRE curves.
* `qubit_2sre_rotation.png` is an example figure output.

## Protocol and figure details

The simulated circuit is the monitored qubit protocol used in the Fig. 6-style comparison: a global random Clifford circuit `U_C`, a single-qubit `Rx(theta_M)` kick on qubit 0, a Z measurement on qubit 0, and the inverse circuit `U_C^\dagger`.

The current scripts use:

* no burn-in; the plotted quantity is the final-time steady-state estimate after `T_max` cycles;
* the same `Nr` and `T_max` for all `theta_M` values within each `N` block;
* standard-error bars, `std / sqrt(Nr)`;
* a small-angle log-log fit over `theta_M < 0.1`.

The trajectory seed is deterministic:

* `seed = N * 100000 + int(theta_M * 10000) + r` in the main script.

The implemented 2-SRE formula is:

* `M2 = N * ln(2) - ln(sum_P |<P>|^4)`;
* `m2 = M2 / N`.

## Usage

Run `qubit_2sre_rotation.py` from within this directory to generate the numerical data. The plotting scripts can then be used to regenerate figure-ready plots from the computed data.
