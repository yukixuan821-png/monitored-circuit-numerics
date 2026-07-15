# Qubit steady-state 2-SRE versus rotation angle

This directory contains the Python code for the qubit monitored-circuit simulation of the steady-state second Rényi entropy (2-SRE) as a function of the rotation angle `theta_M`.

## Files

* `qubit_2sre_rotation.py` contains the main simulation of the monitored qubit protocol and computes the steady-state 2-SRE for different system sizes.
* `plot_qubit_2sre_rotation.py` replots the computed 2-SRE data with connecting lines and a fitted small-angle reference slope.
* `replot_qubit_2sre_rotation.py` produces a relabeled version of the 2-SRE figure with updated notation.
* `qubit_2sre_rotation_data.txt` contains tabulated numerical data associated with the plotted 2-SRE curves.
* `qubit_2sre_rotation.png` is an example figure output.

## Usage

Run `qubit_2sre_rotation.py` from within this directory to generate the numerical data. The plotting scripts can then be used to regenerate figure-ready plots from the computed data.
