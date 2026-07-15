# Numerical codes for monitored-circuit magic and entropy simulations

This repository collects the numerical codes used to produce the main numerical results associated with our work on monitored quantum circuits, steady-state magic, and rotation-angle response.

The current release is organized into three code groups.

## 1. `qubit_2sre_rotation/`

This folder contains the qubit simulation for the steady-state second Rényi entropy (2-SRE) as a function of the rotation angle `theta_M`.

Main files:

- `qubit_2sre_rotation.py`  
  Main simulation script for the monitored qubit protocol. It implements the sequence
  `U_C -> Rx(theta_M) on qubit 0 -> Z measurement -> U_C^dagger`,
  and computes the steady-state 2-SRE for different system sizes.

- `plot_qubit_2sre_rotation.py`  
  Replots the computed 2-SRE data with connecting lines and a fitted small-angle reference slope.

- `replot_qubit_2sre_rotation.py`  
  Produces a relabeled version of the 2-SRE figure with updated notation.

- `qubit_2sre_rotation_data.txt`  
  Tabulated numerical data associated with the plotted 2-SRE curves.

- `qubit_2sre_rotation.png`  
  Example figure output.

This code group is used to:

- compute the qubit monitored-circuit steady-state 2-SRE response;
- verify the small-angle quadratic scaling with `theta_M`;
- provide figure-ready output and tabulated numerical data.

## 2. `qudit_mana_rotation/`

This folder contains the odd-prime qudit simulation for the steady-state Gross mana as a function of the rotation angle `theta_M`.

Main files:

- `qudit_mana_rotation.py`  
  Main qudit monitored-circuit mana simulation.

- `qudit_mana_rotation_strict.py`  
  Version using strict global Clifford sampling over the qudit symplectic group.

- `test_d5_small_theta.py`  
  Dedicated small-angle test for the `d=5` case.

- `mana_strict_clifford.png`  
  Example figure output for steady-state mana versus rotation angle.

This code group is used to:

- study steady-state mana under the monitored qudit protocol;
- compare different odd-prime local dimensions and system sizes;
- verify the weak-angle scaling regime numerically.

Rotation-phase convention used in the qudit mana code:

- for `d = 3`, the diagonal phase data are taken as `tau = (0, 1, -1)`;
- for `d = 5`, the diagonal phase data are taken as `tau = (0, 1, 3, 2, 4)`.

This is the parameter convention used in the local rotation gate for the numerical plots in this folder.

## 3. `four_state_mana_protocol/`

This folder contains a MATLAB implementation for comparing mana dynamics for four different initial-state ensembles under the monitored Clifford protocol.

Main file:

- `Mana_protocol_4states_log2_embedded_0714.m`  
  Self-contained MATLAB script that compares exact Gross mana dynamics, using base-two logarithms, for four initial ensembles at fixed `d=3`, `N=5`, and `theta_M=0.2`.

The four initial ensembles are:

1. global Haar-random state;
2. tensor-product qutrit T-type magic state;
3. GUE-evolved product state;
4. computational-basis product state.

This code group is used to:

- compare how different initial resources affect mana dynamics;
- simulate the monitored Clifford protocol in a public-frame formulation;
- generate figure files and saved numerical data for manuscript use.

## Environment notes

The current release is assembled from working research scripts. Different code groups may require different environments.

Typical dependencies include:

- Python 3
- `numpy`
- `matplotlib`
- `scipy` (optional in some scripts)
- `qiskit` (used in some strict Clifford sampling scripts)
- MATLAB (for the `.m` script)

## Suggested usage

Each code group is intended to be run from within its own folder so that relative paths for numerical data and figure outputs resolve correctly.

## Scope of the release

This repository is intended to present and organize the numerical codes used in the project. A later cleanup may add unified dependency files, command examples, and a more standardized figure/data manifest.
