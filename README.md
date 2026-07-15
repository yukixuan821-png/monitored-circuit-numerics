# Numerical codes for monitored-circuit magic and entropy simulations

This repository provides Python and Matlab code for numerical simulations of monitored quantum circuits, with emphasis on steady-state magic and entropy under local rotation-angle perturbations. The current implementation contains three code groups: a qubit simulation for the steady-state second Rényi entropy (2-SRE), an odd-prime qudit simulation for the steady-state Gross mana, and a Matlab script comparing mana dynamics for four different initial-state ensembles.

## Structure of repository

* `qubit_2sre_rotation` contains the qubit code for the steady-state 2-SRE as a function of the rotation angle `theta_M`.
* `qudit_mana_rotation` contains the odd-prime qudit code for the steady-state Gross mana as a function of the rotation angle `theta_M`.
* `four_state_mana_protocol` contains the Matlab code for comparing exact Gross mana dynamics for four initial-state ensembles.

## Installation

* The code in `qubit_2sre_rotation` and `qudit_mana_rotation` has been developed using Python 3.
* In addition, the Python code uses `numpy` and `matplotlib`.
* Some scripts optionally use `scipy`.
* The strict Clifford sampling scripts use `qiskit`.
* The code in `four_state_mana_protocol` has been developed using Matlab.

## Usage

Each directory contains its own README with a short overview of the files. The main entry points are:

* `qubit_2sre_rotation/qubit_2sre_rotation.py`: computes the steady-state qubit 2-SRE response for different system sizes and rotation angles.
* `qudit_mana_rotation/qudit_mana_rotation.py`: computes the steady-state qudit mana response for odd-prime local dimensions and rotation angles.
* `qudit_mana_rotation/qudit_mana_rotation_strict.py`: computes the same qudit mana response using strict global Clifford sampling over the qudit symplectic group.
* `four_state_mana_protocol/Mana_protocol_4states_log2_embedded_0714.m`: compares exact Gross mana dynamics for four initial-state ensembles in the monitored Clifford protocol.

Run each code group from within its own directory so that relative paths for data and figure outputs resolve correctly.
