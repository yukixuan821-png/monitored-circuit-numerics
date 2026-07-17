# Numerical codes for magic and entropy simulations in random monitored quantum circuits

This repository contains the numerical code for *Invariant Measures and Weak-Magic-Injection Asymptotics in Random Monitored Quantum Circuits*. The code is built around a common random monitored-circuit model and covers the small-angle response of 2-SRE in the qubit setting, the small-angle response of mana in the odd-prime qudit setting, and the exact mana time evolution of representative initial-state ensembles.

## Monitored-circuit model

All three code groups in this repository are based on the same monitored-circuit protocol. One monitored cycle consists of the following steps:

1. Apply a global random Clifford scrambling step.
2. Apply a local magic injection on the first qubit/qudit.
3. Perform a local computational-basis measurement at that site.
4. Apply the corresponding inverse Clifford to return to the original frame.
5. Repeat the cycle with a freshly sampled random Clifford.

The model schematic below represents the common dynamical framework shared by all three code groups: `qubit_2sre_rotation` studies steady-state 2-SRE in the qubit setting, `qudit_mana_rotation` studies steady-state Gross mana in the odd-prime qudit setting, and `four_state_mana_protocol` studies the exact mana time evolution of different initial states at fixed parameters.

![Monitored-circuit model schematic](figures/model_schematic.png)

## Repository structure

* `qubit_2sre_rotation`: **Python** code. Scans the local magic-injection angle `theta_M` and computes the steady-state second Rényi entropy (2-SRE) for different system sizes in the qubit monitored circuit, to study the small-angle weak-magic-injection response in the qubit setting.
* `qudit_mana_rotation`: **Python** code. Scans `theta_M` and computes the steady-state Gross mana for different local dimensions and system sizes in the odd-prime qudit monitored circuit, to study the small-angle weak-magic-injection response in the odd-prime qudit setting. This directory also includes a strict global-Clifford sampling version.
* `four_state_mana_protocol`: **Matlab** code. Compares the exact Gross mana time evolution of representative initial-state ensembles at fixed parameters, to show how different initial states relax toward a common long-time monitored regime.

## Installation

### Python code

The Python code in this repository was tested in a **Python 3.11** environment. The public Python scripts primarily use the following packages:

* `numpy`
* `matplotlib`

Some scripts also use:

* `scipy`, for accelerated computation in the qubit 2-SRE code;
* `qiskit`, for strict or uniform Clifford-sampling-related implementations.

### Matlab code

The Matlab code in `four_state_mana_protocol` was run in **Matlab R2024b**.

## Main entry points

* `qubit_2sre_rotation/qubit_2sre_rotation.py`: computes the steady-state qubit 2-SRE for different system sizes and rotation angles `theta_M`.
* `qudit_mana_rotation/qudit_mana_rotation.py`: computes the steady-state qudit mana for odd-prime local dimensions, different system sizes, and rotation angles `theta_M`.
* `qudit_mana_rotation/qudit_mana_rotation_strict.py`: computes the same qudit mana data using strict global Clifford sampling.
* `four_state_mana_protocol/Mana_protocol_4states_log2_embedded_0714.m`: compares the exact mana time evolution of representative initial-state ensembles in the monitored protocol.
