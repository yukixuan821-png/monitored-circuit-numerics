# Numerical codes for monitored-circuit magic and entropy simulations

This repository contains the numerical code used to illustrate the monitored-circuit dynamics studied in *Invariant Measures and Weak-Magic-Injection Asymptotics in Random Monitored Quantum Circuits*. The repository is organized into three code groups, corresponding to three complementary numerical tasks: steady-state qubit 2-SRE, steady-state odd-prime qudit mana, and exact mana dynamics for four representative initial-state ensembles.

## Model schematic

The monitored-circuit protocol studied in this repository is illustrated below. In each monitored cycle, the circuit applies a fresh global Clifford scrambling step, a local non-Clifford rotation, a local computational-basis measurement on qudit 1, and then the inverse Clifford to return to the original frame.

![Monitored-circuit model schematic](figures/model_schematic.png)

## Repository structure

* `qubit_2sre_rotation`: qubit simulations for the steady-state second Rényi entropy (2-SRE) as a function of the rotation angle `theta_M`.
* `qudit_mana_rotation`: odd-prime qudit simulations for the steady-state Gross mana as a function of `theta_M`, including a strict global-Clifford variant.
* `four_state_mana_protocol`: Matlab code for exact Gross mana dynamics for four initial-state ensembles.

## What each code group computes

### `qubit_2sre_rotation`
This directory computes the steady-state qubit 2-SRE response for several system sizes and small-angle rotation values. It corresponds to the qubit weak-magic-injection asymptotics from the computational-basis product initial state, where the numerical data support the quadratic response law.

### `qudit_mana_rotation`
This directory computes the steady-state Gross mana for odd-prime local dimensions. The main script uses a brick-wall qudit Clifford approximation, while `qudit_mana_rotation_strict.py` uses uniform sampling over the qudit symplectic Clifford group. These simulations support the linear weak-magic-injection response in odd prime dimension.

### `four_state_mana_protocol`
This directory computes the exact Gross mana time evolution for four different initial-state ensembles under the monitored Clifford protocol. It is used as a steady-state diagnostic showing that different initial ensembles approach the same long-time regime.

The four initial-state ensembles are:

* Haar-random states;
* `|T\rangle^{\otimes N}`, with `|T\rangle \propto |0\rangle + e^{2\pi i/9}|1\rangle + e^{-2\pi i/9}|2\rangle`;
* `|\psi\rangle_{\mathrm{GUE}} = \exp(-iH_{\mathrm{GUE}}t)|0\rangle^{\otimes N}` with spectrum normalized to `[-2,2]` and `t=0.1`;
* `|0\rangle^{\otimes N}`.

## Caption-aligned implementation details

The paper captions state the physics and asymptotic interpretation of the plots. The public code additionally fixes the concrete numerical implementation choices summarized below.

### `four_state_mana_protocol`
This script uses `d=3`, `N=5`, `theta_M=0.2`, `Tsteps=2000`, `Nr=1000`, `stride=1`, and `seed=1`.

The current code does not use a separate burn-in. It averages the exact mana at each recorded time step over trajectories, so this is a time-series ensemble average rather than a terminal-time-only average or a time-window average.

All four initial-state ensembles are evolved for the same number of monitored steps and sampled on the same time grid.

### `qudit_mana_rotation`
For each displayed `(d,N)` block, the code uses the same `Nr` and `T_max` for all `theta_M` values. The plotted quantity is obtained from terminal-time sampling after evolution to the chosen numerical stationary regime; the current public scripts do not introduce a separate burn-in counter or a time-window average.

The scripts report standard-error bars `std / sqrt(Nr)` and fit the small-angle slope in log-log scale over `theta_M < 0.1`.

The deterministic seed is

* `seed = N * 100000 + d * 10000 + int(theta * 10000) + r`.

The implemented mana uses phase-point operators and evaluates `log(sum(abs(W.real)))`. The current scripts print the fitted slope but do not save its standard error.

### `qubit_2sre_rotation`
For each displayed `N`, the code uses the same `Nr` and `T_max` for all `theta_M` values. The plotted quantity is obtained from terminal-time sampling after evolution to the chosen numerical stationary regime; the current public scripts do not introduce a separate burn-in counter or a time-window average.

The scripts report standard-error bars `std / sqrt(Nr)` and fit the small-angle slope in log-log scale over `theta_M < 0.1`.

The deterministic seed is

* `seed = N * 100000 + int(theta_M * 10000) + r`.

The implemented 2-SRE formula is

* `M2 = N*ln(2) - ln(sum_P |<P>|^4)`;
* `m2 = M2 / N`.

The current scripts print the fitted slope but do not save its standard error.

## Installation

* `qubit_2sre_rotation` and `qudit_mana_rotation` use Python 3.
* The Python code uses `numpy` and `matplotlib`.
* Some scripts optionally use `scipy`.
* The strict Clifford sampler uses `qiskit`.
* `four_state_mana_protocol` uses Matlab.

## Main entry points

* `qubit_2sre_rotation/qubit_2sre_rotation.py`: computes the steady-state qubit 2-SRE response for different system sizes and rotation angles.
* `qudit_mana_rotation/qudit_mana_rotation.py`: computes the steady-state qudit mana response for odd-prime local dimensions and rotation angles.
* `qudit_mana_rotation/qudit_mana_rotation_strict.py`: computes the same qudit mana response using strict global Clifford sampling over the qudit symplectic group.
* `four_state_mana_protocol/Mana_protocol_4states_log2_embedded_0714.m`: compares exact Gross mana dynamics for four initial-state ensembles in the monitored Clifford protocol.

Run each code group from within its own directory so that relative paths for data and figure outputs resolve correctly.
