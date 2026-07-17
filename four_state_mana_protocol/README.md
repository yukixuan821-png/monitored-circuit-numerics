# Exact mana dynamics for four representative initial states

This directory contains Matlab code for comparing the exact Gross mana dynamics of four representative initial states in a random monitored quantum circuit. The numerical task studies how the magic of different initial states evolves with measurement steps at fixed parameters, and whether those states approach a common long-time behavior under the same monitored dynamics.

The script compares the following four initial states: Haar-random states, tensor-product qutrit `T`-type states, GUE-evolved product states, and computational-basis product states. The single-qutrit `T`-type state is taken as

`|T\rangle \propto |0\rangle + e^{2\pi i/9}|1\rangle + e^{-2\pi i/9}|2\rangle`,

and the GUE-evolved state is

`|\psi\rangle_{\mathrm{GUE}} = \exp(-iH_{\mathrm{GUE}}t)|0\rangle^{\otimes N}`,

with the spectrum normalized to `[-2,2]` and `t = 0.1`.

The simulation uses the public-frame-equivalent implementation of the qutrit monitored Clifford protocol at fixed parameters

* `d = 3`
* `N = 5`
* `theta_M = 0.2`

with `Tsteps = 2000`, `Nr = 1000`, `stride = 1`, and `seed = 1`; the burn-in length is `0`. All four initial states are evolved under the same monitored dynamics and on the same time grid, so this simulation directly compares how different initial magic structures evolve in time and how they behave at long times. The script evaluates the exact mana at every recorded time step and averages over all trajectories at that step, producing an ensemble-averaged exact-mana time-evolution curve.

The exact Gross mana used in this directory is defined as

`\mathcal{M}(\psi) = \log_2\|W_\psi\|_1`.

The main files in this directory are:

* `Mana_protocol_4states_log2_embedded_0714.m`: main simulation script;
* `Mana_0.2_4states.png`: figure output for the four-state comparison;
* `Mana_0.2_4states.pdf`: corresponding PDF output.


## Installation

The Matlab code in this directory was run in **Matlab R2024b**.

## Usage

Run `Mana_protocol_4states_log2_embedded_0714.m` to generate the exact mana time-evolution figure for the four representative initial states and the associated data files.
