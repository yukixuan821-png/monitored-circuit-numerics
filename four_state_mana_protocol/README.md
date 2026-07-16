# Four-state mana protocol comparison

This directory contains the Matlab code for comparing exact Gross mana dynamics for four initial-state ensembles under the monitored Clifford protocol.

## Files

* `Mana_protocol_4states_log2_embedded_0714.m` contains the full simulation. It compares the following initial ensembles:
  * global Haar-random state;
  * tensor-product qutrit T-type magic states;
  * GUE-evolved product state;
  * computational-basis product state.
* `Mana_0.2_4states.png` is the saved figure output for the four-state comparison.
* `Mana_0.2_4states.pdf` is the corresponding PDF figure.

## Protocol and figure details

The script uses the public-frame equivalent of the monitored Clifford protocol at `d = 3`, `N = 5`, and `theta_M = 0.2`. It applies the labelled random Pauli pair, the qutrit T-type kick, and the Born measurement exactly as documented in the script header.

The four initial-state ensembles are:

* Haar-random states;
* `|T\rangle^{\otimes N}`, with `|T\rangle \propto |0\rangle + e^{2\pi i/9}|1\rangle + e^{-2\pi i/9}|2\rangle`;
* `|\psi\rangle_{\mathrm{GUE}} = \exp(-iH_{\mathrm{GUE}}t)|0\rangle^{\otimes N}` with normalized spectrum in `[-2,2]` and `t=0.1`;
* `|0\rangle^{\otimes N}`.

The current settings are:

* `Tsteps = 2000`;
* `Nr = 1000`;
* `stride = 1`;
* `seed = 1`;
* the same measurement depth for all four initial ensembles.

The current code does not use a separate burn-in or terminal-time-only averaging. It averages the exact mana at each recorded time step over trajectories, so every trajectory is sampled at the same set of times.

The paper caption already captures the main numerical meaning; the README adds that the script also saves the figure files and the numerical data.

## Usage

Run `Mana_protocol_4states_log2_embedded_0714.m` from within this directory. The script produces figure files and saved numerical data for the four-state comparison.
