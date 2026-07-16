# Qudit steady-state mana versus rotation angle

This directory contains the Python code for odd-prime qudit simulations of the steady-state Gross mana under the monitored-circuit protocol.

## Files

* `qudit_mana_rotation.py` contains the main qudit mana simulation.
* `qudit_mana_rotation_strict.py` contains the version with strict global Clifford sampling over the qudit symplectic group.
* `test_d5_small_theta.py` contains a focused small-angle test for the case `d = 5`.
* `mana_strict_clifford.png` is an example figure output for steady-state mana versus rotation angle.

## Protocol and figure details

The main script simulates the monitored qudit protocol with:

* a brick-wall qudit Clifford circuit;
* a local non-Clifford rotation `R_X^{(d)}(theta_M)` on qudit 0;
* a computational-basis measurement on qudit 0;
* the inverse Clifford circuit.

The strict variant replaces the approximate brick-wall Clifford draw with uniform sampling over the symplectic Clifford group.

The current scripts use:

* no burn-in; the plotted quantity is the final-time steady-state estimate after `T_max` cycles;
* the same `Nr` and `T_max` for all `theta_M` values within each `(d, N)` block;
* standard-error bars, `std / sqrt(Nr)`;
* a small-angle log-log fit over `theta_M < 0.1`.

The trajectory seed is deterministic:

* `seed = N * 100000 + d * 10000 + int(theta * 10000) + r` in the main script.

The implemented mana uses phase-point operators and the code evaluates `log(sum(abs(W.real)))`.

The local rotation gate uses the following diagonal phase convention:

* `d = 3`: `tau = (0, 1, -1)`
* `d = 5`: `tau = (0, 1, 3, 2, 4)`

## Usage

Run `qudit_mana_rotation.py` or `qudit_mana_rotation_strict.py` from within this directory to generate the steady-state mana data and the corresponding plots.
