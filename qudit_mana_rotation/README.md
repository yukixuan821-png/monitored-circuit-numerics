# Qudit steady-state mana versus rotation angle

This directory contains the Python code for odd-prime qudit simulations of the steady-state Gross mana under the monitored-circuit protocol.

## Files

* `qudit_mana_rotation.py` contains the main qudit mana simulation.
* `qudit_mana_rotation_strict.py` contains the version with strict global Clifford sampling over the qudit symplectic group.
* `test_d5_small_theta.py` contains a focused small-angle test for the case `d = 5`.
* `mana_strict_clifford.png` is an example figure output for steady-state mana versus rotation angle.

The local rotation gate uses the following diagonal phase convention:

* `d = 3`: `tau = (0, 1, -1)`
* `d = 5`: `tau = (0, 1, 3, 2, 4)`

## Usage

Run `qudit_mana_rotation.py` or `qudit_mana_rotation_strict.py` from within this directory to generate the steady-state mana data and the corresponding plots.
