# Qudit steady-state mana versus rotation angle

This folder contains the odd-prime qudit simulations for steady-state Gross mana under the monitored-circuit protocol.

Main files:
- `qudit_mana_rotation.py`: main mana simulation
- `qudit_mana_rotation_strict.py`: strict global Clifford sampling version
- `test_d5_small_theta.py`: focused small-angle test for `d=5`

Figure:
- `mana_strict_clifford.png`

Rotation-phase convention used in the local gate:
- `d = 3`: `tau = (0, 1, -1)`
- `d = 5`: `tau = (0, 1, 3, 2, 4)`

These are the diagonal phase parameters used in the local rotation gate for the numerical curves in this folder.
