# Steady-state mana versus rotation angle in the odd-prime qudit monitored circuit

This directory contains Python code for studying how the steady-state Gross mana in an odd-prime qudit random monitored quantum circuit depends on the local magic-injection angle $\theta_M$. The numerical task scans $\theta_M$ in the small-angle regime, compares the steady-state mana response for different local dimensions and system sizes, and uses those data to analyze the small-angle weak-magic-injection behavior in the odd-prime qudit setting.

The simulation uses the odd-prime qudit monitored circuit in which each monitored cycle applies a random qudit Clifford scrambling step, a local non-Clifford rotation $R_X^{(d)}(\theta_M)$ on the first qudit, a computational-basis measurement, and the corresponding inverse Clifford. The main script uses a brick-wall qudit Clifford approximation, while `qudit_mana_rotation_strict.py` uses strict global sampling over the qudit symplectic Clifford group; both implement the same physical protocol.

The current script uses the following parameter blocks:

* `(d, N) = (3, 2)`
* `(d, N) = (3, 3)`
* `(d, N) = (3, 4)`
* `(d, N) = (5, 2)`

The rotation angles are taken from a logarithmically spaced scan over

* $\theta_M = 10^{-4} \sim 1$.

All four parameter blocks use `T_max = 500`. The corresponding trajectory numbers are:

* `(3, 2)`: `Nr = 50`
* `(3, 3)`: `Nr = 50`
* `(3, 4)`: `Nr = 30`
* `(5, 2)`: `Nr = 50`

For a fixed $(d, N)$, all values of $\theta_M$ use the same `Nr` and `T_max`, so the angle dependence is compared under the same numerical budget. The burn-in length is `0`. Each data point is obtained by averaging the final-time mana over `Nr` independent trajectories, and the error bars are standard errors $\mathrm{std}/\sqrt{N_r}$. The small-angle power-law slope is fitted over the range $\theta_M < 0.1$.

Each trajectory in the main script uses the deterministic seed

* `seed = N * 100000 + d * 10000 + int(theta * 10000) + r`

where `r` is the trajectory index.

The local non-Clifford rotation uses the same phase-lift convention as in the paper. For $d = 3$,

* $\tau_3(0) = 0$
* $\tau_3(1) = 1$
* $\tau_3(2) = -1$

For $d = 5$, the current code fixes $a = 1$ and takes $\tau_{1,5}(k)$ to be the standard representative of $k^3 \in \mathbb{F}_5$, namely

* $\tau_{1,5}(0) = 0$
* $\tau_{1,5}(1) = 1$
* $\tau_{1,5}(2) = 3$
* $\tau_{1,5}(3) = 2$
* $\tau_{1,5}(4) = 4$

Gross mana in this directory is computed through the phase-point-operator representation, implemented as

* $W(u,v) = \langle \psi|A_{u,v}|\psi \rangle / d^N$
* $\mathcal{M}(\psi) = \log \sum_{u,v} |\mathrm{Re}\, W(u,v)|$

The main files in this directory are:

* `qudit_mana_rotation.py`: main simulation script;
* `qudit_mana_rotation_strict.py`: strict global-Clifford sampling version;
* `test_d5_small_theta.py`: focused small-angle test for `d = 5`;
* `mana_strict_clifford.png`: example figure output.

## Installation

The code in this directory was tested in a **Python 3.11** environment.

The main dependencies are:

* `numpy`
* `matplotlib`

Some implementations also use:

* `qiskit`, for strict global-Clifford-sampling-related code.

## Usage

Run `qudit_mana_rotation.py` or `qudit_mana_rotation_strict.py` to generate the steady-state mana data and figure output. The script `test_d5_small_theta.py` can be used for a focused small-angle test in the `d = 5` case.
