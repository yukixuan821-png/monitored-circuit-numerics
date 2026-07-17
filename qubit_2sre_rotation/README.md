# Steady-state 2-SRE versus rotation angle in the qubit monitored circuit

This directory contains Python code for studying how the steady-state second Rényi entropy (2-SRE) in a qubit random monitored quantum circuit depends on the local magic-injection angle $\theta_M$. Starting from the computational-basis product initial state, the numerical task scans $\theta_M$ in the small-angle regime, compares the steady-state 2-SRE response for different system sizes, and uses those data to analyze the small-angle weak-magic-injection behavior in the qubit setting.

The simulation uses the qubit monitored circuit in which each monitored cycle applies a global random Clifford scrambling step, an $R_x(\theta_M)$ rotation on the first qubit, a Z measurement on that qubit, and the corresponding inverse Clifford. The script scans several system sizes and several values of $\theta_M$, and extracts the steady-state 2-SRE at the final evolution time.

The current script uses the following system sizes:

* $N = 6, 8, 10, 12$

The rotation angles are taken from a logarithmically spaced scan over

* $\theta_M = 10^{-4} \sim 1$.

The numbers of trajectories and evolution steps for these system sizes are:

* `N = 6`: `Nr = 50`, `T_max = 200`
* `N = 8`: `Nr = 50`, `T_max = 500`
* `N = 10`: `Nr = 30`, `T_max = 1000`
* `N = 12`: `Nr = 20`, `T_max = 2000`

For a fixed $N$, all values of $\theta_M$ use the same `Nr` and `T_max`, so the small-angle response is compared under the same numerical budget. The burn-in length is `0`. Each data point is obtained by averaging the final-time 2-SRE over `Nr` independent trajectories, and the error bars are standard errors $\mathrm{std}/\sqrt{N_r}$. The small-angle power-law slope is fitted over the range $\theta_M < 0.1$.

Each trajectory in the main script uses the deterministic seed

* `seed = N * 100000 + int(theta_M * 10000) + r`

where `r` is the trajectory index.

The 2-SRE observable used in this directory is implemented as

* $M_2 = N \ln 2 - \ln\!\left(\sum_P |\langle P \rangle|^4\right)$
* $m_2 = M_2 / N$

where the sum runs over all Pauli strings. The main script uses an FWHT-based accelerated implementation for this quantity.

The main files in this directory are:

* `qubit_2sre_rotation.py`: main simulation script;
* `plot_qubit_2sre_rotation.py`: replots the numerical results and adds the reference slope line;
* `replot_qubit_2sre_rotation.py`: generates a relabeled version of the figure;
* `qubit_2sre_rotation_data.txt`: tabulated numerical data;
* `qubit_2sre_rotation.png`: example figure output.

## Installation

The code in this directory was tested in a **Python 3.11** environment.

The main dependencies are:

* `numpy`
* `matplotlib`

Optional dependencies include:

* `scipy`, for accelerated Hadamard-matrix-related computation;
* `qiskit`, for strict random-Clifford-related implementations.

## Usage

Run `qubit_2sre_rotation.py` to generate the steady-state 2-SRE data and figure output. The scripts `plot_qubit_2sre_rotation.py` and `replot_qubit_2sre_rotation.py` can be used to replot existing data.
