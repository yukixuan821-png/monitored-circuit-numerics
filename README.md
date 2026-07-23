# Numerical codes for magic and entropy simulations in random monitored quantum circuits

This repository contains the numerical code for *Invariant Measures and Weak-Magic-Injection Asymptotics in Random Monitored Quantum Circuits*. The code is built around a common random monitored-circuit model and covers the small-angle response of 2-SRE in the qubit setting, the small-angle response of mana in the odd-prime qudit setting, and the exact mana time evolution of representative initial-state ensembles.

## Monitored-circuit model

All three code groups in this repository are based on the same monitored-circuit protocol. One monitored cycle consists of a freshly sampled public Clifford frame, a local magic injection on the first qubit or qudit, a measurement of the conjugate labelled Pauli axis, and the corresponding return to the public frame. The Figure 3 and Figure 4 programs sample the exact protocol-level ordered symplectic Pauli pair, including independent spectral offsets, rather than using a finite-depth brick-wall approximation or a fixed frame across a trajectory.

The public-frame representation makes the equivalence with the original Clifford protocol explicit. For a local prime dimension $d$, write a phase-free Pauli label as

$$
P(a,b)=\bigotimes_{j=1}^{N}X_j^{a_j}Z_j^{b_j},\qquad (a,b)\in\mathbb F_d^{2N}.
$$

If $x=(a_x,b_x)$ and $z=(a_z,b_z)$ are the labels of the kick and measurement axes, the sampled ordered pair satisfies

$$
\langle x,z\rangle_{\mathrm{sp}}=a_x\!\cdot b_z-a_z\!\cdot b_x=1\pmod d.
$$

The labelled operators are $P_X=\omega^{s_X}P(x)$ and $P_Z=\omega^{s_Z}P(z)$, with independent uniformly sampled offsets $s_X,s_Z\in\mathbb F_d$ and $\omega=e^{2\pi i/d}$. Thus $P_ZP_X=\omega P_XP_Z$ under the convention used by the code. Conjugating the original cycle $U_C\to$ local injection along $X_1\to Z_1$ measurement $\to U_C^\dagger$ gives precisely a kick along $P_X=U_C^\dagger X_1U_C$ followed by a measurement along $P_Z=U_C^\dagger Z_1U_C$. The pair, its labels, and its offsets are resampled from one common Clifford frame at every monitored cycle; the kick and measurement are not chosen independently, and a single frame is not reused for an entire trajectory.

The model schematic below represents the common dynamical framework shared by all three code groups: `qubit_2sre_rotation` studies the terminal total base-two 2-SRE in the qubit setting, `qudit_mana_rotation` studies the terminal base-two Gross mana in the odd-prime qudit setting, and `four_state_mana_protocol` studies the exact mana time evolution of different initial states at fixed parameters.

[Monitored-circuit model schematic (vector PDF)](figures/model_schematic.pdf)

## Repository structure

* `qubit_2sre_rotation`: **Python** code. Scans the local magic-injection angle `theta_M` and computes the terminal total base-two stabilizer Rényi entropy for `N = 6, 8, 10, 12`. Its self-contained program saves the tabulated results and generates the PNG and vector PDF figures; `--plot-only` redraws them from the released table.
* `qudit_mana_rotation`: **Python** code. Scans `theta_M` and computes the terminal base-two Gross mana for `(d, N) = (3, 2), (3, 3), (3, 4), (5, 2)`. Its self-contained program saves the tabulated results and generates the PNG and vector PDF figures; `--plot-only` redraws them from the released table.
* `four_state_mana_protocol`: **Matlab** code. Compares the exact Gross mana time evolution of representative initial-state ensembles at fixed parameters, to show how different initial states relax toward a common long-time monitored regime. The supplied $d=3$, $N=5$, $\theta_M=0.2$, $N_r=1000$ run is included in PNG, PDF, EPS, editable FIG, and MAT formats under `four_state_mana_protocol/figures/`.

## Installation

### Python code

The Python programs were tested with **Python 3.11**. The current implementations require Python 3.8 or later, NumPy 1.21 or later, and Matplotlib 3.4 or later.

### Matlab code

The Matlab code in `four_state_mana_protocol` was run in **Matlab R2024b**.

## Main entry points

* `qubit_2sre_rotation/qubit_2sre_rotation.py`: runs the Figure 4 scan, saves the terminal total base-two `S2` data, and generates the PNG and vector PDF outputs; `--plot-only` redraws them from the released table.
* `qudit_mana_rotation/qudit_mana_rotation.py`: runs the Figure 3 scan, saves the terminal base-two mana data, and generates the PNG and vector PDF outputs; `--plot-only` redraws them from the released table.
* `four_state_mana_protocol/Mana_protocol_4states_log2_embedded_0714.m`: compares the exact mana time evolution of representative initial-state ensembles in the monitored protocol.

The released text tables provide the numerical values plotted in Figures 3 and 4 and allow the figures to be checked or redrawn without repeating the full simulations. The committed PNG files are previews; the corresponding PDF files are the vector figures intended for manuscript use. The four-state MAT file contains the time grid, ensemble-mean and representative-trajectory mana arrays, initial-state labels, and run metadata used by the MATLAB figure.
