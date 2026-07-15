"""
Qudit monitored circuit steady-state Mana vs rotation angle.
Protocol same as Fig.6 but for qudits (odd prime d).
Expected: linear response Mana ~ theta at small angles.
Route A: brick-wall approximation for qudit Clifford gates.
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from itertools import product

# ==========================================
# 1. Basic qudit operators
# ==========================================

def omega(d, k):
    return np.exp(2j * np.pi * k / d)


def X_matrix(d):
    X = np.zeros((d, d), dtype=complex)
    for j in range(d):
        X[(j + 1) % d, j] = 1.0
    return X


def Z_matrix(d):
    Z = np.zeros((d, d), dtype=complex)
    for j in range(d):
        Z[j, j] = omega(d, j)
    return Z


def fourier_matrix(d):
    """QFT: F|j> = (1/sqrt(d)) sum_k omega^{jk} |k>"""
    F = np.zeros((d, d), dtype=complex)
    for j in range(d):
        for k in range(d):
            F[k, j] = omega(d, j * k) / np.sqrt(d)
    return F


def phase_matrix(d):
    """Phase gate P|j> = omega^{j(j-1)/2} |j> (Clifford generator for odd d)"""
    P = np.zeros((d, d), dtype=complex)
    for j in range(d):
        P[j, j] = omega(d, j * (j - 1) / 2)
    return P


# ==========================================
# 2. Mana computation (from test_mana.py)
# ==========================================

def D_matrix(d, p, q):
    """Displacement operator w(p,q) = tau^{pq} X^p Z^q, tau = -e^{i pi / d}"""
    X = X_matrix(d)
    Z = Z_matrix(d)
    Xp = np.linalg.matrix_power(X, p)
    Zq = np.linalg.matrix_power(Z, q)
    tau = -np.exp(1j * np.pi / d)
    return Xp @ Zq * (tau ** (p * q))


def precompute_phase_point_operators(d, n):
    A_single = {}
    for u in range(d):
        for v in range(d):
            A = np.zeros((d, d), dtype=complex)
            for p in range(d):
                for q in range(d):
                    A += omega(d, u * q - v * p) * D_matrix(d, p, q)
            A_single[(u, v)] = A / d
    return A_single


def mana(psi, d, A_single=None):
    D = len(psi)
    n = int(round(np.log(D) / np.log(d)))
    if A_single is None:
        A_single = precompute_phase_point_operators(d, n)

    W = np.zeros(d ** (2 * n), dtype=complex)
    idx = 0
    for phase_point in product(range(d), repeat=2 * n):
        A = np.eye(1, dtype=complex)
        for i in range(n):
            u_i = phase_point[2 * i]
            v_i = phase_point[2 * i + 1]
            A = np.kron(A, A_single[(u_i, v_i)])
        W[idx] = np.vdot(psi, A @ psi) / D
        idx += 1

    return np.log(np.sum(np.abs(W.real)))


# ==========================================
# 3. Circuit operations
# ==========================================

def apply_single_qudit(psi, U, i, N, d):
    psi_r = psi.reshape([d] * N)
    psi_r = np.moveaxis(psi_r, i, 0)
    psi_r = np.tensordot(U, psi_r, axes=([1], [0]))
    psi_r = np.moveaxis(psi_r, 0, i)
    return psi_r.reshape(-1)


def apply_SUM(psi, c, t, N, d):
    """Generalized CNOT: SUM|c,t> = |c, c+t mod d>"""
    psi_r = psi.reshape([d] * N)
    psi_r = np.moveaxis(psi_r, [c, t], [0, 1])
    new_psi = np.zeros_like(psi_r)
    for a in range(d):
        for b in range(d):
            new_psi[a, b, ...] = psi_r[a, (b - a) % d, ...]
    psi_r = new_psi
    psi_r = np.moveaxis(psi_r, [0, 1], [c, t])
    return psi_r.reshape(-1)


def apply_SUM_inv(psi, c, t, N, d):
    """Inverse SUM: SUM†|c,t> = |c, t-c mod d>"""
    psi_r = psi.reshape([d] * N)
    psi_r = np.moveaxis(psi_r, [c, t], [0, 1])
    new_psi = np.zeros_like(psi_r)
    for a in range(d):
        for b in range(d):
            new_psi[a, b, ...] = psi_r[a, (b + a) % d, ...]
    psi_r = new_psi
    psi_r = np.moveaxis(psi_r, [0, 1], [c, t])
    return psi_r.reshape(-1)


# ==========================================
# 4. Random qudit Clifford circuit (brick-wall)
# ==========================================

def get_gate_dicts(d):
    F = fourier_matrix(d)
    P = phase_matrix(d)
    X = X_matrix(d)
    Z = Z_matrix(d)
    I = np.eye(d, dtype=complex)

    gates = {
        'F': F, 'Fd': F.conj().T,
        'P': P, 'Pd': P.conj().T,
        'X': X, 'Xd': np.linalg.matrix_power(X, d - 1),
        'Z': Z, 'Zd': np.linalg.matrix_power(Z, d - 1),
        'I': I,
    }
    inv = {
        'F': 'Fd', 'Fd': 'F',
        'P': 'Pd', 'Pd': 'P',
        'X': 'Xd', 'Xd': 'X',
        'Z': 'Zd', 'Zd': 'Z',
        'I': 'I',
    }
    return gates, inv


def random_qudit_clifford_circuit(N, d, depth=None):
    if depth is None:
        depth = N
    gates, _ = get_gate_dicts(d)
    single_choices = list(gates.keys())
    gate_list = []
    for layer in range(depth):
        # Single-qudit gates
        for i in range(N):
            gate_list.append(('single', i, np.random.choice(single_choices)))
        # SUM gates on even or odd pairs
        if layer % 2 == 0:
            for i in range(0, N - 1, 2):
                gate_list.append(('sum', i, i + 1))
        else:
            for i in range(1, N - 1, 2):
                gate_list.append(('sum', i, i + 1))
    return gate_list


def apply_circuit(psi, gate_list, N, d):
    gates, _ = get_gate_dicts(d)
    for g in gate_list:
        if g[0] == 'single':
            psi = apply_single_qudit(psi, gates[g[2]], g[1], N, d)
        else:
            psi = apply_SUM(psi, g[1], g[2], N, d)
    return psi


def apply_inverse_circuit(psi, gate_list, N, d):
    gates, inv = get_gate_dicts(d)
    for g in reversed(gate_list):
        if g[0] == 'single':
            psi = apply_single_qudit(psi, gates[inv[g[2]]], g[1], N, d)
        else:
            psi = apply_SUM_inv(psi, g[1], g[2], N, d)
    return psi


# ==========================================
# 5. Rotation and measurement
# ==========================================

def rotation_gate(d, theta_M, a=1):
    """
    Local non-Clifford rotation R^{(d)}_{X,a}(theta_M).
    Eq.(1): F_d^dagger diag(phases) F_d
    d=2: diag(1, e^{i theta_M pi/4})
    d=3: diag(1, e^{-i theta_M 2pi/9}, e^{i theta_M 2pi/9})
    d>3: diag(e^{-i theta_M (2pi/d) a k^3})
    """
    F = fourier_matrix(d)
    if d == 2:
        phases = np.array([1.0, np.exp(1j * theta_M * np.pi / 4)])
    elif d == 3:
        phases = np.array([1.0,
                           np.exp(-1j * theta_M * 2 * np.pi / 9),
                           np.exp(1j * theta_M * 2 * np.pi / 9)])
    else:
        phases = np.array([np.exp(-1j * theta_M * (2 * np.pi / d) * a * (k ** 3))
                           for k in range(d)])
    return F.conj().T @ np.diag(phases) @ F


def measure_qudit_single(N, psi, d):
    """Measure qudit 0 in computational basis"""
    psi_r = psi.reshape([d] * N)
    # Sum over all other qudits to get probabilities for qudit 0
    probs = np.sum(np.abs(psi_r) ** 2, axis=tuple(range(1, N)))
    probs = probs / np.sum(probs)  # normalize to fix floating point errors
    outcome = np.random.choice(d, p=probs)
    # Project
    new_psi = np.zeros_like(psi_r)
    new_psi[outcome] = psi_r[outcome]
    new_psi = new_psi / np.sqrt(probs[outcome])
    return new_psi.reshape(-1)


def run_trajectory(N, d, theta, T_max, seed=None):
    if seed is not None:
        np.random.seed(seed)

    psi = np.zeros(d ** N, dtype=complex)
    psi[0] = 1.0
    R = rotation_gate(d, theta)
    A_single = precompute_phase_point_operators(d, N)

    for t in range(T_max):
        gates = random_qudit_clifford_circuit(N, d, depth=4 * N)
        psi = apply_circuit(psi, gates, N, d)
        psi = apply_single_qudit(psi, R, 0, N, d)
        psi = measure_qudit_single(N, psi, d)
        psi = apply_inverse_circuit(psi, gates, N, d)

    return mana(psi, d, A_single)


def simulate_for_theta(args):
    N, d, theta, Nr, T_max = args
    A_single = precompute_phase_point_operators(d, N)
    results = []
    for r in range(Nr):
        seed = N * 100000 + d * 10000 + int(theta * 10000) + r
        m = run_trajectory(N, d, theta, T_max, seed)
        results.append(m)
    return np.mean(results), np.std(results) / np.sqrt(Nr)


# ==========================================
# 6. Main
# ==========================================

def main():
    configs = [
        {'d': 3, 'N': 2, 'Nr': 50, 'T_max': 500},
        {'d': 3, 'N': 3, 'Nr': 50, 'T_max': 500},
        {'d': 3, 'N': 4, 'Nr': 30, 'T_max': 500},
        {'d': 5, 'N': 2, 'Nr': 50, 'T_max': 500},
    ]

    theta_list = np.logspace(-4, 0, 10)  # Same as Fig.6: 1e-4 to 1
    all_results = {}

    for cfg in configs:
        d, N = cfg['d'], cfg['N']
        Nr, T_max = cfg['Nr'], cfg['T_max']
        key = f"d={d},N={N}"

        print(f"\n>>> {key} (Nr={Nr}, T_max={T_max})")
        m_mean = np.zeros(len(theta_list))
        m_err = np.zeros(len(theta_list))

        tasks = [(N, d, theta, Nr, T_max) for theta in theta_list]
        print(f"    Total tasks: {len(tasks)}, serial execution")
        t0 = time.time()

        results = [simulate_for_theta(task) for task in tasks]

        elapsed = time.time() - t0
        print(f"    {key} done, time: {elapsed:.1f}s ({elapsed / 60:.1f} min)")

        for th_idx, (theta, (m_m, m_e)) in enumerate(zip(theta_list, results)):
            m_mean[th_idx] = m_m
            m_err[th_idx] = m_e
            print(f"    theta={theta:.3e}: mana={m_m:.3e}, err={m_e:.3e}")

        all_results[key] = {
            'theta': theta_list.copy(),
            'mana': m_mean.copy(),
            'err': m_err.copy(),
            'd': d, 'N': N, 'Nr': Nr, 'T_max': T_max,
        }

        np.savez(f'mana_reproduction_{key.replace("=","").replace(",","_")}_partial.npz',
                 theta=theta_list, mana=m_mean, err=m_err, d=d, N=N, Nr=Nr, T_max=T_max)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {'d=3,N=2': '#e74c3c', 'd=3,N=3': '#2ecc71', 'd=5,N=2': '#3498db'}
    markers = {'d=3,N=2': 's', 'd=3,N=3': 'o', 'd=5,N=2': 'd'}

    for key in all_results:
        theta = all_results[key]['theta']
        mana_vals = all_results[key]['mana']
        err = all_results[key]['err']

        ax.loglog(theta, mana_vals, marker=markers.get(key, 'o'), color=colors.get(key, '#333'),
                  linestyle='-', linewidth=1.5, markersize=6,
                  label=f'{key} (Nr={all_results[key]["Nr"]})')
        ax.errorbar(theta, mana_vals, yerr=err, fmt='none',
                    color=colors.get(key, '#333'), alpha=0.5)

        # Fit linear slope in log-log (small angle)
        mask = theta < 0.1
        if np.sum(mask) > 2:
            coeffs = np.polyfit(np.log(theta[mask]), np.log(mana_vals[mask]), 1)
            print(f"{key}: fitted slope = {coeffs[0]:.3f}")

    # Reference line proportional to theta, aligned with data
    ref_th = np.array([1e-4, 3e-2])
    ref_mana = ref_th * 2.5
    ax.loglog(ref_th, ref_mana, 'k--', linewidth=1.5, label=r'$\propto \theta$')

    ax.set_xlabel(r'$\theta$', fontsize=14)
    ax.set_ylabel(r'$\bar{\mathcal{M}}^{\mathrm{SS}}$', fontsize=14)
    ax.set_title('Steady-state Mana vs rotation angle (qudit monitored circuit)', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, which='both', ls='-', alpha=0.3)
    ax.set_xlim([1e-4, 1])
    ax.set_ylim([5e-5, 1])

    plt.tight_layout()
    plt.savefig('mana_reproduction.png', dpi=300)
    print("\nPlot saved: mana_reproduction.png")

    np.savez('mana_reproduction_data.npz', all_results=all_results, theta_list=theta_list)
    print("Data saved: mana_reproduction_data.npz")


if __name__ == '__main__':
    main()
