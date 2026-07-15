"""
Qudit monitored circuit steady-state Mana vs rotation angle.
Protocol same as Fig.6 but for qudits (odd prime d).
Strict global Clifford: uniform sampling over Sp(2N, F_d).
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from itertools import product
from test_strict_clifford import sample_symplectic, omega as symplectic_omega

# ==========================================
# 1. Basic qudit operators
# ==========================================

def omega_phase(d, k):
    return np.exp(2j * np.pi * k / d)


def X_matrix(d):
    X = np.zeros((d, d), dtype=complex)
    for j in range(d):
        X[(j + 1) % d, j] = 1.0
    return X


def Z_matrix(d):
    Z = np.zeros((d, d), dtype=complex)
    for j in range(d):
        Z[j, j] = omega_phase(d, j)
    return Z


def fourier_matrix(d):
    F = np.zeros((d, d), dtype=complex)
    for j in range(d):
        for k in range(d):
            F[k, j] = omega_phase(d, j * k) / np.sqrt(d)
    return F


def D_matrix(d, p, q):
    """Single-qudit displacement operator w(p,q) = tau^{pq} X^p Z^q."""
    X = X_matrix(d)
    Z = Z_matrix(d)
    Xp = np.linalg.matrix_power(X, p)
    Zq = np.linalg.matrix_power(Z, q)
    tau = -np.exp(1j * np.pi / d)
    return Xp @ Zq * (tau ** (p * q))


def multi_qudit_displacement(d, a, b, N):
    """N-qudit displacement: w(a,b) = \otimes_i w(a_i, b_i)."""
    op = np.eye(1, dtype=complex)
    for i in range(N):
        op = np.kron(op, D_matrix(d, int(a[i]), int(b[i])))
    return op


def index_to_digits(idx, N, d):
    digits = []
    for i in range(N):
        digits.append(idx % d)
        idx //= d
    return digits


def digits_to_index(digits, d):
    idx = 0
    for i in reversed(range(len(digits))):
        idx = idx * d + digits[i]
    return idx


# ==========================================
# 2. Stabilizer state from symplectic matrix
# ==========================================

def project_stabilizer_state(stabilizers, d):
    """
    Find joint +1 eigenstate of commuting stabilizers.
    Uses projector (1/d) sum_{j=0}^{d-1} P^j for each stabilizer P.
    """
    D = stabilizers[0].shape[0]
    rho = np.eye(D, dtype=complex)
    for P in stabilizers:
        proj = sum(np.linalg.matrix_power(P, j) for j in range(d)) / d
        rho = proj @ rho
    # Hermitianize to fix numerical drift
    rho = (rho + rho.conj().T) / 2
    eigvals, eigvecs = np.linalg.eigh(rho)
    idx = np.argmax(eigvals.real)
    psi = eigvecs[:, idx]
    return psi / np.linalg.norm(psi)


def build_clifford_basis(S, N, d):
    """
    From symplectic matrix S (columns [e_1..e_N, f_1..f_N]),
    build the stabilizer state |psi_0> = U_C |0...0> and
    the basis B_j = Q^j |psi_0> where Q_i corresponds to e_i (image of X_i).

    Returns:
        psi_0: U_C |0...0>
        basis_states: list of d^N vectors, basis_states[j] = Q^j |psi_0>
    """
    D = d ** N
    e_cols = [S[:, i] for i in range(N)]      # images of X_i
    f_cols = [S[:, N + i] for i in range(N)]  # images of Z_i

    # Stabilizers of |psi_0> are images of Z_i
    stabilizers = []
    for f in f_cols:
        a = f[:N]
        b = f[N:]
        P = multi_qudit_displacement(d, a, b, N)
        stabilizers.append(P)

    psi_0 = project_stabilizer_state(stabilizers, d)

    # Build basis B_j = Q_1^{j_1} ... Q_N^{j_N} |psi_0>
    basis_states = []
    for j_idx in range(D):
        j_digits = index_to_digits(j_idx, N, d)
        state = psi_0.copy()
        for i in range(N):
            if j_digits[i] != 0:
                a = e_cols[i][:N]
                b = e_cols[i][N:]
                Q = multi_qudit_displacement(d, a, b, N)
                state = np.linalg.matrix_power(Q, j_digits[i]) @ state
        basis_states.append(state)

    return psi_0, basis_states


def apply_UC(state_comp, basis_states):
    """Apply U_C to a state given in computational basis."""
    result = np.zeros_like(basis_states[0])
    for j, coeff in enumerate(state_comp):
        if abs(coeff) > 0:
            result += coeff * basis_states[j]
    return result


def apply_UC_dag(state_cliff, basis_states):
    """Apply U_C^\dagger to a state, returning coefficients in computational basis."""
    return np.array([np.vdot(bj, state_cliff) for bj in basis_states])


# ==========================================
# 3. Mana computation
# ==========================================

def precompute_phase_point_operators(d, n):
    A_single = {}
    for u in range(d):
        for v in range(d):
            A = np.zeros((d, d), dtype=complex)
            for p in range(d):
                for q in range(d):
                    A += omega_phase(d, u * q - v * p) * D_matrix(d, p, q)
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
# 4. Rotation and measurement
# ==========================================

def rotation_gate(d, theta_M, a=1):
    """
    Local non-Clifford rotation R^{(d)}_{X,a}(theta_M).
    Eq.(1): F_d^\dagger diag(phases) F_d
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


def apply_single_qudit(psi, U, i, N, d):
    psi_r = psi.reshape([d] * N)
    psi_r = np.moveaxis(psi_r, i, 0)
    psi_r = np.tensordot(U, psi_r, axes=([1], [0]))
    psi_r = np.moveaxis(psi_r, 0, i)
    return psi_r.reshape(-1)


def measure_qudit_single(N, psi, d):
    """Measure qudit 0 in computational basis."""
    psi_r = psi.reshape([d] * N)
    probs = np.sum(np.abs(psi_r) ** 2, axis=tuple(range(1, N)))
    probs = probs / np.sum(probs)
    outcome = np.random.choice(d, p=probs)
    new_psi = np.zeros_like(psi_r)
    new_psi[outcome] = psi_r[outcome]
    new_psi = new_psi / np.sqrt(probs[outcome])
    return new_psi.reshape(-1)


# ==========================================
# 5. Trajectory with strict Clifford
# ==========================================

def run_trajectory_strict(N, d, theta, T_max, seed=None):
    if seed is not None:
        np.random.seed(seed)

    # Sample random symplectic matrix
    S = sample_symplectic(N, d)

    # Build U_C basis
    psi_0, basis_states = build_clifford_basis(S, N, d)

    # Verify orthonormality (debug only, can remove later)
    # D = d ** N
    # for i in range(min(D, 5)):
    #     for j in range(min(D, 5)):
    #         ov = np.vdot(basis_states[i], basis_states[j])
    #         expected = 1.0 if i == j else 0.0
    #         if abs(ov - expected) > 1e-10:
    #             print(f"Orthonormality fail: <{i}|{j}> = {ov}")

    R = rotation_gate(d, theta)
    A_single = precompute_phase_point_operators(d, N)

    # Start in |0...0>
    state = np.zeros(d ** N, dtype=complex)
    state[0] = 1.0

    for t in range(T_max):
        # Apply U_C
        state = apply_UC(state, basis_states)
        # Apply rotation on qudit 0
        state = apply_single_qudit(state, R, 0, N, d)
        # Measure qudit 0
        state = measure_qudit_single(N, state, d)
        # Apply U_C^\dagger
        state = apply_UC_dag(state, basis_states)

    return mana(state, d, A_single)


def simulate_for_theta(args):
    N, d, theta, Nr, T_max = args
    A_single = precompute_phase_point_operators(d, N)
    results = []
    for r in range(Nr):
        seed = N * 100000 + d * 10000 + int(theta * 10000) + r
        m = run_trajectory_strict(N, d, theta, T_max, seed)
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

    theta_list = np.logspace(-4, 0, 10)
    all_results = {}

    for cfg in configs:
        d, N = cfg['d'], cfg['N']
        Nr, T_max = cfg['Nr'], cfg['T_max']
        key = f"d={d},N={N}"

        print(f"\n>>> {key} (Nr={Nr}, T_max={T_max}, STRICT CLIFFORD)")
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

        np.savez(f'mana_strict_{key.replace("=","").replace(",","_")}_partial.npz',
                 theta=theta_list, mana=m_mean, err=m_err, d=d, N=N, Nr=Nr, T_max=T_max)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {'d=3,N=2': '#e74c3c', 'd=3,N=3': '#2ecc71',
              'd=3,N=4': '#9b59b6', 'd=5,N=2': '#3498db'}
    markers = {'d=3,N=2': 's', 'd=3,N=3': 'o',
               'd=3,N=4': '^', 'd=5,N=2': 'd'}

    for key in all_results:
        theta = all_results[key]['theta']
        mana_vals = all_results[key]['mana']
        err = all_results[key]['err']

        ax.loglog(theta, mana_vals, marker=markers.get(key, 'o'), color=colors.get(key, '#333'),
                  linestyle='-', linewidth=1.5, markersize=6,
                  label=f'{key} (Nr={all_results[key]["Nr"]})')
        ax.errorbar(theta, mana_vals, yerr=err, fmt='none',
                    color=colors.get(key, '#333'), alpha=0.5)

        mask = theta < 0.1
        if np.sum(mask) > 2:
            coeffs = np.polyfit(np.log(theta[mask]), np.log(mana_vals[mask]), 1)
            print(f"{key}: fitted slope = {coeffs[0]:.3f}")

    ref_th = np.array([1e-4, 3e-2])
    ref_mana = ref_th * 2.5
    ax.loglog(ref_th, ref_mana, 'k--', linewidth=1.5, label=r'$\propto \theta_M$')

    ax.set_xlabel(r'$\theta_M$', fontsize=18)
    ax.set_ylabel(r'$\overline{\mathcal{M}}(\theta_M)$', fontsize=18)
    ax.set_title('Steady-state Mana vs rotation angle', fontsize=17)
    ax.legend(fontsize=13)
    ax.grid(True, which='both', ls='-', alpha=0.3)
    ax.set_xlim([1e-4, 1])
    ax.set_ylim([5e-5, 2])

    plt.tight_layout()
    plt.savefig('mana_strict_clifford.png', dpi=300)
    print("\nPlot saved: mana_strict_clifford.png")

    np.savez('mana_strict_data.npz', all_results=all_results, theta_list=theta_list)
    print("Data saved: mana_strict_data.npz")


if __name__ == '__main__':
    main()
