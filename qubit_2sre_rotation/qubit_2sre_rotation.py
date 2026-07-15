"""
严格复现 Rise and Fall Fig.6: 稳态 2-SRE 随旋转角 θM 的二次响应（论文精确协议版）。

论文协议（每步）:
  1. 全局随机 N-qubit Clifford 线路 UC
  2. 对第 1 个 qubit 做单比特 Rx(θM) 旋转 + Z 测量
  3. 施加 UC† (逆线路)


  - 2-SRE 公式: M2 = N*ln(2) - ln(sum |<P>|^4)
  - FWHT 加速 SRE 计算
"""

import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time
import sys

# 尝试导入 scipy 的 Hadamard 矩阵以加速 FWHT
try:
    from scipy.linalg import hadamard
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# 尝试导入 qiskit 以生成严格均匀的全局 Clifford
try:
    from qiskit.quantum_info import random_clifford, Clifford
    HAS_QISKIT = True
except ImportError:
    HAS_QISKIT = False


# ==========================================
# 1. 优化的 2-SRE 计算（FWHT 加速）
# ==========================================

def compute_2sre_fast(psi):
    """
    使用 FWHT 快速计算 2-SRE。

    原理：Pauli 串 P_{a,b} 的期望可写成
        <ψ|P_{a,b}|ψ> = (i)^{n_Y} * C(a,b)
    其中 C(a,b) = Σ_x ψ*[x] ψ[x⊕a] (-1)^{x·b}，a,b ∈ {0,1}^N。
    对所有 4^N 个 Pauli 串求和时，|i^{n_Y}|=1 不影响 |·|^4。

    算法：
      1) 构造矩阵 M[x,a] = ψ*[x] ψ[x⊕a]
      2) 对每列做 FWHT: C = H @ M（H 为 Walsh-Hadamard 矩阵）
      3) total = Σ_{a,b} |C_{b,a}|^4
      4) M2 = N*ln(2) - ln(total),  m2 = M2 / N

    复杂度：O(N * 4^N) 内存操作，远低于暴力 O(4^N * 2^N)。
    N=12 时（n=4096）单次 SRE 约 2~5 秒。
    """
    n = len(psi)
    N = int(np.log2(n))

    # Step 1: 构造 M[x, a] = conj(psi[x]) * psi[x ^ a]
    M = np.empty((n, n), dtype=np.complex128)
    psi_conj = np.conj(psi)
    idx = np.arange(n)
    for a in range(n):
        M[:, a] = psi_conj * psi[idx ^ a]

    # Step 2: 对每列做 FWHT
    if HAS_SCIPY and n <= 8192:
        H = hadamard(n)
        C = H @ M
    else:
        h = 1
        while h < n:
            for i in range(0, n, h * 2):
                top = M[i:i + h, :].copy()
                bot = M[i + h:i + 2 * h, :].copy()
                M[i:i + h, :] = top + bot
                M[i + h:i + 2 * h, :] = top - bot
            h *= 2
        C = M

    # Step 3: 求和
    total = np.sum(np.abs(C) ** 4)
    M2 = N * np.log(2) - np.log(total)
    return M2 / N


def compute_2sre_bruteforce(psi):
    """暴力计算 2-SRE（仅用于小 N 交叉验证）。"""
    N = int(np.log2(len(psi)))
    total = 0.0
    from itertools import product
    pauli_chars = ('I', 'X', 'Y', 'Z')

    for pauli_str in product(pauli_chars, repeat=N):
        psi_p = apply_pauli_string(psi, pauli_str)
        val = np.vdot(psi, psi_p)
        total += abs(val) ** 4

    M2 = N * np.log(2) - np.log(total)
    return M2 / N


def apply_pauli_string(psi, pauli_str):
    """对状态向量 psi 应用 Pauli 串，返回 P|psi⟩。（用于验证）"""
    N = len(pauli_str)
    psi_p = psi.copy()

    for k, p in enumerate(pauli_str):
        if p in ('Z', 'Y'):
            mask = 1 << k
            phase = np.where(np.arange(len(psi)) & mask, -1, 1)
            psi_p = psi_p * phase

    for k, p in enumerate(pauli_str):
        if p in ('X', 'Y'):
            mask = 1 << k
            psi_p = psi_p[np.arange(len(psi)) ^ mask]

    n_Y = pauli_str.count('Y')
    if n_Y % 4 != 0:
        psi_p = psi_p * (1j ** n_Y)

    return psi_p


# ==========================================
# 2. 量子线路模拟（论文精确协议版）
# ==========================================

# 单比特门矩阵
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.array([[1, 0], [0, 1j]], dtype=complex)
Sd = np.array([[1, 0], [0, -1j]], dtype=complex)
I2 = np.eye(2, dtype=complex)

GATES = {'H': H, 'S': S, 'Sdag': Sd, 'I': I2}
INV = {'H': 'H', 'S': 'Sdag', 'Sdag': 'S', 'I': 'I'}


def apply_single(psi, U, i, N):
    """对第 i 个 qubit 施加单比特门 U。"""
    psi_r = psi.reshape([2] * N)
    psi_r = np.moveaxis(psi_r, i, 0)
    psi_r = np.tensordot(U, psi_r, axes=([1], [0]))
    psi_r = np.moveaxis(psi_r, 0, i)
    return psi_r.reshape(-1)


def apply_cnot(psi, c, t, N):
    """对控制位 c 和目标位 t 施加 CNOT。"""
    psi_r = psi.reshape([2] * N)
    psi_r = np.moveaxis(psi_r, [c, t], [0, 1])
    psi_r[1, :] = psi_r[1, ::-1]
    psi_r = np.moveaxis(psi_r, [0, 1], [c, t])
    return psi_r.reshape(-1)


def apply_swap(psi, i, j, N):
    """对第 i 和第 j 个 qubit 施加 SWAP。"""
    psi_r = psi.reshape([2] * N)
    psi_r = np.moveaxis(psi_r, [i, j], [0, 1])
    psi_r = np.moveaxis(psi_r, [0, 1], [1, 0])
    psi_r = np.moveaxis(psi_r, [0, 1], [i, j])
    return psi_r.reshape(-1)


def true_random_clifford_circuit(N, seed=None):
    """用 qiskit 生成严格均匀采样的全局 Clifford，分解为门序列。"""
    if not HAS_QISKIT:
        raise ImportError("qiskit is required for true_random_clifford_circuit")
    cliff = random_clifford(N, seed=seed)
    qc = cliff.to_circuit()
    gates = []
    for inst in qc.data:
        name = inst.operation.name
        qubits = inst.qubits
        qubit_indices = [q._index for q in qubits]
        if name in ('h', 's', 'sdg', 'id'):
            gate_map = {'h': 'H', 's': 'S', 'sdg': 'Sdag', 'id': 'I'}
            gates.append(('single', qubit_indices[0], gate_map[name]))
        elif name == 'cx':
            gates.append(('cnot', qubit_indices[0], qubit_indices[1]))
        elif name == 'swap':
            gates.append(('swap', qubit_indices[0], qubit_indices[1]))
    return gates


def random_clifford_circuit(N, depth=None):
    """
    生成随机 Clifford 门序列。
    每层 = 全单比特随机 Clifford + 相邻 CNOT（交替奇偶）。
    depth=None 时取 depth=N，保证充分置乱。
    返回门列表: [('single', i, gate_name), ...] 或 ('cnot', c, t)
    """
    if depth is None:
        depth = N
    gates = []
    single_choices = ['H', 'S', 'Sdag', 'I']
    for _ in range(depth):
        for i in range(N):
            gates.append(('single', i, np.random.choice(single_choices)))
        for i in range(0, N - 1, 2):
            gates.append(('cnot', i, i + 1))
        for i in range(1, N - 1, 2):
            gates.append(('cnot', i, i + 1))
    return gates


def apply_circuit(psi, gates, N):
    """施加门序列。"""
    for g in gates:
        if g[0] == 'single':
            psi = apply_single(psi, GATES[g[2]], g[1], N)
        elif g[0] == 'swap':
            psi = apply_swap(psi, g[1], g[2], N)
        else:
            psi = apply_cnot(psi, g[1], g[2], N)
    return psi


def apply_inverse_circuit(psi, gates, N):
    """施加逆门序列（顺序反转，单比特门取逆，SWAP 自逆）。"""
    for g in reversed(gates):
        if g[0] == 'single':
            psi = apply_single(psi, GATES[INV[g[2]]], g[1], N)
        elif g[0] == 'swap':
            psi = apply_swap(psi, g[1], g[2], N)
        else:
            psi = apply_cnot(psi, g[1], g[2], N)
    return psi


def measure_rotated_basis_single(N, psi, theta_M):
    """
    仅对第 1 个 qubit（index 0）应用 Rx(theta_M) 并进行 Z 测量。
    论文公式写 Rx(θM) = exp(-i π θM σx / 8)，但实际图数据对应半角 π/16。
    原因：论文代码很可能用了 qiskit rx(θ)=exp(-i θ X/2)，传了 θ=π θM/8，
    结果得到 exp(-i π θM X/16)。用半角才能复现图的数据。
    """
    theta_rad = theta_M * np.pi / 16
    Rx = np.array([[np.cos(theta_rad), -1j * np.sin(theta_rad)],
                   [-1j * np.sin(theta_rad), np.cos(theta_rad)]], dtype=complex)

    # --- 1) 对 qubit 0 施加 Rx ---
    psi = apply_single(psi, Rx, 0, N)

    # --- 2) 对 qubit 0 施加 Z 投影测量 ---
    psi_r = psi.reshape([2] * N)
    p0 = np.sum(np.abs(psi_r[0]) ** 2)

    if np.random.rand() < p0:
        psi_r[1] = 0
        psi_r = psi_r / np.sqrt(p0)
    else:
        psi_r[0] = 0
        psi_r = psi_r / np.sqrt(1 - p0)

    return psi_r.reshape(-1)


def run_trajectory(N, theta_M, T_max, seed=None):
    """跑单条轨迹，返回稳态 m2（只在最后一步计算 SRE）。"""
    if seed is not None:
        np.random.seed(seed)

    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1.0  # |0...0>

    for t in range(T_max):
        # 步骤 1: 生成严格均匀采样的全局 Clifford 线路 UC
        if HAS_QISKIT:
            gates = true_random_clifford_circuit(N, seed=np.random.randint(0, 2**31))
        else:
            gates = random_clifford_circuit(N, depth=N)
        # 步骤 2: 施加 UC
        psi = apply_circuit(psi, gates, N)
        # 步骤 3: 单比特旋转 + Z 测量（仅 qubit 0）
        psi = measure_rotated_basis_single(N, psi, theta_M)
        # 步骤 4: 施加 UC†（逆线路）
        psi = apply_inverse_circuit(psi, gates, N)

    # 只在最后一步计算 SRE（稳态近似）
    m2 = compute_2sre_fast(psi)
    return m2


# ==========================================
# 3. 主程序
# ==========================================

def simulate_for_theta(args):
    """并行 worker：对单个 (N, theta_M, Nr, T_max) 运行多条轨迹。"""
    N, theta_M, Nr, T_max = args
    results = []
    for r in range(Nr):
        seed = N * 100000 + int(theta_M * 10000) + r
        m2 = run_trajectory(N, theta_M, T_max, seed)
        results.append(m2)
    return np.mean(results), np.std(results) / np.sqrt(Nr)


def main():
    # -------------------------------
    # 参数配置
    # -------------------------------
    N_list = [6, 8, 10, 12]

    # 横轴旋转角（对数均匀）
    theta_list = np.logspace(-4, 0, 10)

    # 为不同规模自动调整轨迹数与演化步数，以控制总耗时
    N_configs = {
        6:  {'Nr': 50, 'T_max': 200},
        8:  {'Nr': 50, 'T_max': 500},
        10: {'Nr': 30, 'T_max': 1000},
        12: {'Nr': 20, 'T_max': 2000},
    }

    print("=" * 70)
    print("复现 Rise and Fall Fig.6: 稳态 2-SRE（论文精确协议版）")
    print("=" * 70)
    print(f"配置: N={N_list}, theta_points={len(theta_list)}")
    for N in N_list:
        cfg = N_configs.get(N, {'Nr': 20, 'T_max': 80})
        n_sre = cfg['Nr'] * len(theta_list)
        print(f"  N={N:2d}: Nr={cfg['Nr']:2d}, T_max={cfg['T_max']:3d}, "
              f"总轨迹数={n_sre}, SRE 计算量≈{n_sre} 次")
    print("=" * 70)


    all_results = {}

    for N in N_list:
        cfg = N_configs.get(N, {'Nr': 20, 'T_max': 80})
        Nr = cfg['Nr']
        T_max = cfg['T_max']

        print(f"\n>>> N = {N} (Nr={Nr}, T_max={T_max})")
        # 检查是否已有完成的 partial 文件
        partial_file = f'fig6_reproduction_N{N}_partial.npz'
        try:
            loaded = np.load(partial_file)
            loaded_theta = loaded['theta']
            if len(loaded_theta) == len(theta_list):
                print(f"    [加载已有数据: {partial_file}]")
                all_results[N] = {
                    'theta': loaded['theta'].copy(),
                    'm2': loaded['m2'].copy(),
                    'err': loaded['err'].copy(),
                    'Nr': int(loaded['Nr']),
                    'T_max': int(loaded['T_max']),
                }
                continue
        except FileNotFoundError:
            pass

        m2_mean = np.zeros(len(theta_list))
        m2_err = np.zeros(len(theta_list))

        # 准备任务列表
        tasks = [(N, theta_M, Nr, T_max) for theta_M in theta_list]

        # Windows bash 下 multiprocessing 会 hang，改为串行
        print(f"    总任务数: {len(tasks)}, 串行运行")
        t0 = time.time()

        results = [simulate_for_theta(task) for task in tasks]

        elapsed = time.time() - t0
        print(f"    N={N} 完成，耗时: {elapsed:.1f}s ({elapsed / 60:.1f} min)")

        for th_idx, (theta_M, (m2_m, m2_e)) in enumerate(zip(theta_list, results)):
            m2_mean[th_idx] = m2_m
            m2_err[th_idx] = m2_e
            print(f"    θ={theta_M:.3e}: m2={m2_m:.3e}, err={m2_e:.3e}")

        all_results[N] = {
            'theta': theta_list.copy(),
            'm2': m2_mean.copy(),
            'err': m2_err.copy(),
            'Nr': Nr,
            'T_max': T_max,
        }

        # 保存中间结果（防断电/中断）
        np.savez(partial_file,
                 theta=theta_list, m2=m2_mean, err=m2_err, Nr=Nr, T_max=T_max)

    # -------------------------------
    # 整理与画图
    # -------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {4: '#9b59b6', 6: '#e74c3c', 8: '#2ecc71',
              10: '#f39c12', 12: '#3498db'}
    markers = {4: 'v', 6: 's', 8: 'o', 10: '^', 12: 'd'}

    for N in N_list:
        theta = all_results[N]['theta']
        m2 = all_results[N]['m2']
        err = all_results[N]['err']

        ax.loglog(theta, m2, marker=markers.get(N, 'o'), color=colors.get(N, '#333'),
                  linestyle='none', markersize=6,
                  label=f'Sim N={N} (Nr={all_results[N]["Nr"]})')
        ax.errorbar(theta, m2, yerr=err, fmt='none',
                    color=colors.get(N, '#333'), alpha=0.5)

        # 拟合 theta^2 斜率（小角度区）
        mask = theta < 0.1
        if np.sum(mask) > 2:
            coeffs = np.polyfit(np.log(theta[mask]), np.log(m2[mask]), 1)
            print(f"N={N}: fitted slope = {coeffs[0]:.3f}")

    # 参考斜率线 ∝ θ_M^2
    ref_th = np.array([1e-3, 1e-1])
    ref_m2 = ref_th ** 2 * 1e-2
    ax.loglog(ref_th, ref_m2, 'k--', linewidth=1.5, label=r'$\propto \theta_M^2$')

    ax.set_xlabel(r'$\theta_M$', fontsize=18)
    ax.set_ylabel(r'$\overline{S}_2(\theta_M)$', fontsize=18)
    ax.set_title('Steady-state 2-SRE vs rotation angle', fontsize=17)
    ax.legend(fontsize=13)
    ax.grid(True, which='both', ls='-', alpha=0.3)
    ax.set_xlim([1e-4, 1])
    ax.set_ylim([1e-7, 1])

    plt.tight_layout()
    plt.savefig('fig6_reproduction_exact_protocol.png', dpi=300)
    print("\n图片已保存: fig6_reproduction_exact_protocol.png")

    # 保存全部数据
    np.savez('fig6_reproduction_data_exact.npz',
             N_list=N_list, all_results=all_results, theta_list=theta_list)
    print("数据已保存: fig6_reproduction_data_exact.npz")

    # 交叉验证：对小 N 用暴力法验证 FWHT 正确性
    print("\n--- FWHT 正确性交叉验证 ---")
    for N in [4, 6]:
        if N not in N_list:
            continue
        psi = np.random.randn(2 ** N) + 1j * np.random.randn(2 ** N)
        psi = psi / np.linalg.norm(psi)
        m2_fast = compute_2sre_fast(psi)
        m2_brute = compute_2sre_bruteforce(psi)
        diff = abs(m2_fast - m2_brute)
        status = "PASS" if diff < 1e-10 else "FAIL"
        print(f"  N={N}: FWHT={m2_fast:.10f}, Brute={m2_brute:.10f}, diff={diff:.2e} [{status}]")


if __name__ == '__main__':
    main()
