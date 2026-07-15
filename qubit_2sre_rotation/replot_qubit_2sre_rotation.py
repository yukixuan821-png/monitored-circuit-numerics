"""
Replot fig6 (2-SRE) and mana_strict (Mana) with updated notation:
  - 2-SRE: ylabel = $\overline{S}_2(\theta_M)$
  - Mana : ylabel = $\overline{\mathcal{M}}(\theta_M)$
Reads existing data (fig6_data.txt for fig6, mana_strict_data.npz for mana).
"""

import numpy as np
import matplotlib.pyplot as plt


def read_fig6_table(path='fig6_data.txt'):
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith('-') or s.startswith('theta_M'):
                continue
            parts = s.split()
            if len(parts) < 9:
                continue
            try:
                vals = [float(x) for x in parts[:9]]
            except ValueError:
                continue
            rows.append(vals)
    arr = np.array(rows)
    theta = arr[:, 0]
    data = {
        6:  {'m2': arr[:, 1], 'err': arr[:, 2], 'Nr': 50,  'T_max': 200},
        8:  {'m2': arr[:, 3], 'err': arr[:, 4], 'Nr': 50,  'T_max': 500},
        10: {'m2': arr[:, 5], 'err': arr[:, 6], 'Nr': 30,  'T_max': 1000},
        12: {'m2': arr[:, 7], 'err': arr[:, 8], 'Nr': 20,  'T_max': 2000},
    }
    for N in data:
        data[N]['theta'] = theta
    return theta, data


def plot_fig6():
    theta, data = read_fig6_table()
    N_list = [6, 8, 10, 12]
    colors = {6: '#e74c3c', 8: '#2ecc71', 10: '#f39c12', 12: '#3498db'}
    markers = {6: 's', 8: 'o', 10: '^', 12: 'd'}

    fig, ax = plt.subplots(figsize=(8, 8))
    for N in N_list:
        m2 = data[N]['m2']
        err = data[N]['err']
        ax.loglog(theta, m2, marker=markers[N], color=colors[N],
                  linestyle='-', linewidth=1.5, markersize=6,
                  label=f'N={N}')
        ax.errorbar(theta, m2, yerr=err, fmt='none',
                    color=colors[N], alpha=0.5)
        mask = theta < 0.1
        if np.sum(mask) > 2:
            c = np.polyfit(np.log(theta[mask]), np.log(m2[mask]), 1)
            print(f"N={N}: fitted slope = {c[0]:.3f}")

    ref_th = np.array([1e-3, 1e-1])
    ref_m2 = ref_th ** 2 * 1e-2
    ax.loglog(ref_th, ref_m2, 'k--', linewidth=1.5,
              label=r'$\propto \theta_M^2$')

    ax.set_xlabel(r'$\theta_M$', fontsize=18)
    ax.set_ylabel(r'$\overline{S}_2(\theta_M)$', fontsize=18)
    ax.set_title('Steady-state 2-SRE vs rotation angle', fontsize=17)
    ax.legend(fontsize=13)
    ax.grid(True, which='both', ls='-', alpha=0.3)
    ax.set_xlim([1e-4, 1])
    ax.set_ylim([1e-9, 1])
    plt.tight_layout()
    plt.savefig('fig6_reproduction_exact_protocol.png', dpi=300)
    print("Saved: fig6_reproduction_exact_protocol.png")


def plot_mana():
    d = np.load('mana_strict_data.npz', allow_pickle=True)
    all_results = d['all_results'].item()

    colors = {'d=3,N=2': '#e74c3c', 'd=3,N=3': '#2ecc71',
              'd=3,N=4': '#9b59b6', 'd=5,N=2': '#3498db'}
    markers = {'d=3,N=2': 's', 'd=3,N=3': 'o',
               'd=3,N=4': '^', 'd=5,N=2': 'd'}

    fig, ax = plt.subplots(figsize=(8, 6))
    for key in all_results:
        theta = all_results[key]['theta']
        mana_vals = all_results[key]['mana']
        err = all_results[key]['err']
        ax.loglog(theta, mana_vals, marker=markers.get(key, 'o'),
                  color=colors.get(key, '#333'),
                  linestyle='-', linewidth=1.5, markersize=6,
                  label=f'{key} (Nr={all_results[key]["Nr"]})')
        ax.errorbar(theta, mana_vals, yerr=err, fmt='none',
                    color=colors.get(key, '#333'), alpha=0.5)
        mask = theta < 0.1
        if np.sum(mask) > 2:
            c = np.polyfit(np.log(theta[mask]), np.log(mana_vals[mask]), 1)
            print(f"{key}: fitted slope = {c[0]:.3f}")

    ref_th = np.array([1e-4, 3e-2])
    ref_mana = ref_th * 2.5
    ax.loglog(ref_th, ref_mana, 'k--', linewidth=1.5,
              label=r'$\propto \theta_M$')

    ax.set_xlabel(r'$\theta_M$', fontsize=18)
    ax.set_ylabel(r'$\overline{\mathcal{M}}(\theta_M)$', fontsize=18)
    ax.set_title('Steady-state Mana vs rotation angle', fontsize=17)
    ax.legend(fontsize=13)
    ax.grid(True, which='both', ls='-', alpha=0.3)
    ax.set_xlim([1e-4, 1])
    ax.set_ylim([5e-5, 2])
    plt.tight_layout()
    plt.savefig('mana_strict_clifford.png', dpi=300)
    print("Saved: mana_strict_clifford.png")


if __name__ == '__main__':
    print("--- fig6 (2-SRE) ---")
    plot_fig6()
    print("\n--- mana (strict Clifford) ---")
    plot_mana()
