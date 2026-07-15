"""
Load existing Fig.6 data and replot with lines connecting the points.
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.load('fig6_reproduction_data_exact.npz', allow_pickle=True)
N_list = data['N_list']
theta_list = data['theta_list']
all_results = data['all_results'].item()

fig, ax = plt.subplots(figsize=(7, 7))
colors = {4: '#9b59b6', 6: '#e74c3c', 8: '#2ecc71',
          10: '#f39c12', 12: '#3498db'}
markers = {4: 'v', 6: 's', 8: 'o', 10: '^', 12: 'd'}

for N in N_list:
    theta = all_results[N]['theta']
    m2 = all_results[N]['m2']
    err = all_results[N]['err']

    # Plot with both markers and connecting lines
    ax.loglog(theta, m2, marker=markers.get(N, 'o'), color=colors.get(N, '#333'),
              linestyle='-', linewidth=1.5, markersize=6,
              label=f'N={N}')
    ax.errorbar(theta, m2, yerr=err, fmt='none',
                color=colors.get(N, '#333'), alpha=0.5)

    # Fit slope for small angles
    mask = theta < 0.1
    if np.sum(mask) > 2:
        coeffs = np.polyfit(np.log(theta[mask]), np.log(m2[mask]), 1)
        print(f"N={N}: fitted slope = {coeffs[0]:.3f}")

# Reference line ~ theta^2
ref_th = np.array([1e-3, 1e-1])
ref_m2 = ref_th ** 2 * 1e-2
ax.loglog(ref_th, ref_m2, 'k--', linewidth=1.5, label=r'$\propto \theta_M^2$')

ax.set_xlabel(r'$\theta_M$', fontsize=18)
ax.set_ylabel(r'$\bar{m}_2^{\mathrm{SS}}$', fontsize=18)
ax.set_title('Reproduction of Fig.6 (exact paper protocol)', fontsize=17)
ax.legend(fontsize=13)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.grid(True, which='both', ls='-', alpha=0.3)
ax.set_xlim([1e-4, 1])
ax.set_ylim([1e-9, 1])

plt.tight_layout()
plt.savefig('fig6_reproduction_exact_protocol.png', dpi=300)
print("\n图片已保存: fig6_reproduction_exact_protocol.png")
