"""
Quick test: d=5,N=2 with smaller theta range to check linear slope.
"""

import numpy as np
from reproduce_mana import run_trajectory

def quick_test():
    d, N = 5, 2
    Nr = 50
    T_max = 500
    theta_list = np.logspace(-4, -3, 8)

    print(f"d={d}, N={N}, Nr={Nr}, T_max={T_max}, depth=4N={4*N}")
    print(f"theta range: {theta_list[0]:.2e} to {theta_list[-1]:.2e}")
    print("-" * 50)

    means = []
    errs = []
    for theta_M in theta_list:
        results = []
        for r in range(Nr):
            seed = N * 100000 + d * 10000 + int(theta_M * 100000) + r
            m = run_trajectory(N, d, theta_M, T_max, seed)
            results.append(m)
        mean = np.mean(results)
        err = np.std(results) / np.sqrt(Nr)
        means.append(mean)
        errs.append(err)
        print(f"theta_M={theta_M:.3e}: mana={mean:.3e} +/- {err:.3e}")

    # Fit slope in log-log
    coeffs = np.polyfit(np.log(theta_list), np.log(np.array(means)), 1)
    print(f"\nFitted slope (log-log): {coeffs[0]:.3f}")

    # Local slopes
    log_t = np.log(theta_list)
    log_m = np.log(np.array(means))
    for i in range(1, len(theta_list)):
        local = (log_m[i] - log_m[i-1]) / (log_t[i] - log_t[i-1])
        print(f"  local slope [{i-1}->{i}]: {local:.3f}")

    print("-" * 50)

if __name__ == '__main__':
    quick_test()
