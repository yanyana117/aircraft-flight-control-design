"""Classical control vs reinforcement learning on the same F-16A plant.

Trains the from-scratch REINFORCE agent and tabular Q-learning, then
compares against the discrete LQR baseline on the nominal 10-degree
bank-angle disturbance episode. Saves learning curves and trajectory
plots to ../results and prints a summary table.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from f16_env import F16LateralEnv, T_SAMPLE
from lqr_baseline import lqr_gain, rollout
import reinforce
import q_learning


def main():
    # --- baselines and agents -------------------------------------------
    k = lqr_gain()
    lqr_ret, lqr_xs, lqr_us = rollout(lambda x: -k @ x)
    print(f"LQR baseline return: {lqr_ret:.2f}")

    policy, pg_hist = reinforce.train()
    rl_ret, rl_xs, rl_us = rollout(policy.act_greedy)
    print(f"REINFORCE greedy return: {rl_ret:.2f}")

    q_table, ql_hist = q_learning.train()
    ql_ret, ql_xf = q_learning.greedy_eval(q_table)
    print(f"Q-learning greedy eval: {ql_ret:.2f}, final phi = {np.rad2deg(ql_xf[1]):.2f} deg")

    # --- learning curve ---------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    smooth = np.convolve(pg_hist, np.ones(100) / 100, mode="valid")
    ax[0].plot(pg_hist, alpha=0.25, label="episode return")
    ax[0].plot(np.arange(len(smooth)) + 99, smooth, label="100-ep moving avg")
    ax[0].axhline(lqr_ret, color="k", ls="--", label=f"LQR baseline ({lqr_ret:.0f})")
    ax[0].set(title="REINFORCE on F-16A lateral model", xlabel="episode",
              ylabel="return (neg. LQR cost)")
    ax[0].legend()
    ql_smooth = np.convolve(ql_hist, np.ones(200) / 200, mode="valid")
    ax[1].plot(ql_hist, alpha=0.2, label="episode return")
    ax[1].plot(np.arange(len(ql_smooth)) + 199, ql_smooth, label="200-ep moving avg")
    ax[1].set(title="Tabular Q-learning (roll subsystem)", xlabel="episode",
              ylabel="return")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig("../results/learning_curves.png", dpi=150)

    # --- trajectories ------------------------------------------------------
    t_lqr = np.arange(len(lqr_xs)) * T_SAMPLE
    t_rl = np.arange(len(rl_xs)) * T_SAMPLE
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for i, (name, idx) in enumerate([("bank angle $\\phi$ (deg)", 3),
                                     ("sideslip $\\beta$ (deg)", 0)]):
        ax[i].plot(t_lqr, np.rad2deg(lqr_xs[:, idx]), label="LQR")
        ax[i].plot(t_rl, np.rad2deg(rl_xs[:, idx]), "--", label="REINFORCE")
        ax[i].set(title=name, xlabel="time (s)")
        ax[i].legend()
    ax[2].plot(t_lqr[:-1], np.rad2deg(lqr_us[:, 0]), label="LQR aileron")
    ax[2].plot(t_rl[:-1], np.rad2deg(rl_us[:, 0]), "--", label="RL aileron")
    ax[2].set(title="aileron command (deg)", xlabel="time (s)")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig("../results/trajectories.png", dpi=150)

    # --- summary -----------------------------------------------------------
    gap = 100 * (rl_ret - lqr_ret) / abs(lqr_ret)
    print("\nSummary (10-degree bank-angle disturbance, 10 s episode)")
    print(f"{'controller':<22}{'return':>12}")
    print(f"{'discrete LQR':<22}{lqr_ret:>12.2f}")
    print(f"{'REINFORCE (greedy)':<22}{rl_ret:>12.2f}   ({gap:+.1f}% vs LQR)")
    print(f"{'Q-learning (greedy, roll subsystem)':<36}{ql_ret:>12.2f}   final phi {np.rad2deg(ql_xf[1]):+.2f} deg")


if __name__ == "__main__":
    main()
