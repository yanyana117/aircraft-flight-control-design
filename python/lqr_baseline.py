"""Discrete LQR baseline for the F-16A lateral model.

Python counterpart of the MATLAB lqrdjv-based SDR design: solve the
discrete algebraic Riccati equation for the ZOH-discretized plant and
simulate the closed loop under the same episode/cost definition as the
RL agents, so the returns are directly comparable.
"""

import numpy as np
from scipy.linalg import solve_discrete_are

from f16_env import AD, BD, Q_COST, R_COST, U_MAX, F16LateralEnv


def lqr_gain():
    p = solve_discrete_are(AD, BD, Q_COST, R_COST)
    return np.linalg.solve(R_COST + BD.T @ p @ BD, BD.T @ p @ AD)


def rollout(policy, env=None, seed=None):
    """Run one episode; policy maps state -> control. Returns (return, states, controls)."""
    env = env or F16LateralEnv(seed=seed)
    x = env.reset()
    states, controls, total = [x], [], 0.0
    done = False
    while not done:
        u = np.clip(policy(x), -U_MAX, U_MAX)  # record the applied (limited) command
        x, r, done, _ = env.step(u)
        states.append(x)
        controls.append(u)
        total += r
    return total, np.array(states), np.array(controls)


if __name__ == "__main__":
    k = lqr_gain()
    ret, xs, us = rollout(lambda x: -k @ x)
    print("LQR gain K =\n", np.round(k, 4))
    print(f"LQR episode return: {ret:.2f}")
    print(f"final |phi| = {np.rad2deg(abs(xs[-1, 3])):.3f} deg")
