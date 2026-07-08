"""Tabular Q-learning from scratch on the roll-axis subsystem.

Reduced-order plant (p, phi) extracted from the same F-16A lateral model
(roll-rate damping and aileron effectiveness terms), discretized on a
coarse state grid with a small set of discrete aileron commands --
the classic tabular setting, kept honest by naming it a reduced model.
"""

import numpy as np
from scipy.linalg import expm

from f16_env import T_SAMPLE

# roll subsystem from A_LAT/B_LAT: p_dot = -1.179 p - 5.935 da, phi_dot = p
A2 = np.array([[-1.179, 0.0], [1.0, 0.0]])
B2 = np.array([[-5.935], [0.0]])
aug = np.zeros((3, 3))
aug[:2, :2] = A2
aug[:2, 2:] = B2
PHI = expm(aug * T_SAMPLE)
AD2, BD2 = PHI[:2, :2], PHI[:2, 2:]

P_BINS = np.linspace(-1.2, 1.2, 17)      # roll rate grid (rad/s)
PHI_BINS = np.linspace(-0.6, 0.6, 17)    # bank angle grid (rad)
ACTIONS = np.deg2rad([-10.0, -5.0, 0.0, 5.0, 10.0])  # aileron commands
Q_W, R_W = 50.0, 1.0                     # cost weights: phi^2 and u^2


def to_idx(x):
    p_i = np.clip(np.digitize(x[0], P_BINS), 0, len(P_BINS) - 1)
    phi_i = np.clip(np.digitize(x[1], PHI_BINS), 0, len(PHI_BINS) - 1)
    return p_i, phi_i


def step(x, u):
    x_next = AD2 @ x + BD2 @ [u]
    reward = -(Q_W * x[1] ** 2 + R_W * u ** 2)
    return x_next, reward


def train(episodes=10000, horizon=200, alpha0=0.1, alpha1=0.02, gamma=0.99,
          eps_start=1.0, eps_end=0.1, seed=0, verbose=True):
    rng = np.random.default_rng(seed)
    q = np.zeros((len(P_BINS), len(PHI_BINS), len(ACTIONS)))
    history = []
    for ep in range(episodes):
        frac = ep / episodes
        eps = eps_start + (eps_end - eps_start) * frac
        alpha = alpha0 + (alpha1 - alpha0) * frac
        x = np.array([0.0, rng.uniform(0.1, 0.35) * rng.choice([-1, 1])])
        total = 0.0
        for _ in range(horizon):
            s = to_idx(x)
            a = rng.integers(len(ACTIONS)) if rng.random() < eps else np.argmax(q[s])
            x_next, r = step(x, ACTIONS[a])
            s_next = to_idx(x_next)
            q[s + (a,)] += alpha * (r + gamma * np.max(q[s_next]) - q[s + (a,)])
            x, total = x_next, total + r
        history.append(total)
        if verbose and (ep + 1) % 2000 == 0:
            print(f"episode {ep + 1:5d}  mean return (last 200): {np.mean(history[-200:]):9.2f}")
    return q, np.array(history)


def greedy_eval(q, phi0=0.175, horizon=200):
    """Greedy rollout from the nominal 10-degree bank disturbance."""
    x = np.array([0.0, phi0])
    total = 0.0
    for _ in range(horizon):
        a = int(np.argmax(q[to_idx(x)]))
        x, r = step(x, ACTIONS[a])
        total += r
    return total, x


if __name__ == "__main__":
    q, hist = train()
    ret, xf = greedy_eval(q)
    print(f"greedy eval return: {ret:.2f}, final phi = {np.rad2deg(xf[1]):.2f} deg")
