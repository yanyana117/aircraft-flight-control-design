"""Gym-style environment for the F-16A lateral/directional model.

The plant is the same linear model used in the AERO 625 MATLAB controller
designs in ../controllers (flight condition M = 0.18, U1 = 205 ft/s):
states [beta, p, r, phi] (sideslip, roll rate, yaw rate, bank angle),
inputs [aileron, rudder]. The heading state (psi, a pure integrator) is
dropped, matching the reduced-model option noted in the original scripts.

Dynamics are discretized exactly (zero-order hold) at the same T = 0.05 s
sample time as the MATLAB digital designs. Reward is the negative LQR
stage cost -(x'Qx + u'Ru), so "return" is directly comparable with the
LQR baseline cost.
"""

import numpy as np
from scipy.linalg import expm

# F-16A lateral/directional model from controllers/sdr_regulator.m
A_LAT = np.array([
    [-0.132,   0.324,   -0.94,    0.149],
    [-10.614, -1.179,    1.0023,  0.0],
    [0.997,   -0.00182, -0.259,   0.0],
    [0.0,      1.0,      0.34,    0.0],
])
B_LAT = np.array([
    [0.0069,  0.0189],
    [-5.935,  1.203],
    [-0.122, -0.614],
    [0.0,     0.0],
])

T_SAMPLE = 0.05          # s, same as the MATLAB designs
# Cost weights follow the AERO 625 SDR design intent: penalize sideslip
# and bank angle, light weight on the rates, identity on the controls.
Q_COST = np.diag([25.0, 1.0, 1.0, 50.0])
R_COST = np.eye(2)
# Surface position limits (rad): +-21.5 deg aileron, +-30 deg rudder
U_MAX = np.array([np.deg2rad(21.5), np.deg2rad(30.0)])
X_LIMIT = 5.0            # abort episode if any |state| exceeds this (rad, rad/s)


def discretize(a, b, t):
    """Exact ZOH discretization via the augmented matrix exponential."""
    n, m = b.shape
    aug = np.zeros((n + m, n + m))
    aug[:n, :n] = a
    aug[:n, n:] = b
    phi = expm(aug * t)
    return phi[:n, :n], phi[:n, n:]


AD, BD = discretize(A_LAT, B_LAT, T_SAMPLE)


class F16LateralEnv:
    """Bank-angle disturbance regulation, 10 s horizon (200 steps)."""

    def __init__(self, horizon=200, phi0_deg=10.0, random_init=False, seed=None):
        self.horizon = horizon
        self.phi0 = np.deg2rad(phi0_deg)
        self.random_init = random_init
        self.rng = np.random.default_rng(seed)
        self.x = None
        self.k = 0

    def reset(self):
        phi0 = self.phi0
        if self.random_init:
            phi0 = self.rng.uniform(0.5, 1.5) * self.phi0 * self.rng.choice([-1, 1])
        self.x = np.array([0.0, 0.0, 0.0, phi0])
        self.k = 0
        return self.x.copy()

    def step(self, u):
        u = np.clip(np.asarray(u, dtype=float), -U_MAX, U_MAX)
        cost = self.x @ Q_COST @ self.x + u @ R_COST @ u
        self.x = AD @ self.x + BD @ u
        self.k += 1
        diverged = np.any(np.abs(self.x) > X_LIMIT)
        done = self.k >= self.horizon or diverged
        reward = -cost - (500.0 if diverged else 0.0)
        return self.x.copy(), reward, done, {"diverged": diverged}
