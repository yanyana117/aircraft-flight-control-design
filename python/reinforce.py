"""REINFORCE from scratch (NumPy only) on the F-16A lateral model.

Gaussian policy u ~ N(Wx, diag(exp(2*log_std))), linear in the state --
the same structure as the LQR law, so the agent can in principle recover
it. Policy-gradient update uses reward-to-go with a moving-average
baseline and advantage normalization. No RL library is used.
"""

import numpy as np

from f16_env import F16LateralEnv, U_MAX


class GaussianLinearPolicy:
    def __init__(self, n_state=4, n_act=2, seed=0):
        rng = np.random.default_rng(seed)
        self.w = rng.normal(0.0, 0.01, size=(n_act, n_state))
        self.log_std = np.log(0.05 * U_MAX)  # exploration ~5% of surface limits

    def sample(self, x, rng):
        mean = self.w @ x
        std = np.exp(self.log_std)
        u = mean + std * rng.standard_normal(mean.shape)
        return u, mean

    def act_greedy(self, x):
        return self.w @ x


def train(episodes=4000, batch=16, lr=2e-3, gamma=0.995, seed=0, verbose=True):
    rng = np.random.default_rng(seed)
    env = F16LateralEnv(random_init=True, seed=seed)
    policy = GaussianLinearPolicy(seed=seed)
    history = []
    batch_grads, batch_returns = [], []

    for ep in range(episodes):
        x = env.reset()
        traj = []          # (x, u, mean, reward)
        done = False
        while not done:
            u, mean = policy.sample(x, rng)
            x_next, r, done, _ = env.step(u)
            traj.append((x, u, mean, r))
            x = x_next

        rewards = np.array([t[3] for t in traj])
        # reward-to-go
        g = np.zeros_like(rewards)
        acc = 0.0
        for i in range(len(rewards) - 1, -1, -1):
            acc = rewards[i] + gamma * acc
            g[i] = acc

        ep_return = rewards.sum()
        history.append(ep_return)
        batch_returns.append(g)
        batch_grads.append(traj)

        if len(batch_grads) == batch:
            all_g = np.concatenate(batch_returns)
            adv_mean, adv_std = all_g.mean(), all_g.std() + 1e-8
            grad_w = np.zeros_like(policy.w)
            var = np.exp(2 * policy.log_std)
            for trj, gs in zip(batch_grads, batch_returns):
                for (xs, us, ms, _), gt in zip(trj, gs):
                    adv = (gt - adv_mean) / adv_std
                    # d logpi / dW for Gaussian with mean Wx
                    grad_w += np.outer((us - ms) / var, xs) * adv
            policy.w += lr * grad_w / batch
            batch_grads, batch_returns = [], []

        if verbose and (ep + 1) % 500 == 0:
            recent = np.mean(history[-100:])
            print(f"episode {ep + 1:5d}  mean return (last 100): {recent:10.2f}")

    return policy, np.array(history)


if __name__ == "__main__":
    policy, hist = train()
    print("learned W =\n", np.round(policy.w, 4))
