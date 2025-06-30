import numpy as np

from scipy.special import gammaln
from threadpoolctl import threadpool_limits


class ZTP:
    def __init__(self, x):
        '''
        "Equivalence of Truncated Count Mixture Distributions and Mixtures of Truncated Count Distributions"
        https://doi.org/10.1111/j.1541-0420.2006.00565.x
        '''
        self.x = np.asarray(x)
        self.n = len(x)
        self.e = np.finfo(float).eps

    @staticmethod
    def _pmf(x, lmd):
        return np.exp(x * np.log(lmd) - lmd - gammaln(x + 1)) / (1 - np.exp(-lmd))

    @threadpool_limits.wrap(limits=1)
    def fit(self, components, max_iter=np.inf, tol=1e-5):
        '''
        Fit finite mixture of ZTP distributions with EM.
        '''
        lmds = np.quantile(self.x, [0.5] if components == 1 else [i / (components - 1) for i in range(components)])
        weights = np.repeat(1 / components, components)

        old_logl = -np.inf
        old_lmds = lmds.copy()

        max_iter = np.iinfo(np.int32).max if np.isinf(max_iter) else max_iter
        for _ in range(max_iter):

            ## get probabilities
            com_probs = self._pmf(self.x[:, None], lmds[None, :]) * weights
            mix_probs = com_probs.sum(axis=1) + self.e

            ## test for log-likelihood convergence
            logl = np.sum(np.log(mix_probs))
            if np.isclose(logl, old_logl, atol=tol, rtol=0):
                break
            old_logl = logl

            ## update lmds
            R = com_probs / mix_probs[:, None] + self.e
            Z = R.T.dot(self.x) / R.sum(axis=0)
            for _ in range(max_iter):
                lmds = Z * (1 - np.exp(-old_lmds))
                if np.allclose(lmds, old_lmds, atol=tol, rtol=0):
                    break
                old_lmds = lmds.copy()

            ## update weights
            weights = R.sum(axis=0) / self.n

        bic = -2 * logl + (2 * components - 1) * np.log(len(self.x))
        return lmds, weights / weights.sum(), bic
