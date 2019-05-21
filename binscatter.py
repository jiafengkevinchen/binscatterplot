import warnings
import copy
from math import factorial

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from numpy.linalg import LinAlgError
from scipy.stats import norm
from scipy.special import binom

from seaborn.regression import _LinearPlotter
from patsy import bs
from statsmodels.distributions.empirical_distribution import ECDF


class _BinnedRegressionPlotter(_LinearPlotter):
    def __init__(self, x, y, spline_order=3, n_bins="rot", color=None,
                 data=None, covariates=None):
        """
        Parameters
        ----------
        x : str
            The independent variable column name in data
        y : str
            The dependent variable column name in data
        spline_order : int, optional
            Order of the B-spline fitting, default 3 for cubic splines
        n_bins : int or str, optional
            If int, then overrides automatic bin selection. Otherwise use
            IMSE rule-of-thumb selector
        color : str, optional
            Matplotlib color codes
        data : pandas.DataFrame, optional
            DataFrame for data
        covariates : iterable of str, optional
            List of column names for covariates
        """
        self.covariate_names = covariates
        if data is not None and type(covariates) is list:
            covariates = data[covariates]
        self.establish_variables(data, x=x, y=y, covariates=covariates)
        self.dropna("x", "y", "covariates")
        self.x_range = self.x.min(), self.x.max()
        self.spline_order = spline_order
        self.n_bins = self.optimal_imse_bins(
            n_bins) if type(n_bins) is not int else n_bins
        self.color = color

        self.get_bins()

    def optimal_imse_bins(self, n_bins):
        max_bins = len(set(ECDF(self.x)(self.x)))
        if n_bins.lower() == "rot":
            powers = np.array(
                [self.x ** i for i in range(self.spline_order + 3)]).T
            rhs = self._concat_with_covariates(powers)

            coefs1, _, _, _ = np.linalg.lstsq(rhs, self.y ** 2, rcond=None)
            coefs2, _, _, _ = np.linalg.lstsq(rhs, self.y, rcond=None)
            sigma2 = rhs @ coefs1 - (rhs @ coefs2) ** 2

            V = (self.spline_order + 1) * sigma2.mean()

            mu_sum = np.array([factorial(self.spline_order + t) / factorial(t)
                               * coefs2[self.spline_order + t]
                               * (self.x ** t)
                               for t in range(3)]).sum(axis=0)

            fz = np.clip(norm.pdf(self.x, loc=self.x.mean(), scale=self.x.std()),
                         a_min=norm.pdf(2, loc=0, scale=1),
                         a_max=None)

            B = (mu_sum ** 2 / (fz ** (2 * (self.spline_order + 1)))).mean()
            b_cons = (1 / (2 * (self.spline_order + 1) + 1)
                      / factorial(self.spline_order + 1) ** 2
                      / binom(2*(self.spline_order + 1), self.spline_order + 1) ** 2)
            B *= b_cons
            bins = int(np.ceil((2 * (self.spline_order + 1) * B / V) ** (1/(2 * self.spline_order + 3)) * len(self.x)
                               ** (1/(2 * self.spline_order + 3))))
        else:
            raise NotImplementedError
        if bins > max_bins // 2:
            warnings.warn("Optimal number of bins suggested is too large, reverting to half the maximum number of bins",
                          RuntimeWarning)
        return min(bins, max_bins // 2)

    def get_bins(self, x=None, y=None):
        J = self.n_bins
        if x is None:
            x = self.x
        if y is None:
            y = self.y

        percentile_targets = np.arange(J + 1) / J
        bin_boundaries = np.quantile(x, percentile_targets)
        bin_boundaries[-1] += 1e-10
        self.bin_boundaries = bin_boundaries

        group = np.greater_equal.outer(x, bin_boundaries).sum(axis=1)
        self.x_group = group

        x_centers = pd.Series(x).groupby(group).mean()
        y_centers_mean = pd.Series(y).groupby(group).mean()
        return x_centers, y_centers_mean

    def _spline_transform(self, x):
        splines = bs(x, knots=self.bin_boundaries[1:-1],
                     degree=self.spline_order, include_intercept=True,
                     lower_bound=self.bin_boundaries[0],
                     upper_bound=self.bin_boundaries[-1])
        return splines

    def _concat_with_covariates(self, columns):
        if self.covariates is not None:
            return np.hstack([columns, self.covariates])
        else:
            return columns

    def _get_spline_design(self):
        splines = self._spline_transform(self.x)
        n, k = splines.shape
        rhs = np.array(self._concat_with_covariates(splines))
        lhs = np.array(self.y)
        return lhs, rhs, (n, k)

    def fit_spline(self):
        lhs, rhs, (n, k) = self._get_spline_design()
        coefs, _, _, _ = np.linalg.lstsq(rhs, lhs, rcond=None)
        self.coefs = coefs[:k].copy()

        if self.covariates is not None:
            self.covariate_coefs = coefs[k:].copy()
            self.covariate_mean = np.array(self.covariates) @ self.covariate_coefs
        else:
            self.covariate_mean = np.zeros(len(self.x))

    def predict_spline(self, x):
        return np.array(self._spline_transform(x)) @ self.coefs

    def fit_predict_spline(self, x):
        self.fit_spline()
        return self.predict_spline(x)

    def fit_linear(self):
        if self.covariates is not None:
            rhs = np.array(self.covariates)
            rhs = np.hstack([np.ones(len(rhs))[:, np.newaxis], rhs])

            coefs_yz, _, _, _ = np.linalg.lstsq(rhs, self.y, rcond=None)
            coefs_xz, _, _, _ = np.linalg.lstsq(rhs, self.x, rcond=None)

            self.fw_residualized_y = self.y - rhs @ coefs_yz + self.y.mean()
            self.fw_residualized_x = self.x - rhs @ coefs_xz + self.x.mean()

    def spline_se(self, x):
        basis = self._spline_transform(x)
        basis_data = self._spline_transform(self.x).values
        Q_inv = np.linalg.inv(basis_data.T @ basis_data / len(self.x))
        Sigma = np.array([np.outer(b, b) * ((self.y[i] - np.dot(b, self.coefs) - self.covariate_mean[i]) ** 2)
                          for i, b in enumerate(basis_data)]) / len(self.x)
        Sigma = Sigma.sum(axis=0)
        middle = Q_inv @  Sigma @ Q_inv
        var = np.diag(basis @ middle @ basis.T)
        return (var / len(self.x)) ** .5

    def plot(self, ax, bin_boundary=False, alpha_data=0,
             plot_residualize=False, scatter_kws=dict(), ci_kws=dict()):

        if plot_residualize:
            self.fit_linear()
            x_centers_fw, y_centers_fw = self.get_bins(x=self.fw_residualized_x,
                                                       y=self.fw_residualized_y)
            ax.scatter(x_centers_fw, y_centers_fw, marker='x')

        x_centers, y_centers_mean = self.get_bins()

        if bin_boundary:
            for b in self.bin_boundaries[1:-1]:
                ax.axvline(b, ls='--', linewidth=0.5)

        if self.spline_order is not None:

            if self.color is None:
                lines, = ax.plot(self.x.mean(), self.y.mean())
                color = lines.get_color()
                lines.remove()
            else:
                color = self.color

            # Ensure that color is hex to avoid matplotlib weidness
            color = mpl.colors.rgb2hex(mpl.colors.colorConverter.to_rgb(color))

            # Plot fit line
            xs = np.linspace(
                self.bin_boundaries[0], self.bin_boundaries[-1], 1000)
            ys = self.fit_predict_spline(xs)
            ax.plot(xs, ys, color=color)

            try:
                se = self.spline_se(xs)
                ax.fill_between(xs, ys + 1.96 * se, ys - 1.96 * se,
                                facecolor=color, alpha=0.15)
            except LinAlgError:
                warnings.warn(
                    "Computation of error band encountered LinAlgError",
                    RuntimeWarning
                )

            y_centers = self.predict_spline(x_centers)
            ax.scatter(x_centers, y_centers, color=color)

            if alpha_data > 0:

                if "marker" not in scatter_kws:
                    scatter_kws["marker"] = '.'
                if "color" not in scatter_kws:
                    scatter_color = np.array(mpl.colors.hex2color(color)) / 1.5
                    scatter_color = mpl.colors.rgb2hex(scatter_color)
                    scatter_kws["color"] = scatter_color
                scatter_kws["alpha"] = alpha_data

                ax.scatter(np.array(self.x), np.array(self.y) - self.covariate_mean,
                           **scatter_kws)

            groups = pd.Series(np.array(self.y) -
                               self.covariate_mean).groupby(self.x_group)
            ci = (2 * groups.std().values) / (groups.size().values ** .5)
            ci = zip(y_centers - ci, y_centers + ci)

            if "color" not in ci_kws:
                ci_kws['color'] = color
            if "linewidth" not in ci_kws:
                ci_kws["linewidth"] = mpl.rcParams["lines.linewidth"] * 1.75

            for x, ci in zip(x_centers, ci):
                ax.plot([x, x], ci, **ci_kws)

        # Label the axes
        if hasattr(self.x, "name"):
            ax.set_xlabel(self.x.name)
        if hasattr(self.y, "name"):
            if self.covariates is not None:
                ax.set_ylabel(
                    f"{self.y.name} (adjusted for {', '.join(self.covariate_names)})")


def binscatterplot(x, y, covariates=None, spline_order=3,
                   n_bins="rot", color=None,
                   data=None, ax=None, scatter_kws=None,
                   bin_boundary=False,
                   plot_residualize=False,
                   ci_kws=None,
                   alpha_data=0):
    plotter = _BinnedRegressionPlotter(x, y, covariates=covariates,
                                      spline_order=spline_order,
                                      n_bins=n_bins, color=color,
                                      data=data)
    if ax is None:
        ax = plt.gca()

    scatter_kws = dict() if scatter_kws is None else copy.copy(scatter_kws)
    ci_kws = dict() if ci_kws is None else copy.copy(ci_kws)

    plotter.plot(bin_boundary=bin_boundary, alpha_data=alpha_data, ci_kws=ci_kws,
                 plot_residualize=plot_residualize, scatter_kw=scatter_kws)
