import copy
import warnings
from typing import Any, Callable, Dict, Iterable, List, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rpy2
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from seaborn.regression import _LinearPlotter

numpy2ri.activate()
pandas2ri.activate()

binsreg = importr("binsreg")


def unwrap(robject, depth=1):

    if depth == 0:
        return robject
    else:
        return {
            n: unwrap(robject.rx[str(n)][0], depth=depth - 1) for n in robject.names
        }


def binsreg_data(
    x: np.array,
    y: np.array,
    w: np.array = rpy2.rinterface.NULL,
    by: np.array = rpy2.rinterface.NULL,
    dots: tuple = (0, 0),
    line: tuple = rpy2.rinterface.NULL,
    ci: tuple = rpy2.rinterface.NULL,
    cb: tuple = rpy2.rinterface.NULL,
    nbins: int = rpy2.rinterface.NULL,
    **kwargs,
) -> dict:

    if type(w) is not rpy2.rinterface.NULLType and w.ndim == 1:
        w = w.reshape((-1, 1))

    dots = np.array(dots)

    if type(line) is not rpy2.rinterface.NULLType:
        line = np.array(line)

    if type(ci) is not rpy2.rinterface.NULLType:
        ci = np.array(ci)

    if type(cb) is not rpy2.rinterface.NULLType:
        cb = np.array(cb)

    if type(by) is not rpy2.rinterface.NULLType:
        by = by.reshape((-1, 1))

    result = binsreg.binsreg(
        y.reshape((-1, 1)),
        x.reshape((-1, 1)),
        w,
        by=by,
        dots=dots,
        line=line,
        ci=ci,
        cb=cb,
        nbins=nbins,
        **kwargs,
    )

    dataplot = result.rx["data.plot"][0]

    return unwrap(dataplot, depth=2)


def binsplot(
    data_dict,
    ax,
    color=None,
    truncate=False,
    label=None,
    scatter_args=dict(),
    line_args=dict(),
    ci_args=dict(),
    poly_args=dict(),
):
    if color is None:
        lines, = ax.plot(0, 0)
        color = lines.get_color()
        lines.remove()
    else:
        color = color

    dots = data_dict["data.dots"]
    sc = ax.scatter(dots["x"], dots["fit"], color=color, label=label, **scatter_args)

    if "data.ci" in data_dict:
        ci = data_dict["data.ci"]
        ax.errorbar(
            x=ci["x"].values,
            y=dots["fit"].values,
            yerr=np.abs(
                ci[["ci.l", "ci.r"]].values - dots["fit"].values[:, np.newaxis]
            ).T,
            marker="",
            ls="",
            capsize=2,
            color=color,
            **ci_args,
        )

    if "data.line" in data_dict:
        line = data_dict["data.line"]
        if truncate:
            line = line.query("x >= @dots.x.min() and x <= @dots.x.max()")
        if "data.poly" not in data_dict:
            ax.plot(
                line["x"].values,
                line["fit"].values,
                color=color,
                lw=mpl.rcParams["lines.linewidth"] * 1.5,
                **line_args,
            )

    if "data.cb" in data_dict:
        cb = data_dict["data.cb"].query("x >= @line.x.min() and x <= @line.x.max()")
        ax.fill_between(
            cb["x"].values,
            cb["cb.l"].values,
            cb["cb.r"].values,
            color=color,
            alpha=0.15,
        )

    if "data.poly" in data_dict:
        poly = data_dict["data.poly"]
        if truncate:
            poly = poly.query("x >= @dots.x.min() and x <= @dots.x.max()")
        ax.plot(
            poly["x"].values,
            poly["fit"].values,
            color=color,
            lw=mpl.rcParams["lines.linewidth"] * 1.5,
            **poly_args,
        )

    if "data.polyci" in data_dict:
        polyci = data_dict["data.polyci"].query(
            "x >= @poly.x.min() and x <= @poly.x.max()"
        )
        ax.fill_between(
            polyci["x"].values,
            polyci["polyci.l"].values,
            polyci["polyci.r"].values,
            color=color,
            alpha=0.15,
        )
    return ax


class _BinnedRegressionPlotter(_LinearPlotter):
    def __init__(
        self,
        x: Union[str, np.array],
        y: Union[str, np.array],
        covariates: Union[str, List] = None,
        hue: Union[str, np.array] = None,
        data: pd.DataFrame = None,
        dots: tuple = (0, 0),
        line: tuple = rpy2.rinterface.NULL,
        ci: tuple = rpy2.rinterface.NULL,
        cb: tuple = rpy2.rinterface.NULL,
        nbins: int = rpy2.rinterface.NULL,
        color=None,
        **kwargs,
    ):
        if type(covariates) is str:
            covariates = [covariates]
        self.covariate_names = covariates
        if data is not None and type(covariates) is list:
            covariates = data[covariates]
        self.establish_variables(data, x=x, y=y, covariates=covariates, hue=hue)
        self.dropna("x", "y", "covariates", "hue")
        self.color = color

        self.data_plot = binsreg_data(
            x=np.array(self.x),
            y=np.array(self.y),
            w=np.array(self.covariates)
            if self.covariates is not None
            else rpy2.rinterface.NULL,
            by=np.array(self.hue).astype(str)
            if self.hue is not None
            else rpy2.rinterface.NULL,
            dots=dots,
            line=line,
            ci=ci,
            cb=cb,
            nbins=nbins,
            **kwargs,
        )

    def plot(
        self,
        ax,
        truncate=False,
        legend_args=dict(frameon=False),
        scatter_args=dict(),
        line_args=dict(),
        poly_args=dict(),
        ci_args=dict(),
        raw=True,
    ):

        for grp in self.data_plot:
            if raw:

                x = self.x[self.hue == grp[6:]] if self.hue is not None else self.x
                x = np.array(x)
                y = self.y[self.hue == grp[6:]] if self.hue is not None else self.y
                y = np.array(y)

                if self.covariates is not None:
                    mat = np.array(
                        self.covariates[self.hue == grp[6:]]
                        if self.hue is not None
                        else self.covariates
                    )
                    mat = np.hstack([np.ones((mat.shape[0], 1)), mat])
                    residual_proj = (
                        np.eye(mat.shape[0]) - mat @ np.linalg.inv(mat.T @ mat) @ mat.T
                    )
                    res_x = residual_proj @ x + x.mean()
                    res_y = residual_proj @ y + y.mean()
                    ax.scatter(res_x, res_y, alpha=0.1, marker=".")
                else:
                    ax.scatter(x, y, alpha=0.1, marker=".")

            if len(self.data_plot) > 1:
                label = grp[6:]
            else:
                label = None
            binsplot(
                self.data_plot[grp],
                ax=ax,
                truncate=truncate,
                label=label,
                scatter_args=scatter_args,
                line_args=line_args,
                poly_args=poly_args,
                ci_args=ci_args,
            )
        if len(self.data_plot) > 1:
            ax.legend(**legend_args)

        # Label the axes
        if hasattr(self.x, "name"):
            ax.set_xlabel(self.x.name)
        if hasattr(self.y, "name"):
            if self.covariates is not None:
                ax.set_ylabel(
                    f"{self.y.name} (adjusted for {', '.join(self.covariate_names)})"
                )
            else:
                ax.set_ylabel(f"{self.y.name}")
        return ax


def binscatterplot(
    x: Union[str, np.array],
    y: Union[str, np.array],
    fit_spline: bool = True,
    fit_reg: bool = False,
    disp_ci: bool = True,
    covariates: Union[str, List] = None,
    hue: Union[str, np.array] = None,
    data: pd.DataFrame = None,
    dots: tuple = (0, 0),
    line: tuple = rpy2.rinterface.NULL,
    ci: tuple = rpy2.rinterface.NULL,
    cb: tuple = rpy2.rinterface.NULL,
    nbins: int = rpy2.rinterface.NULL,
    ax: plt.Axes = None,
    color: str = None,
    truncate: bool = False,
    legend_args: dict = dict(frameon=False),
    scatter_args: dict = dict(),
    line_args: dict = dict(),
    poly_args: dict = dict(),
    ci_args: dict = dict(),
    raw: bool = True,
    **kwargs,
) -> plt.Axes:
    if ax is None:
        ax = plt.gca()

    if fit_spline:
        line = (3, 3)
        cb = (3, 3)

    if disp_ci:
        ci = (0, 0)

    if fit_reg:
        kwargs["polyreg"] = 1

    plotter = _BinnedRegressionPlotter(
        x,
        y,
        covariates=covariates,
        hue=hue,
        data=data,
        dots=dots,
        line=line,
        ci=ci,
        cb=cb,
        nbins=nbins,
        color=color,
        **kwargs,
    )

    ax = plotter.plot(
        ax,
        truncate=truncate,
        legend_args=legend_args,
        scatter_args=scatter_args,
        line_args=line_args,
        poly_args=poly_args,
        ci_args=ci_args,
        raw=raw,
    )

    return ax


binscatterplot.__doc__ = f"""\
Plot data via binned scatterplot.

See https://arxiv.org/pdf/1902.09608.pdf for details.
Site and documentation for binsreg: https://sites.google.com/site/nppackages/binsreg/.

Parameters
----------
x : Union[str, np.array]
    Independent variable
y : Union[str, np.array]
    Dependent variable
covariates : Union[str, List], optional
    List of additional covariates to control for, by default None
hue : Union[str, np.array], optional
    Variable for groupings of data, by default None
fit_spline : bool, optional
    Equivalent to setting line=(3, 3) and cb=(3, 3) for fitting cubic splines,
    If fit_spline is True, line and cb arguments are ignored; by default True
fit_reg : bool, optional
    Equivalent to setting polyreg=1. If True, polyreg is ignored, by default False
disp_ci : bool, optional
    Equivalent to setting ci=(3, 3). If True, ci is ignored, by default True
data : pd.DataFrame, optional
    DataFrame for the data, must be specified if x, y are not iterables, by default None
raw : bool, optional
    Whether or not raw data is plotted in the background, by default True.
    Mean-shifted residualized data is plotted in the presence of additional covariates:
    The residuals of z (for z being x, y) regressed on intercept and covariates are
    recentered at the original mean of z (instead of zero) and plotted.
dots : tuple, optional
    Argument passed to binsreg in R, by default (0, 0)
line : tuple, optional
    Argument passed to binsreg in R, by default rpy2.rinterface.NULL
ci : tuple, optional
    Argument passed to binsreg in R, by default rpy2.rinterface.NULL
cb : tuple, optional
    Argument passed to binsreg in R, by default rpy2.rinterface.NULL
nbins : int, optional
    Argument passed to binsreg in R, by default rpy2.rinterface.NULL
ax : matplotlib.pyplot.Axes, optional
    Axes object for the plot, by default None
color : str, optional
    Color of the plot, by default None
truncate : bool, optional
    Whether to truncate the line fits (spline or global regression) at the
    left and right most bins, by default False
legend_args : dict, optional
    Arguments passed to ax.legend(), by default dict(frameon=False)
scatter_args : dict, optional
    Arguments passed to ax.scatter for the binned dots, by default dict()
line_args : dict, optional
    Arguments passed to ax.plot for the spline fit, by default dict()
poly_args : dict, optional
    Arguments passed to ax.plot for the global regression fit, by default dict()
ci_args : dict, optional
    Arguments passed to ax.errorbar for the CI on binned points, by default dict()
**kwargs:
    Additional arguments passed to binsreg in R

Returns
-------
ax: plt.Axes
    The Axes for the plot

`binsreg` documentation in R
----------------------------
{binsreg.binsreg.__doc__}
"""
