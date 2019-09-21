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
):

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
                if self.covariates is not None:
                    mat = self.covariates.values[self.hue == grp[6:]]
                    mat = np.hstack([np.ones((mat.shape[0], 1)), mat])
                    residual_proj = (
                        np.eye(mat.shape[0]) - mat @ np.linalg.inv(mat.T @ mat) @ mat.T
                    )
                    res_x = (
                        residual_proj @ self.x.values[self.hue == grp[6:]]
                        + self.x.values[self.hue == grp[6:]].mean()
                    )
                    res_y = (
                        residual_proj @ self.y.values[self.hue == grp[6:]]
                        + self.y.values[self.hue == grp[6:]].mean()
                    )
                    ax.scatter(
                        res_x, res_y,
                        alpha=0.1,
                        marker=".",
                    )
                else:
                    ax.scatter(
                        self.x[self.hue == grp[6:]].values,
                        self.y[self.hue == grp[6:]].values,
                        alpha=0.1,
                        marker=".",
                    )
            binsplot(
                self.data_plot[grp],
                ax=ax,
                truncate=truncate,
                label=grp[6:],
                scatter_args=scatter_args,
                line_args=line_args,
                poly_args=poly_args,
                ci_args=ci_args,
            )

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
    ax=None,
    color=None,
    truncate=False,
    legend_args=dict(frameon=False),
    scatter_args=dict(),
    line_args=dict(),
    poly_args=dict(),
    ci_args=dict(),
    raw: bool = True,
    **kwargs,
):
    if ax is None:
        ax = plt.gca()

    if fit_spline:
        line = (3, 3)
        cb = (3, 3)

    if disp_ci:
        ci = (3, 3)

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
