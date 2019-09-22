# Binscatter implementation in Python

A Python wrapper of `binsreg` in R for binned scatterplots with automatic bandwidth selection and
nonparametric fitting (See [Cattaneo, Crump, Farrell, and Feng](https://arxiv.org/pdf/1902.09608.pdf)). 

- Uses `rpy2` and handles the input and output, so the user wouldn't have to worry about various R objects
- Uses Pythonic plotting capabilities like `seaborn` 
- Mimicks `seaborn.regplot` in usage 

## Installation

`pip install git+https://github.com/jiafengkevinchen/binscatterplot`


## Minimal working example

![image](https://user-images.githubusercontent.com/24930289/65379164-aea61780-dc91-11e9-9d0b-f497d0d917ca.png)

## Documentation

```
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
Python representation of an R function.
description
-----------


 binsreg  implements binscatter estimation with robust inference proposed and plots, following the
             results in  https://arxiv.org/abs/1902.09608 Cattaneo, Crump, Farrell and Feng (2019a) .
             Binscatter provides a flexible way of describing
             the mean relationship between two variables, after possibly adjusting for other covariates, based on
             partitioning/binning of the independent variable of interest. The main purpose of this function is to
             generate binned scatter plots with curve estimation with robust pointwise confidence intervals and
             uniform confidence band.  If the binning scheme is not set by the user, the companion function
              binsregselect  is used to implement binscatter in a data-driven (optimal)
             way. Hypothesis testing about the regression function can also be conducted via the companion
             function  binsregtest .



binsreg(
    y,
    x,
    w = rinterface.NULL,
    deriv = 0.0,
    dots = c,
    dotsgrid = 0.0,
    dotsgridmean = <rpy2.rinterface.ListSexpVector object at 0x113ef9048> [RTYPES.VECSXP],
    line = rinterface.NULL,
    linegrid = 20.0,
    ci = rinterface.NULL,
    cigrid = 0.0,
    cigridmean = <rpy2.rinterface.ListSexpVector object at 0x113ef9048> [RTYPES.VECSXP],
    cb = rinterface.NULL,
    cbgrid = 20.0,
    polyreg = rinterface.NULL,
    polyreggrid = 20.0,
    polyregcigrid = 0.0,
    by = rinterface.NULL,
    bycolors = rinterface.NULL,
    bysymbols = rinterface.NULL,
    bylpatterns = rinterface.NULL,
    legendTitle = rinterface.NULL,
    legendoff = <rpy2.rinterface.ListSexpVector object at 0x113ef9708> [RTYPES.VECSXP],
    testmodel = c,
    testmodelparfit = rinterface.NULL,
    testmodelpoly = rinterface.NULL,
    testshape = c,
    testshapel = rinterface.NULL,
    testshaper = rinterface.NULL,
    testshape2 = rinterface.NULL,
    nbins = rinterface.NULL,
    binspos = qs,
    binsmethod = dpi,
    nbinsrot = rinterface.NULL,
    samebinsby = <rpy2.rinterface.ListSexpVector object at 0x113ef9708> [RTYPES.VECSXP],
    nsims = 500.0,
    simsgrid = 20.0,
    simsseed = 666.0,
    vce = HC1,
    cluster = rinterface.NULL,
    level = 95.0,
    noplot = <rpy2.rinterface.ListSexpVector object at 0x113ef9408> [RTYPES.VECSXP],
    dfcheck = c,
    masspoints = on,
    weights = rinterface.NULL,
    subset = rinterface.NULL,
)

y :  outcome variable. A vector. ,

x :  independent variable of interest. A vector. ,

w :  control variables. A matrix or a vector. ,

deriv :  derivative order of the regression function for estimation, testing and plotting.The default is `deriv=0`, which corresponds to the function itself. ,

dots :  a vector. `dots=c(p,s)` sets a piecewise polynomial of degree `p` with `s` smoothness constraints forpoint estimation and plotting as "dots". The default is `dots=c(0,0)`, which corresponds topiecewise constant (canonical binscatter) ,

dotsgrid :  number of dots within each bin to be plotted. Given the choice, these dots are point estimatesevaluated over an evenly-spaced grid within each bin. The default is `dotsgrid=0`, and onlythe point estimates at the mean of `x` within each bin are presented. ,

dotsgridmean :  If true, the dots corresponding to the point estimates evaluated at the mean of `x` within each binare presented. By default, they are presented, i.e., `dotsgridmean=T`. ,

line :  a vector. `line=c(p,s)` sets a piecewise polynomial of degree `p` with `s` smoothness constraintsfor plotting as a "line".  By default, the line is not included in the plot unless explicitlyspecified.  Recommended specification is `line=c(3,3)`, which adds a cubic B-spline estimateof the regression function of interest to the binned scatter plot. ,

linegrid :  number of evaluation points of an evenly-spaced grid within each bin used for evaluation ofthe point estimate set by the `line=c(p,s)` option. The default is `linegrid=20`,which corresponds to 20 evenly-spaced evaluation points within each bin for fitting/plotting the line. ,

ci :  a vector. `ci=c(p,s)` sets a piecewise polynomial of degree `p` with `s` smoothness constraints used forconstructing confidence intervals. By default, the confidence intervals are not included in the plotunless explicitly specified.  Recommended specification is `ci=c(3,3)`, which adds confidenceintervals based on cubic B-spline estimate of the regression function of interest to the binned scatter plot. ,

cigrid :  number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the pointestimate set by the `ci=c(p,s)` option. The default is `cigrid=1`, which corresponds to 1evenly-spaced evaluation point within each bin for confidence interval construction. ,

cigridmean :  If true, the confidence intervals corresponding to the point estimates evaluated at the mean of `x` within each binare presented. The default is `cigridmean=T`. ,

cb :  a vector. `cb=c(p,s)` sets a the piecewise polynomial of degree `p` with `s` smoothness constraints used forconstructing the confidence band. By default, the confidence band is not included in the plot unlessexplicitly specified.  Recommended specification is `cb=c(3,3)`, which adds a confidence bandbased on cubic B-spline estimate of the regression function of interest to the binned scatter plot. ,

cbgrid :  number of evaluation points of an evenly-spaced grid within each bin used for evaluation of the pointestimate set by the `cb=c(p,s)` option. The default is `cbgrid=20`, which correspondsto 20 evenly-spaced evaluation points within each bin for confidence interval construction. ,

polyreg :  degree of a global polynomial regression model for plotting. By default, this fit is not includedin the plot unless explicitly specified. Recommended specification is `polyreg=3`, whichadds a cubic (global) polynomial fit of the regression function of interest to the binned scatter plot. ,

polyreggrid :  number of evaluation points of an evenly-spaced grid within each bin used for evaluation ofthe point estimate set by the `polyreg=p` option. The default is `polyreggrid=20`,which corresponds to 20 evenly-spaced evaluation points within each bin for confidenceinterval construction. ,

polyregcigrid :  number of evaluation points of an evenly-spaced grid within each bin used for constructingconfidence intervals based on polynomial regression set by the `polyreg=p` option.The default is `polyregcigrid=0`, which corresponds to not plotting confidenceintervals for the global polynomial regression approximation. ,

by :  a vector containing the group indicator for subgroup analysis; both numeric and string variablesare supported. When `by` is specified, `binsreg` implements estimation and inference by each subgroupseparately, but produces a common binned scatter plot. By default, the binning structure is selected for eachsubgroup separately, but see the option `samebinsby` below for imposing a common binning structure across subgroups. ,

bycolors :  an ordered list of colors for plotting each subgroup series defined by the option `by`. ,

bysymbols :  an ordered list of symbols for plotting each subgroup series defined by the option `by`. ,

bylpatterns :  an ordered list of line patterns for plotting each subgroup series defined by the option `by`. ,

legendTitle :  String, title of legend. ,

legendoff :  If true, no legend is added. ,

testmodel :  a vector. `testmodel=c(p,s)` sets a piecewise polynomial of degree `p` with `s`smoothness constraints for parametric model specification testing.  The default is`testmodel=c(3,3)`, which corresponds to a cubic B-spline estimate of the regressionfunction of interest for testing against the fitting from a parametric model specification. ,

testmodelparfit :  a data frame or matrix which contains the evaluation grid and fitted values of the model(s) to betested against.  The first column contains a series of evaluation pointsat which the binscatter model and the parametric model of interest are compared witheach other.  Each parametric model is represented by other columns, which mustcontain the fitted values at the corresponding evaluation points. ,

testmodelpoly :  degree of a global polynomial model to be tested against. ,

testshape :  a vector. `testshape=c(p,s)` sets a piecewise polynomial of degree `p` with `s`smoothness constraints for nonparametric shape restriction testing. The default is`testshape=c(3,3)`, which corresponds to a cubic B-spline estimate of the regressionfunction of interest for one-sided or two-sided testing. ,

testshapel :  a vector of null boundary values for hypothesis testing. Each number `a` in the vectorcorresponds to one boundary of a one-sided hypothesis test to the left of the form`H0: sup_x mu(x)<=a`. ,

testshaper :  a vector of null boundary values for hypothesis testing. Each number `a` in the vectorcorresponds to one boundary of a one-sided hypothesis test to the right of the form`H0: inf_x mu(x)>=a`. ,

testshape2 :  a vector of null boundary values for hypothesis testing. Each number `a` in the vectorcorresponds to one boundary of a two-sided hypothesis test ofthe form`H0: sup_x |mu(x)-a|=0`. ,

nbins :  number of bins for partitioning/binning of `x`.  If not specified, the number of bins isselected via the companion function `binsregselect` in a data-driven, optimal way whenever possible. ,

binspos :  position of binning knots. The default is `binspos="qs"`, which corresponds to quantile-spacedbinning (canonical binscatter).  The other options are `"es"` for evenly-spaced binning, ora vector for manual specification of the positions of inner knots (which must be within the range of`x`). ,

binsmethod :  method for data-driven selection of the number of bins. The default is `binsmethod="dpi"`,which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: `"rot"`for rule of thumb implementation. ,

nbinsrot :  initial number of bins value used to construct the DPI number of bins selector.If not specified, the data-driven ROT selector is used instead. ,

samebinsby :  if true, a common partitioning/binning structure across all subgroups specified by the option `by` is forced.The knots positions are selected according to the option `binspos` and using the full sample. If `nbins`is not specified, then the number of bins is selected via the companion command `binsregselect` andusing the full sample. ,

nsims :  number of random draws for constructing confidence bands and hypothesis testing. The default is`nsims=500`, which corresponds to 500 draws from a standard Gaussian random vector of size`[(p+1)*J - (J-1)*s]`. ,

simsgrid :  number of evaluation points of an evenly-spaced grid within each bin used for evaluation ofthe supremum (or infimum) operation needed to construct confidence bands and hypothesis testingprocedures. The default is `simsgrid=20`, which corresponds to 20 evenly-spacedevaluation points within each bin for approximating the supremum (or infimum) operator. ,

simsseed :  seed for simulation. ,

vce :  Procedure to compute the variance-covariance matrix estimator. Options are\itemize{\item `"const"` homoskedastic variance estimator.\item `"HC0"` heteroskedasticity-robust plug-in residuals variance estimator                   without weights.\item `"HC1"` heteroskedasticity-robust plug-in residuals variance estimator                   with hc1 weights. Default.\item `"HC2"` heteroskedasticity-robust plug-in residuals variance estimator                   with hc2 weights.\item `"HC3"` heteroskedasticity-robust plug-in residuals variance estimator                   with hc3 weights.} ,

cluster :  cluster ID. Used for compute cluster-robust standard errors. ,

level :  nominal confidence level for confidence interval and confidence band estimation. Default is `level=95`. ,

noplot :  If true, no plot produced. ,

dfcheck :  adjustments for minimum effective sample size checks, which take into account number of uniquevalues of `x` (i.e., number of mass points), number of clusters, and degrees of freedom ofthe different stat models considered. The default is `dfcheck=c(20, 30)`.See \href{{https://sites.google.com/site/nppackages/binsreg/Cattaneo-Crump-Farrell-Feng_2019_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2019b)}} for more details. ,

masspoints :  how mass points in `x` are handled. Available options:\itemize{\item `"on"` all mass point and degrees of freedom checks are implemented. Default.\item `"noadjust"` mass point checks and the corresponding effective sample size adjustments are omitted.\item `"nolocalcheck"` within-bin mass point and degrees of freedom checks are omitted.\item `"off"` "noadjust" and "nolocalcheck" are set simultaneously.\item `"veryfew"` forces the function to proceed as if `x` has only a few number of mass points (i.e., distinct values).                       In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.} ,

weights :  an optional vector of weights to be used in the fitting process. Should be `NULL` ora numeric vector. For more details, see `lm`. ,

subset :  Optional rule specifying a subset of observations to be used. ,
```



# Author

Jiafeng Chen

