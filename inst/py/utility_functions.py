# (C) 2019 Potsdam Institute for Climate Impact Research (PIK)
# 
# This file is part of ISIMIP3BASD.
#
# ISIMIP3BASD is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ISIMIP3BASD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with ISIMIP3BASD. If not, see <http://www.gnu.org/licenses/>.



"""
Utility functions
=================

Provides auxiliary functions used by the modules bias_adjustment and
statistical_downscaling.

"""



import os
import collections
import warnings
import iris
import pickle
import numpy as np
import pandas as pd
import datetime as dt
import scipy.stats as sps
import scipy.linalg as spl
import scipy.interpolate as spi
from scipy.signal import convolve



def assert_validity_of_months(months):
    """
    Raises an assertion error if any of the numbers in months is not in
    {1,...,12}.

    Parameters
    ----------
    months : array_like
        Sequence of ints representing calendar months.

    """
    months_allowed = np.arange(1, 13)
    for month in months:
        assert month in months_allowed, f'found {month} in months'



def assert_coord_axis(cube, coord='time', axis=0):
    """
    Raises an assertion error if the axis of the specified coordinate in the
    specified iris cube differs from the specified one.

    Parameters
    ----------
    cube : iris cube
        Cube for which the assertion shall be done.
    coord : str, optional
        Coordinate identifier.
    axis : int, optional
        Dimension number.

    """
    msg = coord+' must be the coordinate of axis %i'%axis
    assert cube.coord_dims(coord) == (axis,), msg



def assert_consistency_of_bounds_and_thresholds(
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None):
    """
    Raises an assertion error if the pattern of specified and
    unspecified bounds and thresholds is not valid or if
    lower_bound < lower_threshold < upper_threshold < upper_bound
    does not hold.

    Parameters
    ----------
    lower_bound : float, optional
        Lower bound of values in time series.
    lower_threshold : float, optional
        Lower threshold of values in time series. All values below this
        threshold will be replaced by random numbers between lower_bound and
        lower_threshold before bias adjustment.
    upper_bound : float, optional
        Upper bound of values in time series.
    upper_threshold : float, optional
        Upper threshold of values in time series. All values above this
        threshold will be replaced by random numbers between upper_threshold and
        upper_bound before bias adjustment.

    """
    lower = lower_bound is not None and lower_threshold is not None
    upper = upper_bound is not None and upper_threshold is not None

    if not lower:
       msg = 'lower_bound is not None and lower_threshold is None'
       assert lower_bound is None, msg
       msg = 'lower_bound is None and lower_threshold is not None'
       assert lower_threshold is None, msg
    if not upper:
       msg = 'upper_bound is not None and upper_threshold is None'
       assert upper_bound is None, msg
       msg = 'upper_bound is None and upper_threshold is not None'
       assert upper_threshold is None, msg

    if lower:
        assert lower_bound < lower_threshold, 'lower_bound >= lower_threshold'
    if upper:
        assert upper_bound > upper_threshold, 'upper_bound <= upper_threshold'
    if lower and upper:
        msg = 'lower_threshold >= upper_threshold'
        assert lower_threshold < upper_threshold, msg



def assert_consistency_of_distribution_and_bounds(
        distribution,
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None):
    """
    Raises an assertion error if the the distribution is not consistent with the
    pattern of specified and unspecified bounds and thresholds.

    Parameters
    ----------
    distribution : str
        Kind of distribution used for parametric quantile mapping:
        [None, 'normal', 'weibull', 'gamma', 'beta', 'rice'].
    lower_bound : float, optional
        Lower bound of values in time series.
    lower_threshold : float, optional
        Lower threshold of values in time series. All values below this
        threshold will be replaced by random numbers between lower_bound and
        lower_threshold before bias adjustment.
    upper_bound : float, optional
        Upper bound of values in time series.
    upper_threshold : float, optional
        Upper threshold of values in time series. All values above this
        threshold will be replaced by random numbers between upper_threshold and
        upper_bound before bias adjustment.

    """
    if distribution is not None:
        lower = lower_bound is not None and lower_threshold is not None
        upper = upper_bound is not None and upper_threshold is not None
    
        msg = distribution+' distribution '
        if distribution == 'normal':
            assert not lower and not upper, msg+'can not have bounds'
        elif distribution in ['weibull', 'gamma', 'rice']:
            assert lower and not upper, msg+'must only have lower bound'
        elif distribution == 'beta':
            assert lower and upper, msg+'must have lower and upper bound'
        else:
            raise AssertionError(msg+'not supported')



def assert_calendar(cube, calendar='proleptic_gregorian'):
    """
    Raises an assertion error if the calendar of the given iris cube differs
    from the one given.

    Parameters
    ----------
    cube : iris cube
        Cube for which the assertion shall be done.
    calendar : string, optional
        Calendar identifier.

    """
    msg = 'cube calendar != '+calendar
    assert cube.coord('time').units.calendar == calendar, msg



def assert_no_infs_or_nans(x_before, x_after):
    """
    Raises a value error if there are infs or nans in x_after. Prints the
    corresponding values in x_before.

    Parameters
    ----------
    x_before : ndarray
        Array before bias adjustement or statistical downscaling.
    x_after : ndarray
        Array after bias adjustement or statistical downscaling.

    """
    is_invalid = np.logical_or(np.isinf(x_after), np.isnan(x_after))
    if np.any(is_invalid):
        print(x_before[is_invalid])
        print(x_after[is_invalid], flush=True)
        msg = 'found infs or nans in x_after'
        raise ValueError(msg)



def split(s, n=None, converter=str, empty=None, delimiter=','):
    """
    Splits string s into a list of strings using the specified delimiter,
    replaces empty strings in the resulting list by empty, and applies the
    specified converter to all other values in the resulting list.

    If the split of s results in a list of length 1 and n > 1 then s is repeated
    to obtain a list of length n. Otherwise the split of s is asserted to result
    in a list of length n.

    Parameters
    ----------
    s : str
        String to be split.
    n : int, optional
        Target list length.
    converter : function, optional
        Function to change the data type of the list elements
    empty : any type, optional
        Value that empty list elements are mapped to.
    delimiter : str, optional
        Delimiter used to split s.

    Returns
    -------
    l : list
        Resulting list.

    """
    l = s.split(delimiter)
    if isinstance(n, int):
        m = len(l)
        if m == 1 and n > 1:
            l = [s for i in range(n)]
        else:
            msg = 'list length = %i != %i'%(m, n)
            assert m == n, msg
    return [empty if s == '' else converter(s) for s in l]



def add_basd_attributes(cube, options, prefix='', index=None):
    """
    Adds keys and values of all command line options as attributes to iris cube.

    Parameters
    ----------
    cube : iris cube
        Cube to which attributes shall be added. Is changed in-place.
    options : optparse.Values
        Command line options parsed by optparse.OptionParser.
    prefix : str, optional
        Prefix attached to all attributes before they are put into the cube.
    index : int, optional
        If provided then only save the options relevant for the variable with
        the given index.

    """
    a = cube.attributes
    a[prefix+'version'] = 'ISIMIP3BASD v2.4.1'
    for key, value in options.__dict__.items():
        if index is not None and key != 'months' and isinstance(value, str):
            v = value.split(',')[index] if ',' in value else value
        else:
            v = value
        a[prefix+key] = str(v)
    cube.attributes = a



def aggregate_periodic(a, halfwin, aggregator='mean'):
    """
    Aggregates a using the given aggregator and a running window of length
    2 * halfwin + 1 assuming that a is periodic.

    Parameters
    ----------
    a : array
        Array to be aggregated.
    halfwin : int
        Determines length of running window used for aggregation.
    aggregator : str, optional
        Determines how a is aggregated along axis 0 for every running window.

    Returns
    -------
    rm : ndarray
        Result of aggregation. Same shape as a.

    """
    assert halfwin >= 0, 'halfwin < 0'
    if not halfwin: return a

    # extend a periodically
    n = a.size
    assert n >= halfwin, 'length of a along axis 0 less than halfwin'
    b = np.concatenate((a[-halfwin:], a, a[:halfwin]))

    # aggregate using algorithm for max inspired by
    # <http://p-nand-q.com/python/algorithms/searching/max-sliding-window.html>
    window = 2 * halfwin + 1
    if aggregator == 'max':
        c = list(np.maximum.accumulate(b[:window][::-1]))
        rm = np.empty_like(a)
        rm[0] = c[-1]
        for i in range(n-1):
            c_new = b[i+window]
            del c[-1]
            for j in range(window-1):
                if c_new > c[j]: c[j] = c_new
                else: break
            c.insert(0, c_new)
            rm[i+1] = c[-1]
    elif aggregator == 'mean':
        rm = convolve(b, np.repeat(1./window, window), 'valid')
    else:
        raise ValueError(f'aggregator {aggregator} not supported')

    return rm



def get_upper_bound_climatology(d, doys, halfwin):
    """
    Estimates an annual cycle of upper bounds as running mean values of running
    maximum values of multi-year daily maximum values.

    Parameters
    ----------
    d : array
        Time series for which annual cycles of upper bounds shall be estimated.
    doys : array
        Day of the year time series corresponding to d.
    halfwin : int
        Determines length of running windows used for estimation.

    Returns
    -------
    ubc : array
        Upper bound climatology.
    doys_unique : array
        Days of the year of upper bound climatology.

    """
    assert d.shape == doys.shape, 'd and doys differ in shape' 

    # check length of time axis of resulting array
    doys_unique, counts = np.unique(doys, return_counts=True)
    n = doys_unique.size
    if n != 366:
        msg = (f'upper bound climatology only defined for {n} days of the year:'
                ' this may imply an invalid computation of the climatology')
        warnings.warn(msg)

    # compute multi year daily maximum
    d_sorted = d[np.argsort(doys)]
    mydm = np.empty(n, dtype=d.dtype)
    if np.unique(counts[:-1]).size == 1:
        # fast version which applies in the usual case
        if counts[0] == counts[-1]:
            d_stacked = d_sorted.reshape(n, counts[0])
            mydm = np.max(d_stacked, axis=1) 
        else:
            mydm[-1] =  np.max(d_sorted[-counts[-1]:])
            d_stacked = d_sorted[:-counts[-1]].reshape(n-1, counts[0])
            mydm[:-1] = np.max(d_stacked, axis=1) 
    else:
        # slow version which always works
        j = 0
        for i in range(n):
            k = j + counts[i]
            mydm[i] = np.max(d_sorted[j:k])
            j = k

    # smooth multi year daily maximum
    mydmrm = aggregate_periodic(mydm, halfwin, aggregator='max')
    ubc = aggregate_periodic(mydmrm, halfwin, aggregator='mean')

    return ubc, doys_unique



def ccs_transfer_sim2obs_upper_bound_climatology(obs_hist, sim_hist, sim_fut):
    """
    Multiplicatively transfers simulated climate change signal from sim_hist,
    sim_fut to obs_hist.

    Parameters
    ----------
    obs_hist : array
        Upper bound climatology of observed climate data representing the
        historical or training time period.
    sim_hist : array
        Upper bound climatology of simulated climate data representing the
        historical or training time period.
    sim_fut : array
        Upper bound climatology of simulated climate data representing the
        future or application time period.

    Returns
    -------
    sim_fut_ba : array
        Result of climate change signal transfer.

    """
    assert obs_hist.shape == sim_hist.shape == sim_fut.shape, \
        'obs_hist, sim_hist, sim_fut differ in shape'
    with np.errstate(divide='ignore', invalid='ignore'):
        change_factor = np.where(sim_hist == 0, 1, sim_fut / sim_hist)
    change_factor = np.maximum(.1, np.minimum(10., change_factor))
    sim_fut_ba = obs_hist * change_factor
    return sim_fut_ba



def scale_by_upper_bound_climatology(
        d, ubc, d_doys, ubc_doys, divide=True):
    """
    Scales all values in d using the annual cycle of upper bounds.

    Parameters
    ----------
    d : array
        Time series to be scaled. Is changed in-place.
    ubc : array
        Upper bound climatology used for scaling.
    d_doys : array
        Days of the year corresponding to d.
    ubc_doys : array
        Days of the year corresponding to ubc.
    divide : boolean, optional
        If True the cube is divided by upper_bound_climatology, otherwise they
        are multiplied.

    """
    assert d.shape == d_doys.shape, 'd and d_doys differ in shape' 
    assert ubc.shape == ubc_doys.shape, 'ubc and ubc_doys differ in shape' 

    if divide:
        with np.errstate(divide='ignore', invalid='ignore'):
            scaling_factors = np.where(ubc == 0, 1., 1./ubc)
    else:
        scaling_factors = ubc

    # use fast solution if ubc covers all days of the year
    # this fast solution assumes that ubc_doys is sorted
    scaling_factors_broadcasted = scaling_factors[d_doys-1] \
        if ubc.size == 366 else \
        np.array([scaling_factors[ubc_doys == doy][0] for doy in d_doys])

    d *= scaling_factors_broadcasted

    # make sure d does not exceed the upper bound climatology
    if not divide:
        d_too_large = d > scaling_factors_broadcasted
        n = np.sum(d_too_large)
        msg = f'capping {n} values exceeding the upper bound climatology'
        if n:
            warnings.warn(msg)
            d[d_too_large] = scaling_factors_broadcasted[d_too_large]



def subtract_or_add_trend(x, years, trend=None):
    """
    Subtracts or adds trend from or to x.

    Parameters
    ----------
    x : array
        Time series.
    years : array
        Years of time points of x used to subtract or add trend at annual
        temporal resolution.
    trend : array, optional
        Trend line. If provided then this is the trend line added to x.
        Otherwise, a trend line is computed and subtracted from x

    Returns
    -------
    y : array
        Result of trend subtraction or addition from or to x.
    trend : array, optional
        Trend line. Is only returned if the parameter trend is None.

    """
    assert x.size == years.size, 'size of x != size of years'
    unique_years = np.unique(years)

    # compute trend
    if trend is None:
        annual_means = np.array([np.mean(x[years == y]) for y in unique_years])
        r = sps.linregress(unique_years, annual_means)
        if r.pvalue < .05:  # detrend preserving multi-year mean value
            trend = r.slope * (unique_years - np.mean(unique_years))
        else:  # do not detrend because trend is insignificant
            trend = np.zeros(unique_years.size, dtype=x.dtype)
        return_trend = True
    else:
        msg = 'size of trend array != number of unique years'
        assert trend.size == unique_years.size, msg
        trend = -trend
        return_trend = False

    # subtract or add trend
    if np.any(trend):
        y = np.empty_like(x)
        for i, year in enumerate(unique_years):
            ## is_year = years == year
            is_year = np.reshape(years == year, y.shape)
            y[is_year] = x[is_year] - trend[i]
    else:
        y = x.copy()

    # return result(s)
    if return_trend:
        return y, trend
    else:
        return y



def percentile1d(a, p):
    """
    Fast version of np.percentile with linear interpolation for 1d arrays
    inspired by
    <https://krstn.eu/np.nanpercentile()-there-has-to-be-a-faster-way/>.

    Parameters
    ----------
    a : array
        Input array.
    p : array
        Percentages expressed as real numbers in [0, 1] for which percentiles
        are computed.

    Returns
    -------
    percentiles : array
        Percentiles

    """
    n = a.size - 1
    b = np.sort(a)
    i = n * p
    i_below = np.floor(i).astype(int)
    w_above = i - i_below
    return b[i_below] * (1. - w_above) + b[i_below + (i_below < n)] * w_above



def map_quantiles_non_parametric_trend_preserving(
        x_obs_hist, x_sim_hist, x_sim_fut, 
        trend_preservation='additive', n_quantiles=50,
        max_change_factor=100., max_adjustment_factor=9.,
        adjust_obs=False, lower_bound=None, upper_bound=None):
    """
    Adjusts biases with a modified version of the quantile delta mapping by
    Cannon (2015) <https://doi.org/10.1175/JCLI-D-14-00754.1> or uses this
    method to transfer a simulated climate change signal to observations.

    Parameters
    ----------
    x_obs_hist : array
        Time series of observed climate data representing the historical or
        training time period.
    x_sim_hist : array
        Time series of simulated climate data representing the historical or
        training time period.
    x_sim_fut : array
        Time series of simulated climate data representing the future or
        application time period.
    trend_preservation : str, optional
        Kind of trend preservation:
        'additive'       # Preserve additive trend.
        'multiplicative' # Preserve multiplicative trend, ensuring
                         # 1/max_change_factor <= change factor
                         #                     <= max_change_factor.
        'mixed'          # Preserve multiplicative or additive trend or mix of
                         # both depending on sign and magnitude of bias. Purely
                         # additive trends are preserved if adjustment factors
                         # of a multiplicative adjustment would be greater then
                         # max_adjustment_factor.
        'bounded'        # Preserve trend of bounded variable. Requires
                         # specification of lower_bound and upper_bound. It is
                         # ensured that the resulting values stay within these
                         # bounds.
    n_quantiles : int, optional
        Number of quantile-quantile pairs used for non-parametric quantile
        mapping.
    max_change_factor : float, optional
        Maximum change factor applied in non-parametric quantile mapping with
        multiplicative or mixed trend preservation.
    max_adjustment_factor : float, optional
        Maximum adjustment factor applied in non-parametric quantile mapping
        with mixed trend preservation.
    adjust_obs : boolean, optional
        If True then transfer simulated climate change signal to x_obs_hist,
        otherwise apply non-parametric quantile mapping to x_sim_fut.
    lower_bound : float, optional
        Lower bound of values in x_obs_hist, x_sim_hist, and x_sim_fut. Used
        for bounded trend preservation.
    upper_bound : float, optional
        Upper bound of values in x_obs_hist, x_sim_hist, and x_sim_fut. Used
        for bounded trend preservation.

    Returns
    -------
    y : array
        Result of quantile mapping or climate change signal transfer.

    """
    # make sure there are enough input data for quantile delta mapping
    # reduce n_quantiles if necessary
    assert n_quantiles > 0, 'n_quantiles <= 0'
    n = min([n_quantiles + 1, x_obs_hist.size, x_sim_hist.size, x_sim_fut.size])
    if n < 2:
        msg = 'not enough input data: returning x_obs_hist'
        warnings.warn(msg)
        return x_obs_hist
    elif n < n_quantiles + 1:
        msg = 'due to little input data: reducing n_quantiles to %i'%(n-1)
        warnings.warn(msg)
    p_zeroone = np.linspace(0., 1., n)

    # compute quantiles of input data
    q_obs_hist = percentile1d(x_obs_hist, p_zeroone)
    q_sim_hist = percentile1d(x_sim_hist, p_zeroone)
    q_sim_fut = percentile1d(x_sim_fut, p_zeroone)

    # compute quantiles needed for quantile delta mapping
    if adjust_obs: p = np.interp(x_obs_hist, q_obs_hist, p_zeroone)
    else: p = np.interp(x_sim_fut, q_sim_fut, p_zeroone)
    F_sim_fut_inv  = np.interp(p, p_zeroone, q_sim_fut)
    F_sim_hist_inv = np.interp(p, p_zeroone, q_sim_hist)
    F_obs_hist_inv = np.interp(p, p_zeroone, q_obs_hist)

    # do augmented quantile delta mapping
    if trend_preservation == 'bounded':
        msg = 'lower_bound or upper_bound not specified'
        assert lower_bound is not None and upper_bound is not None, msg
        assert lower_bound < upper_bound, 'lower_bound >= upper_bound'
        y = ccs_transfer_sim2obs(
            F_obs_hist_inv, F_sim_hist_inv, F_sim_fut_inv,
            lower_bound, upper_bound)
    elif trend_preservation in ['mixed', 'multiplicative']:
        assert max_change_factor > 1, 'max_change_factor <= 1'
        with np.errstate(divide='ignore', invalid='ignore'):
            y = np.where(F_sim_hist_inv == 0, 1., F_sim_fut_inv/F_sim_hist_inv)
            y[y > max_change_factor] = max_change_factor
            y[y < 1. / max_change_factor] = 1. / max_change_factor
        y *= F_obs_hist_inv
        if trend_preservation == 'mixed':  # if not then we are done here
            assert max_adjustment_factor > 1, 'max_adjustment_factor <= 1'
            y_additive = F_obs_hist_inv + F_sim_fut_inv - F_sim_hist_inv
            fraction_multiplicative = np.zeros_like(y)
            fraction_multiplicative[F_sim_hist_inv >= F_obs_hist_inv] = 1.
            i_transition = np.logical_and(F_sim_hist_inv < F_obs_hist_inv,
                F_obs_hist_inv < max_adjustment_factor * F_sim_hist_inv)
            fraction_multiplicative[i_transition] = .5 * (1. + 
                np.cos((F_obs_hist_inv[i_transition] /
                F_sim_hist_inv[i_transition] - 1.) * 
                np.pi / (max_adjustment_factor - 1.)))
            y = fraction_multiplicative * y + (1. -
                fraction_multiplicative) * y_additive
    elif trend_preservation == 'additive':
        y = F_obs_hist_inv + F_sim_fut_inv - F_sim_hist_inv
    else:
        msg = 'trend_preservation = '+trend_preservation+' not supported'
        raise AssertionError(msg)

    return y



def map_quantiles_non_parametric_with_constant_extrapolation(x, q_sim, q_obs):
    """
    Uses quantile-quantile pairs represented by values in q_sim and q_obs
    for quantile mapping of x.

    Values in x beyond the range of q_sim are mapped following the constant
    extrapolation approach, see Boe et al. (2007)
    <https://doi.org/10.1002/joc.1602>.

    Parameters
    ----------
    x : array
        Simulated time series.
    q_sim : array
        Simulated quantiles.
    q_obs : array
        Observed quantiles.
    
    Returns
    -------
    y : array
        Result of quantile mapping.

    """
    assert q_sim.size == q_obs.size
    lunder = x < q_sim[0]
    lover = x > q_sim[-1]
    y = np.interp(x, q_sim, q_obs)
    y[lunder] = x[lunder] + (q_obs[0] - q_sim[0])
    y[lover] = x[lover] + (q_obs[-1] - q_sim[-1])
    return y



def map_quantiles_non_parametric_brute_force(x, y):
    """
    Quantile-map x to y using the empirical CDFs of x and y.

    Parameters
    ----------
    x : array
        Simulated time series.
    y : array
        Observed time series.
    
    Returns
    -------
    z : array
        Result of quantile mapping.

    """
    if x.size == 0:
        msg = 'found no values in x: returning x'
        warnings.warn(msg)
        return x

    if np.unique(y).size < 2:
        msg = 'found fewer then 2 different values in y: returning x'
        warnings.warn(msg)
        return x

    p_x = (sps.rankdata(x) - 1.) / x.size  # percent points of x
    p_y = np.linspace(0., 1., y.size)  # percent points of sorted y
    z = np.interp(p_x, p_y, np.sort(y))  # quantile mapping
    return z



def ccs_transfer_sim2obs(
        x_obs_hist, x_sim_hist, x_sim_fut,
        lower_bound=0., upper_bound=1.):
    """
    Generates pseudo future observation(s) by transfering a simulated climate
    change signal to historical observation(s) respecting the given bounds.

    Parameters
    ----------
    x_obs_hist : float or array
        Historical observation(s).
    x_sim_hist : float or array
        Historical simulation(s).
    x_sim_fut : float or array
        Future simulation(s).
    lower_bound : float, optional
        Lower bound of values in input and output data.
    upper_bound : float, optional
        Upper bound of values in input and output data.
    
    Returns
    -------
    x_obs_fut : float or array
        Pseudo future observation(s).

    """
    # change scalar inputs to arrays
    if np.isscalar(x_obs_hist): x_obs_hist = np.array([x_obs_hist])
    if np.isscalar(x_sim_hist): x_sim_hist = np.array([x_sim_hist])
    if np.isscalar(x_sim_fut): x_sim_fut = np.array([x_sim_fut])

    # check input
    assert lower_bound < upper_bound, 'lower_bound >= upper_bound'
    for x_name, x in zip(['x_obs_hist', 'x_sim_hist', 'x_sim_fut'],
                         [x_obs_hist, x_sim_hist, x_sim_fut]):
        assert np.all(x >= lower_bound), 'found '+x_name+' < lower_bound'
        assert np.all(x <= upper_bound), 'found '+x_name+' > upper_bound'

    # compute x_obs_fut
    i_neg_bias = x_sim_hist < x_obs_hist
    i_zero_bias = x_sim_hist == x_obs_hist
    i_pos_bias = x_sim_hist > x_obs_hist
    i_additive = np.logical_or(
        np.logical_and(i_neg_bias, x_sim_fut < x_sim_hist),
        np.logical_and(i_pos_bias, x_sim_fut > x_sim_hist))
    x_obs_fut = np.empty_like(x_obs_hist)
    x_obs_fut[i_neg_bias] = upper_bound - \
                            (upper_bound - x_obs_hist[i_neg_bias]) * \
                            (upper_bound - x_sim_fut[i_neg_bias]) / \
                            (upper_bound - x_sim_hist[i_neg_bias])
    x_obs_fut[i_zero_bias] = x_sim_fut[i_zero_bias]
    x_obs_fut[i_pos_bias] = lower_bound + \
                            (x_obs_hist[i_pos_bias] - lower_bound) * \
                            (x_sim_fut[i_pos_bias] - lower_bound) / \
                            (x_sim_hist[i_pos_bias] - lower_bound)
    x_obs_fut[i_additive] = x_obs_hist[i_additive] + \
                            x_sim_fut[i_additive] - x_sim_hist[i_additive]

    # make sure x_obs_fut is within bounds
    x_obs_fut = np.maximum(lower_bound, np.minimum(upper_bound, x_obs_fut))

    return x_obs_fut[0] if x_obs_fut.size == 1 else x_obs_fut



def transfer_odds_ratio(p_obs_hist, p_sim_hist, p_sim_fut):
    """
    Transfers simulated changes in event likelihood to historical observations
    by multiplying the historical odds by simulated future-over-historical odds
    ratio. The method is inspired by the return interval scaling proposed by
    Switanek et al. (2017) <https://doi.org/10.5194/hess-21-2649-2017>.

    Parameters
    ----------
    p_obs_hist : array
        Culmulative probabbilities of historical observations.
    p_sim_hist : array
        Culmulative probabbilities of historical simulations.
    p_sim_fut : array
        Culmulative probabbilities of future simulations.

    Returns
    -------
    p_obs_fut : array
        Culmulative probabbilities of pseudo future observations.

    """
    x = np.sort(p_obs_hist)
    y = np.sort(p_sim_hist)
    z = np.sort(p_sim_fut)

    # interpolate x and y if necessary
    if x.size != z.size or y.size != z.size:
        p_x = np.linspace(0, 1, x.size)
        p_y = np.linspace(0, 1, y.size)
        p_z = np.linspace(0, 1, z.size)
        ppf_x = spi.interp1d(p_x, x)
        ppf_y = spi.interp1d(p_y, y)
        x = ppf_x(p_z)
        y = ppf_y(p_z)

    # transfer
    A = x * (1. - y) * z
    B = (1. - x) * y * (1. - z)
    z_scaled = 1. / (1. + B / A)

    # avoid the generation of unrealistically extreme p-values
    z_min = 1. / (1. + np.power(10.,  1. - np.log10(x / (1. - x))))
    z_max = 1. / (1. + np.power(10., -1. - np.log10(x / (1. - x))))
    z_scaled = np.maximum(z_min, np.minimum(z_max, z_scaled))

    return z_scaled[np.argsort(np.argsort(p_sim_fut))]



def randomize_censored_values_core(y, bound, threshold, inverse, power, lower):
    """
    Randomizes values beyond threshold in y or de-randomizes such formerly
    randomized values. Note that y is changed in-place. The randomization
    algorithm is inspired by <https://stackoverflow.com/questions/47429845/
    rank-with-ties-in-python-when-tie-breaker-is-random>

    Parameters
    ----------
    y : array
        Time series to be (de-)randomized.
    bound : float
        Lower or upper bound of values in time series.
    threshold : float
        Lower or upper threshold of values in time series.
    inverse : boolean
        If True, values beyond threshold in y are set to bound.
        If False, values beyond threshold in y are randomized.
    power : float
        Numbers for randomizing values are drawn from a uniform distribution
        and then taken to this power.
    lower : boolean
        If True/False, consider bound and threshold to be lower/upper bound and
        lower/upper threshold, respectively.

    """
    if lower: i = y < threshold
    else: i = y > threshold
    if inverse:
        y[i] = bound
    else:
        n = np.sum(i)
        if n:
            p = np.power(np.random.uniform(0, 1, n), power)
            v = bound + p * (threshold - bound)
            s = pd.Series(y[i])
            r = s.sample(frac=1).rank(method='first').reindex_like(s)
            y[i] = np.sort(v)[r.values.astype(int) - 1]



def randomize_censored_values(x,
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None,
        inplace=False, inverse=False,
        seed=None, lower_power=1., upper_power=1.):
    """
    Randomizes values beyond threshold in x or de-randomizes such formerly
    randomized values.

    Parameters
    ----------
    x : array
        Time series to be (de-)randomized.
    lower_bound : float, optional
        Lower bound of values in time series.
    lower_threshold : float, optional
        Lower threshold of values in time series.
    upper_bound : float, optional
        Upper bound of values in time series.
    upper_threshold : float, optional
        Upper threshold of values in time series.
    inplace : boolean, optional
        If True, change x in-place. If False, change a copy of x.
    inverse : boolean, optional
        If True, values beyond thresholds in x are set to the respective bound.
        If False, values beyond thresholds in x are randomized, i.e. values that
        exceed upper_threshold are replaced by random numbers from the
        interval [lower_bound, lower_threshold), and values that fall short
        of lower_threshold are replaced by random numbers from the interval
        (upper_threshold, upper_bound]. The ranks of the censored values are
        preserved using a random tie breaker. 
    seed : int, optional
        Used to seed the random number generator before replacing values beyond
        threshold.
    lower_power : float, optional
        Numbers for randomizing values that fall short of lower_threshold are
        drawn from a uniform distribution and then taken to this power.
    upper_power : float, optional
        Numbers for randomizing values that exceed upper_threshold are drawn
        from a uniform distribution and then taken to this power.

    Returns
    -------
    x : array
        Randomized or de-randomized time series.

    """
    y = x if inplace else x.copy()
    if seed is not None:
        np.random.seed(seed)

    # randomize lower values
    if lower_bound is not None and lower_threshold is not None:
        randomize_censored_values_core(
            y, lower_bound, lower_threshold, inverse, lower_power, True)

    # randomize upper values
    if upper_bound is not None and upper_threshold is not None:
        randomize_censored_values_core(
            y, upper_bound, upper_threshold, inverse, upper_power, False)

    return y



def check_shape_loc_scale(spsdotwhat, shape_loc_scale):
    """
    Analyzes how distribution fitting has worked.

    Parameters
    ----------
    spsdotwhat : sps distribution class
        Known classes are [sps.norm, sps.weibull_min, sps.gamma, sps.rice,
        sps.beta].
    shape_loc_scale : tuple
        Fitted shape, location, and scale parameter values.

    Returns
    -------
    i : int
        0 if everything is fine,
        1 if there are infs or nans in shape_loc_scale,
        2 if at least one value in shape_loc_scale is out of bounds,
        3 if spsdotwhat is unknown.

    """
    if np.any(np.isnan(shape_loc_scale)) or np.any(np.isinf(shape_loc_scale)):
        return 1
    elif spsdotwhat == sps.norm:
        return 2 if shape_loc_scale[1] <= 0 else 0
    elif spsdotwhat in [sps.weibull_min, sps.gamma, sps.rice]:
        return 2 if shape_loc_scale[0] <= 0 or shape_loc_scale[2] <= 0 else 0
    elif spsdotwhat == sps.beta:
        return 2 if shape_loc_scale[0] <= 0 or shape_loc_scale[1] <= 0 \
            or shape_loc_scale[0] > 1e10 or shape_loc_scale[1] > 1e10 else 0
    else:
        return 3



def fit(spsdotwhat, x, fwords):
    """
    Attempts to fit a distribution from the family defined through spsdotwhat
    to the data represented by x, holding parameters fixed according to fwords.

    A maximum likelihood estimation of distribution parameter values is tried
    first. If that fails the method of moments is tried for some distributions.

    Parameters
    ----------
    spsdotwhat : sps distribution class
        Known classes are [sps.norm, sps.weibull_min, sps.gamma, sps.rice,
        sps.beta].
    x : array
        Data to be fitted.
    fwords : dict
        Dictionary with keys 'floc' and (optinally) 'fscale' specifying location
        and scale parmeter values, respectively, that are to be held fixed when
        fitting.

    Returns
    -------
    shape_loc_scale : tuple
        Fitted shape, location, and scale parameter values if fitting worked,
        otherwise None.

    """
    # make sure that there are at least two distinct data points because
    # otherwise it is impossible to fit more than 1 parameter
    if np.unique(x).size < 2:
        msg = 'found fewer then 2 different values in x: returning None'
        warnings.warn(msg)
        return None

    # try maximum likelihood estimation
    try:
        shape_loc_scale = spsdotwhat.fit(x, **fwords)
    except:
        shape_loc_scale = (np.nan,)


    # try method of moment estimation
    if check_shape_loc_scale(spsdotwhat, shape_loc_scale):
        msg = 'maximum likelihood estimation'
        if spsdotwhat == sps.gamma:
            msg += ' failed: method of moments'
            x_mean = np.mean(x) - fwords['floc']
            x_var = np.var(x)
            scale = x_var / x_mean
            shape = x_mean / scale
            shape_loc_scale = (shape, fwords['floc'], scale)
        elif spsdotwhat == sps.beta:
            msg += ' failed: method of moments'
            y = (x - fwords['floc']) / fwords['fscale']
            y_mean = np.mean(y)
            y_var = np.var(y)
            p = np.square(y_mean) * (1. - y_mean) / y_var - y_mean
            q = p * (1. - y_mean) / y_mean
            shape_loc_scale = (p, q, fwords['floc'], fwords['fscale'])
    else:
        msg = ''

    # return result and utter warning if necessary
    if check_shape_loc_scale(spsdotwhat, shape_loc_scale):
        msg += ' failed: returning None'
        warnings.warn(msg)
        return None
    elif msg != '':
        msg += ' succeeded'

    # do rough goodness of fit test to filter out worst fits using KS test
    ks_stat = sps.kstest(x, spsdotwhat.name, args=shape_loc_scale)[0]
    if ks_stat > .5:
        if msg == '': msg = 'maximum likelihood estimation succeeded'
        msg += ' but fit is not good: returning None'
        warnings.warn(msg)
        return None
    else:
        if msg != '':
            warnings.warn(msg)
        return shape_loc_scale



def output_chunksizes(output_shape, space_chunksize=10):
    """
    Returns chunk sizes for output netcdf file(s) with chunking in space. This
    means that (1) the chunk size of the 0th or time dimension is set equal to
    the length of the time dimension as specified by output_shape and (2) the
    chunk size of all other dimensions is set equal to the minimum of the length
    of the respective dimension as specified by output_shape and
    space_chunksize.

    Parameters
    ----------
    output_shape : tuple
        Shape of data array of output netcdf file(s).
    space_chunksize : int, optional
        Chunk size in spatial dimensions.

    Returns
    -------
    chunksizes : tuple
        Chunk sizes for output netcdf file(s).

    """
    space_shape = np.array(output_shape[1:])
    space_chunksizes = np.repeat(space_chunksize, space_shape.size)
    space_chunksizes = np.minimum(space_chunksizes, space_shape)
    chunksizes = output_shape[:1] + tuple(space_chunksizes)
    return chunksizes



def npy_stack_dir(path):
    """
    Turns path to output netcdf file into path to directory that contains the
    associated npy stack.

    Parameters
    ----------
    path : str
        Path to file.

    Returns
    -------
    result : str
        Path to directory.

    """
    return path + '.npy_stack/'



def setup_npy_stack(path, shape):
    """
    Prepares storage of local bias adjustment or statistical downscaling results
    such that these can be loaded using dask.array.from_npy_stack by creating a
    directory for the npy files and putting a suitable info file into that
    directory, see also dask.array.to_npy_stack.

    Parameters
    ----------
    path : str
        Path to associated output netcdf file.
    shape : tuple
        Shape of data array of output netcdf file.

    """
    chunks = (shape[:1], tuple(np.ones(np.prod(np.array(shape[1:])), int)))
    npy_stack_meta = {'chunks': chunks, 'dtype': np.float32, 'axis': 1}
    if not os.path.exists(npy_stack_dir(path)):
        os.makedirs(npy_stack_dir(path))
    with open(npy_stack_dir(path)+'info', 'wb') as f:
        pickle.dump(npy_stack_meta, f)



def get_fine_location_indices(i_loc_coarse, downscaling_factors):
    """
    Returns a list of index tuples describing all fine locations associated with
    the the given coarse location according to the given downscaling factors.
    The code is inspired by <https://stackoverflow.com/questions/1208118/
    using-numpy-to-build-an-array-of-all-combinations-of-two-arrays>.

    Parameters
    ----------
    i_loc_coarse : tuple
        Coarse location index.
    downscaling_factors : array of ints
        Downscaling factors for all grid dimensions.

    Returns
    -------
    c : list of tuples
        Fine location index tuples.

    """
    assert len(i_loc_coarse) == len(downscaling_factors)
    a = [np.arange(df * i_loc_coarse[i], df * (i_loc_coarse[i] + 1)) 
        for i, df in enumerate(downscaling_factors)]
    b = np.array(np.meshgrid(*a)).T.reshape(-1, len(downscaling_factors))
    c = [tuple(a) for a in b]
    return c



def only_missing_values_in_at_least_one_dataset(data):
    """
    Tests whether there are only missing values in at least one of the datasets
    included in data.

    Parameters
    ----------
    data : dict of lists of arrays with keys 'obs_hist', 'sim_hist', 'sim_fut'
        Every list represents one climate dataset. Every array represents one
        climate variable.

    Returns
    -------
    result : bool
        Test result.

    """
    for key, array_list in data.items():
        only_missing_values = True
        for a in array_list:
            if not isinstance(a, np.ma.MaskedArray):
                # a is not masked
                only_missing_values = False
                break
            m = a.mask
            if not isinstance(m, np.ndarray):
                # m is a scalar
                if m:
                    # all values in a are masked
                    continue
                else:
                    # no value in a is masked
                    only_missing_values = False
                    break
            if not np.all(m):
                # at least one value in a is not masked
                only_missing_values = False
                break
        if only_missing_values:
            return True
    return False



def only_missing_values_in_at_least_one_time_series(data):
    """
    Tests whether there are only missing values in at least one time series
    included in data.

    Parameters
    ----------
    data : dict of arrays
        Dictionary keys are 'obs_fine', 'sim_coarse', 'sim_coarse_remapbil'.
        Every array represents one climate dataset. The 0th axis of every array
        is considered the time axis.

    Returns
    -------
    result : bool
        Test result.

    """
    for key, a in data.items():
        assert a.ndim in [1, 2], f'{key} array has {a.ndim} dimensions'
        if isinstance(a, np.ma.MaskedArray):
            m = a.mask
            if isinstance(m, np.ndarray):
                if a.ndim == 1:
                    if np.all(m): return True
                else:
                    if np.any(np.all(m, axis=0)): return True
            else:
                if m: return True
    return False



def average_respecting_bounds(x,
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None):
    """
    Average values in x after values <= lower_threshold have been set to
    lower_bound and values >= upper_threshold have been set to upper_bound.

    Parameters
    ----------
    x : array
        Time series to be (de-)randomized.
    lower_bound : float, optional
        Lower bound of values in time series.
    lower_threshold : float, optional
        Lower threshold of values in time series.
    upper_bound : float, optional
        Upper bound of values in time series.
    upper_threshold : float, optional
        Upper threshold of values in time series.

    Returns
    -------
    a : float
        Average.

    """
    y = x.copy()
    if lower_bound is not None and lower_threshold is not None:
        y[y <= lower_threshold] = lower_bound
    if upper_bound is not None and upper_threshold is not None:
        y[y >= upper_threshold] = upper_bound
    return np.mean(y)



def average_valid_values(a, if_all_invalid_use=np.nan,
        lower_bound=None, lower_threshold=None,
        upper_bound=None, upper_threshold=None):
    """
    Returns the average over all valid values in a, where missing/inf/nan values
    are considered invalid, unless there are only invalid values in a, in which
    case if_all_invalid_use is returned. Prior to averaging, values beyond
    threshold are set to the respective bound.

    Parameters
    ----------
    a : array or masked array
        If this is an array then infs and nans in a are replaced.
        If this is a masked array then infs, nans, and missing values in a.data
        are replaced using a.mask to indicate missing values.
    if_all_invalid_use : float, optional
        Used as replacement of invalid values if no valid values can be found.
    lower_bound : float, optional
        Lower bound of values in time series.
    lower_threshold : float, optional
        Lower threshold of values in time series.
    upper_bound : float, optional
        Upper bound of values in time series.
    upper_threshold : float, optional
        Upper threshold of values in time series.

    Returns
    -------
    average : float or array of floats
        Result of averaging. The result is scalar if a is one-dimensional.
        Otherwise the result is an array containing averages for every location.

    """
    # assert that a is a masked array
    if isinstance(a, np.ma.MaskedArray):
        d = a.data
        m = a.mask
        if not isinstance(m, np.ndarray):
            m = np.empty(a.shape, dtype=bool)
            m[:] = a.mask
    else:
        d = a
        m = np.zeros(a.shape, dtype=bool)

    # look for missing values, infs and nans
    l_invalid = np.logical_or(np.logical_or(m, np.isinf(d)), np.isnan(d))

    # compute mean value of all valid values per location
    average1d = lambda x,l : \
        if_all_invalid_use if np.all(l) else average_respecting_bounds(
        x[np.logical_not(l)], lower_bound, lower_threshold,
        upper_bound, upper_threshold)
    space_shape = a.shape[1:]
    if len(space_shape):
        average = np.empty(space_shape, dtype=float)
        for i in np.ndindex(space_shape): 
            j = (slice(None, None),) + i
            average[i] = average1d(d[j], l_invalid[j])
    else:
        average = average1d(d, l_invalid)

    return average



def sample_invalid_values(a, seed=None, if_all_invalid_use=np.nan, warn=False):
    """
    Replaces missing/inf/nan values in a by if_all_invalid_use or by sampling
    from all other values from the same location.

    Parameters
    ----------
    a : array or masked array
        If this is an array then infs and nans in a are replaced.
        If this is a masked array then infs, nans, and missing values in a.data
        are replaced using a.mask to indicate missing values.
    seed : int, optional
        Used to seed the random number generator before replacing invalid
        values.
    if_all_invalid_use : float or array of floats, optional
        Used as replacement of invalid values if no valid values can be found.
    warn : boolean, optional
        Warn user about replacements being made.

    Returns
    -------
    d_replaced : array
        Result of invalid data replacement.
    l_invalid : array
        Boolean array indicating indices of replacement.

    """
    ## # make sure types and shapes of a and if_all_invalid_use fit
    ## space_shape = a.shape[1:]
    ## if len(space_shape):
    ##     msg = 'expected if_all_invalid_use to be an array'
    ##     assert isinstance(if_all_invalid_use, np.ndarray), msg
    ##     msg = 'shapes of a and if_all_invalid_use do not fit'
    ##     assert if_all_invalid_use.shape == space_shape, msg
    ## else:
    ##     msg = 'expected if_all_invalid_use to be scalar'
    ##     assert np.isscalar(if_all_invalid_use), msg

    # assert that a is a masked array
    if isinstance(a, np.ma.MaskedArray):
        d = a.data
        m = a.mask
        if not isinstance(m, np.ndarray):
            m = np.empty(a.shape, dtype=bool)
            m[:] = a.mask
    else:
        d = a
        m = np.zeros(a.shape, dtype=bool)
    
    # look for missing values
    l_invalid = m
    n_missing = np.sum(l_invalid)
    if n_missing:
        msg = 'found %i missing value(s)'%n_missing
        if warn: warnings.warn(msg)
    
    # look for infs
    l_inf = np.isinf(d)
    n_inf = np.sum(l_inf)
    if n_inf:
        msg = 'found %i inf(s)'%n_inf
        if warn: warnings.warn(msg)
        l_invalid = np.logical_or(l_inf, l_invalid)
    
    # look for nans
    l_nan = np.isnan(d)
    n_nan = np.sum(l_nan)
    if n_nan:
        msg = 'found %i nan(s)'%n_nan
        if warn: warnings.warn(msg)
        l_invalid = np.logical_or(l_nan, l_invalid)
    
    # return d if all values are valid
    n_invalid = np.sum(l_invalid)
    if not n_invalid:
        return d, None
    
    ## # otherwise replace invalid values location by location
    ## if len(space_shape):
    ##     d_replaced = np.empty_like(d)
    ##     for i in np.ndindex(space_shape): 
    ##         j = (slice(None, None),) + i
    ##         d_replaced[j] = sample_invalid_values_core(
    ##             d[j], seed, if_all_invalid_use[i], warn, l_invalid[j])
    ## else:
    d_replaced = sample_invalid_values_core(
        d, seed, if_all_invalid_use, warn, l_invalid)

    return d_replaced, l_invalid



def sample_invalid_values_core(d, seed, if_all_invalid_use, warn, l_invalid):
    """
    Replaces missing/inf/nan values in d by if_all_invalid_use or by sampling
    from all other values.

    Parameters
    ----------
    d : array
        Containing values to be replaced.
    seed : int
        Used to seed the random number generator before sampling.
    if_all_invalid_use : float
        Used as replacement of invalid values if no valid values can be found.
    warn : boolean
        Warn user about replacements being made.
    l_invalid : array
        Indicating which values in a are invalid and hence to be replaced.

    Returns
    -------
    d_replaced : array
        Result of invalid data replacement.

    """
    # return d if all values in d are valid
    n_invalid = np.sum(l_invalid)
    if not n_invalid:
        return d

    # no sampling possible if there are no valid values in d
    n_valid = d.size - n_invalid
    if not n_valid:
        msg = 'found no valid value(s)'
        if np.isnan(if_all_invalid_use):
            raise ValueError(msg)
        else:
            msg += ': setting them all to %f'%if_all_invalid_use
            if warn: warnings.warn(msg)
            d_replaced = np.empty_like(d)
            d_replaced[:] = if_all_invalid_use
            return d_replaced

    # replace invalid values by sampling from valid values
    # shuffle sampled values to mimic trend in valid values
    msg = 'replacing %i invalid value(s)'%n_invalid + \
    ' by sampling from %i valid value(s)'%n_valid
    if warn: warnings.warn(msg)
    l_valid = np.logical_not(l_invalid)
    d_valid = d[l_valid]
    if seed is not None: np.random.seed(seed)
    p_sampled = np.random.random_sample(n_invalid)
    d_sampled = percentile1d(d_valid, p_sampled)
    d_replaced = d.copy()
    if n_valid == 1:
        d_replaced[l_invalid] = d_sampled
    else:
        i_valid = np.where(l_valid)[0]
        r_valid = np.argsort(np.argsort(d_valid))
        r_valid_interp1d = spi.interp1d(i_valid, r_valid, fill_value='extrapolate')
        i_sampled = np.where(l_invalid)[0]
        r_sampled = np.argsort(np.argsort(r_valid_interp1d(i_sampled)))
        d_replaced[l_invalid] = np.sort(d_sampled)[r_sampled]
    return d_replaced



def convert_datetimes(datetimes, to):
    """
    Converts a sequence of datetime objects.

    Parameters
    ----------
    datetimes : sequence of datetime objects
        Conversion source.
    to : str
        Conversion target.

    Returns
    -------
    converted_datetimes : array of ints
        Conversion result.

    """
    if to == 'month_number':
        return np.array([d.month for d in datetimes]).astype(np.uint8)
    elif to == 'year':
        return np.array([d.year for d in datetimes]).astype(np.int16)
    elif to == 'day_of_year':
        return np.array([d.timetuple()[7] for d in datetimes]).astype(np.uint16)
    else:
        raise ValueError(f'cannot convert to {to}')



def time_range_to_iris_constraints(time_range):
    """
    Transforms a time range to an iris time constraint.

    Parameters
    ----------
    time_range : str
        Time range specified in format start+'-'+end, where start and end are
        each of format %Y%m%dT%H%M%S or parts thereof.

    Returns
    -------
    c : iris.Constraint
        Time constraint for constraining iris cubes.

    """
    if time_range is None:
      return None
    time_tuples = []
    for t in time_range.split('-'):
      try:
        tt = dt.datetime.strptime(t, '%Y').timetuple()[:1]
      except:
        try:
          tt = dt.datetime.strptime(t, '%Y%m').timetuple()[:2]
        except:
          try:
            tt = dt.datetime.strptime(t, '%Y%m%d').timetuple()[:3]
          except:
            try:
              tt = dt.datetime.strptime(t, '%Y%m%dT%H').timetuple()[:4]
            except:
              try:
                tt = dt.datetime.strptime(t, '%Y%m%dT%H%M').timetuple()[:5]
              except:
                try:
                  tt = dt.datetime.strptime(t, '%Y%m%dT%H%M%S').timetuple()[:6]
                except:
                  tt = None
      if tt is None:
        raise ValueError('unable to turn '+t+' into datetime object')
      time_tuples.append(tt)
    pdts = [iris.time.PartialDateTime(*tt) for tt in time_tuples]
    return iris.Constraint(time=lambda cell: pdts[0] <= cell <= pdts[1])



def load_cube(path, var, time_range=None):
    """
    Loads iris cube.

    Parameters
    ----------
    path : str
        Path to netcdf file.
    var : str
        Standard name of variable to be loaded.
    time_range : str, optional
        Time range to be used for constrained loading.

    Returns
    -------
    cube : iris cube
        Requested cube.

    """
    constraint = time_range_to_iris_constraints(time_range)
    cube = iris.load_cube(path, var if constraint is None else var & constraint)
    return cube



def get_anonymous_dimension_indices(c):
    """
    Returns a list of indices of anonymous dimensions of a given iris cube.

    Parameters
    ----------
    c : iris cube
        Cube whose anonymous dimensions are to be found.

    Returns
    -------
    l : list
        List of anonymous dimension indices.

    """
    all_dims = set(range(c.ndim))
    covered_dims = set(c.coord_dims(coord)[0] for coord in c.dim_coords)
    anonymous_dims = all_dims - covered_dims
    return list(anonymous_dims)



def name_first_anonymous_dimension(c, standard_name):
    """
    Adds a new dimension coordinate to an iris cube. The new coordinate has
    integer values from 0 to the length of the first anonymous dimension of 
    the iris cube minus 1.

    Parameters
    ----------
    c : iris cube
        Cube whose first anonymous dimension is to be named. Is changed
        in-place.
    standard_name : str
        Standard name of the new dimension coordinate.

    """
    if standard_name is not None:
        indices = get_anonymous_dimension_indices(c)
        if len(indices):
            i = indices[0]
            c.add_dim_coord(iris.coords.DimCoord(np.arange(c.shape[i]),
                standard_name=standard_name), i)



def adjust_copula_mbcn(x, rotation_matrices=[], n_quantiles=50):
    """
    Applies the MBCn algorithm for an adjustment of the multivariate rank
    distribution of x['sim_fut'].

    Parameters
    ----------
    x : dict of lists of arrays with keys 'obs_hist', 'sim_hist', 'sim_fut'
        Every list holds n arrays and every array represents the time series for
        one climate variable.
    rotation_matrices : list of (n,n) ndarrays, optional
        List of orthogonal matrices defining a sequence of rotations in variable
        space.
    n_quantiles : int, optional
        Number of quantile-quantile pairs used for non-parametric quantile
        mapping.

    Returns
    -------
    x_sim_fut_ba : list of arrays
        Result of copula adjustment.

    """
    # transform values to standard normal distributions per variable
    # stack resulting arrays row wise
    y = {}
    for key in x:
        y[key] = np.stack([sps.norm.ppf((np.argsort(np.argsort(xi))+.5)/xi.size)
                           for xi in x[key]])

    # initialize total rotation matrix
    n_variables = len(x['sim_fut'])
    o_total = np.diag(np.ones(n_variables))

    # iterate
    for o in rotation_matrices:
        o_total = np.dot(o, o_total)

        # rotate data
        for key in y:
            y[key] = np.dot(o, y[key])

        # do univariate non-parametric quantile delta mapping for every variable
        for i in range(n_variables):
            y_sim_hist_old = y['sim_hist'][i].copy()
            y['sim_hist'][i] = map_quantiles_non_parametric_trend_preserving(
                y['obs_hist'][i], y_sim_hist_old, y_sim_hist_old,
                'additive', n_quantiles)
            y['sim_fut'][i] = map_quantiles_non_parametric_trend_preserving(
                y['obs_hist'][i], y_sim_hist_old, y['sim_fut'][i],
                'additive', n_quantiles)

    # rotate back to original axes
    y['sim_fut'] = np.dot(o_total.T, y['sim_fut'])

    # shuffle x_sim_fut according to the result of the copula adjustment
    x_sim_fut_ba = []
    for i in range(n_variables):
        r_sim_fut_ba = np.argsort(np.argsort(y['sim_fut'][i]))
        x_sim_fut_ba.append(np.sort(x['sim_fut'][i])[r_sim_fut_ba])

    return x_sim_fut_ba



def remapbil(coarse, fine):
    """
    Bilinearly interpolates coarse iris cube to grid of fine iris cube.
    Where bilinear interpolation is not possible due to missing values in
    neighbouring coarse grid cells, values obtained with nearest neighbour
    interpolation are filled in. In the end, every fine grid cell fully
    contained in a coarse grid cell with available data should also have
    available data.

    Parameters
    ----------
    coarse : iris cube
        Cube whose data shall be interpolated.
    fine : iris cube
        Cube which defines the target grid of the interpolation.

    Returns
    -------
    bil : iris cube
        Cube with interpolated data.

    """
    # make sure bilinear interplation works well across 180th meridian
    try:
        coarse.coord('longitude').circular = True
    except iris.exceptions.CoordinateNotFoundError:
        warnings.warn('could not find longitude coordinate')

    # interpolate bilinearly
    bil = coarse.regrid(fine, iris.analysis.Linear(extrapolation_mode='mask'))

    # use coarse values where bilinear interpolation is not possible
    if np.ma.is_masked(bil.data):
        nn = coarse.regrid(fine, iris.analysis.Nearest(
            extrapolation_mode='extrapolate'))
        # nn has the same values as coarse but on the fine grid
        mask = np.logical_and(bil.data.mask, np.logical_not(nn.data.mask))
        # mask is True where bil is missing and nn is not missing
        bil.data.data[mask] = nn.data.data[mask]
        bil.data.mask = nn.data.mask

    return bil



def generateCREmatrix(n):
    """
    Returns a random orthogonal n x n matrix from the circular real ensemble
    (CRE), see Mezzadri (2007) <http://arxiv.org/abs/math-ph/0609050v2>

    Parameters
    ----------
    n : int
        Number of rows and columns of the CRE matrix.

    Returns
    -------
    m : (n,n) ndarray
        CRE matrix.

    """
    z = np.random.randn(n, n)
    q, r = spl.qr(z)  # QR decomposition
    d = np.diagonal(r)
    return q * (d / np.abs(d))



def generate_rotation_matrix_fixed_first_axis(v, transpose=False):
    """
    Generates an n x n orthogonal matrix whose first row or column is equal to
     v/|v|, and whose other rows or columns are found by Gram-Schmidt
    orthogonalisation of v and the standard unit vectors except the first.

    Parameters
    ----------
    v : (n,) array
        Array of n non-zero numbers.
    transpose : boolean, optional
        If True/False generate an n x n orthogonal matrix whose first row/column
        is equal to v/|v|.

    Returns
    -------
    m : (n,n) ndarray
        Rotation matrix.

    """
    assert np.all(v > 0), 'all elements of v have to be positive'

    # generate matrix of vectors that span the R^n with v being the first vector
    a = np.diag(np.ones_like(v))
    a[:,0] = v

    # use QR decomposition for Gram-Schmidt orthogonalisation of these vectors
    q, r = spl.qr(a)

    return -q.T if transpose else -q



def get_downscaling_factors(shape_fine, shape_coarse):
    """
    Computes the downscaling factor for every grid dimension.

    Parameters
    ----------
    shape_fine : tuple
        Shape of the fine resolution grid.
    shape_coarse : tuple
        Shape of the coarse resolution grid.

    Returns
    -------
    downscaling_factors : array of ints
        The downscaling factors.

    """
    msg = 'number of spatial dimensions differs between fine and coarse grid'
    assert len(shape_fine) == len(shape_coarse), msg

    msg = 'downscaling factors are not all integers'
    assert np.all(np.array(shape_fine) % np.array(shape_coarse) == 0), msg

    downscaling_factors = np.array(shape_fine) // np.array(shape_coarse)
    return downscaling_factors



def flatten_all_dimensions_but_first(a):
    """
    Flattens all dimensions but the first of a multidimensional array.

    Parameters
    ----------
    a : ndarray
        Array to be flattened.

    Returns
    -------
    b : ndarray
        Result of flattening, two-dimensional.

    """
    s = a.shape
    s_flattened = (s[0], np.prod(s[1:]))
    return a.reshape(*s_flattened)
