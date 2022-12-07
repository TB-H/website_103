"""
Contains all the functions necessary to perform calculations and generate the figure
"""
import warnings

import numpy as np
import pandas as pd
import scipy.optimize
from bokeh.layouts import column
import scipy.special
import scipy.stats as st
import os

import bebi103

import bokeh.io

rg = np.random.default_rng()
def line_area(a0, k, t):
    """Compute area using linear model"""

    return a0 * (1 + k * t)

def line_resid(params, t, a):
    """Residual for linear model."""
    a0, k = params
    return a - line_area(a0, k, t)

def line_area_mle_lstq(data):
    """Compute MLE for parameters in linear area model."""

    tdata = data[0]
    adata = data[1]

    res = scipy.optimize.least_squares(
        line_resid,
        np.array([35, 1]),
        args=(tdata, adata),
        bounds=([0, 0], [np.inf, np.inf]),
    )

    # Compute residual sum of squares from optimal params
    rss_mle = np.sum(line_resid(res.x, tdata, adata) ** 2)

    # Compute MLE for sigma
    sigma_mle = np.sqrt(rss_mle / len(tdata))

    return tuple([x for x in res.x] + [sigma_mle])

def exp_area(a0, k, t):
    """Compute area using exponential model"""

    return a0 * np.exp(k * t)

def exp_resid(params, t, a):
    """Residual for exp model."""
    a0, k = params
    return a - exp_area(a0, k, t)

def exp_area_mle_lstq(data):
    """Compute MLE for parameters in exponential area model."""

    tdata = data[0]
    adata = data[1]

    res = scipy.optimize.least_squares(
        exp_resid,
        np.array([1, 1]),
        args=(tdata, adata),
        bounds=([0, 0], [np.inf, np.inf]),
    )

    # Compute residual sum of squares from optimal params
    rss_mle = np.sum(exp_resid(res.x, tdata, adata) ** 2)

    # Compute MLE for sigma
    sigma_mle = np.sqrt(rss_mle / len(tdata))

    return tuple([x for x in res.x] + [sigma_mle])

def sample_line(a0, k,sigma, t, size=1):
    """Generate samples of area vs time for linear model."""
    samples = np.empty((size, len(t)))

    for i in range(size):
        mu = line_area(a0, k, t)
        samples[i] = np.maximum(0, rg.normal(mu, sigma))

    return samples

def sample_exp(a0, k,sigma, t, size=1):
    """Generate samples of area vs time for linear model."""
    samples = np.empty((size, len(t)))

    for i in range(size):
        mu = exp_area(a0, k, t)
        samples[i] = np.maximum(0, rg.normal(mu, sigma))

    return samples