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
from homework9_2 import *

# Data path
data_path = "../datasets"

# Import dataset
df = pd.read_csv(os.path.join(data_path, "caulobacter_growth_events.csv"))

# Save the times and area for each bacterium and each growth event, resetting time as list of numpy arrays indexed by event
times1 = [
    df.loc[(df["bacterium"] == 1) & (df["growth event"] == i), "time (min)"].values - df.loc[(df["bacterium"] == 1) & (df["growth event"] == i), "time (min)"].values[0]+1 for i in range(19)
]
areas1 = [
    df.loc[(df["bacterium"] == 1) & (df["growth event"] == i), "area (µm²)"].values for i in range(19)
]
times2 = [
    df.loc[(df["bacterium"] == 2) & (df["growth event"] == i), "time (min)"].values - df.loc[(df["bacterium"] == 2) & (df["growth event"] == i), "time (min)"].values[0]+1 for i in range(max(df["growth event"])+1)
]
areas2 = [
    df.loc[(df["bacterium"] == 2) & (df["growth event"] == i), "area (µm²)"].values for i in range(max(df["growth event"])+1)
]

# Define empty list for all eight sets to calculate MLEs for linear + exp
lin_mles1 = []
lin_mles2 = []
exp_mles1 = []
exp_mles2 = []

# Calculate MLEs for bacterium 1
for i in range(19):
    # Extract the times and areas
    time1 = times1[i]
    area1 = areas1[i]
    
    # Append MLES
    lin_mles1.append(line_area_mle_lstq((time1, area1)))
    exp_mles1.append(exp_area_mle_lstq((time1, area1)))

# and 2
for i in range(max(df["growth event"])):

    # Extract times and areas
    time2 = times2[i]
    area2 = areas2[i]

    # Append MLEs
    lin_mles2.append(line_area_mle_lstq((time2, area2)))
    exp_mles2.append(exp_area_mle_lstq((time2, area2)))


# Generate samples
lin_samples1_d = [
    sample_line(*lin_mles1[i], times1[i], size=10000)
    for i in range(len(lin_mles1))
]
lin_samples2_d = [
    sample_line(*lin_mles2[i], times2[i], size=10000)
    for i in range(len(lin_mles2))
]
exp_samples1_d = [
    sample_exp(*exp_mles1[i], times1[i], size=10000)
    for i in range(len(exp_mles1))
]
exp_samples2_d = [
    sample_exp(*exp_mles2[i], times2[i], size=10000)
    for i in range(len(exp_mles2))
]

plots = []

for i in range(len(lin_samples1_d)):
    p = bokeh.plotting.figure(
        height=200,
        width=200,
        x_axis_label="time (min)",
        y_axis_label="diff. area (µm²)",
        title="Bact 1 Lin For GE {}".format(i),
    )
    bebi103.viz.predictive_regression(
        samples=lin_samples1_d[i] + i,
        samples_x=times1[i],
        data=(times1[i], areas1[i] + i),
        diff=True,
        p=p,
    )
    plots.append(p)
    p = bokeh.plotting.figure(
        height=200,
        width=200,
        x_axis_label="time (min)",
        y_axis_label="diff. area (µm²)",
        title="Bact 1 Exp For GE {}".format(i),
    )
    bebi103.viz.predictive_regression(
        samples=exp_samples1_d[i] + i,
        samples_x=times1[i],
        data=(times1[i], areas1[i] + i),
        diff=True,
        p=p,
    )
    plots.append(p)


bokeh.io.show(bokeh.layouts.gridplot(plots, ncols=6))