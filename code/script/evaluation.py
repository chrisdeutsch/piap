import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("publication")

data = pd.read_csv("observables.csv")
grouped = data[["beta", "obs"]].groupby("beta")


# Mean, stddev and quantiles
avg = grouped.mean()
stddev = grouped.std().rename(columns={"obs": "stddev"})

quantiles = grouped.quantile([0.05, 0.95]).unstack()["obs"]
quantiles = quantiles.rename(columns={0.05: "obs_5", 0.95: "obs_95"})

result = pd.concat([avg, stddev, quantiles], axis=1)

# Plot of the mean
line_avg, = plt.plot(result.obs, label="mean")

# 1 sigma band
std_lower = result.obs.values - result.stddev.values
std_upper = result.obs.values + result.stddev.values

plt.fill_between(result.index, std_lower, std_upper, 
                 facecolor=line_avg.get_color(), alpha=0.25,
                 label="$1\sigma$-range")

# 90% confidence band
plt.fill_between(result.index, std_upper, result.obs_95, facecolor='0.5',
                 alpha=0.25)
plt.fill_between(result.index, result.obs_5, std_lower, facecolor='0.5',
                 alpha=0.25, label="$90\%$-confidence")

# Logscale
plt.xscale('log')

# Labeling
plt.legend(loc="lower left")
plt.xlabel(r"Temperature~$\beta$")
plt.ylabel(r"average pair-distance~$d$")
plt.title(r"Temperature dependence")

plt.savefig("temperature_dep.pdf")
