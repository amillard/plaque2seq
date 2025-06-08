


import pandas as pd 


df = pd.read_csv("Coli_Run1.tsv", sep="\t")




#######
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Data
x = df["Yeild MB"]
y = df ["Qubit [ng/uL]"]


# Linear regression
slope, intercept, r_value, _, _ = linregress(x, y)
line = slope * x + intercept
r_squared = r_value**2

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(x, y, color="darkred", s=55)            # Red points

plt.plot(x, line, color="black")  # Black line

# R² annotation
plt.text(0.05, 0.95, f"$R^2 = {r_squared:.3f}$", transform=plt.gca().transAxes,
        fontsize=14, verticalalignment='top')

# Labels and legend
plt.xlabel("DNA input (ng)", fontsize=12)
plt.ylabel("Sequencing Yield (Mb) ", fontsize=12)
plt.title("Squencing Yield verus DNA input", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.legend(fontsize=12)

plt.tight_layout()
plt.show()





###################


import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import linregress

# Sample DataFrame assumed to already be loaded as `df`

# Setup the figure
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# --- Panel 1: Boxplot ---
sns.boxplot(x=df["Host Contamination"], ax=axs[0], color="darkred")
axs[0].set_title("Bacterial DNA contamination per library", fontsize=14)
axs[0].set_xlabel(r"% of reads mapping to $ \it{E. coli} $ MG1655", fontsize=12)
axs[0].set_ylabel("Sample", fontsize=12)
axs[0].tick_params(axis='both', labelsize=12)

# --- Panel 2: Scatter plot with regression ---
x = df["Yeild MB"]
y = df["Qubit [ng/uL]"]

slope, intercept, r_value, _, _ = linregress(x, y)
line = slope * x + intercept
r_squared = r_value**2

axs[1].scatter(x, y, color="darkred", s=55, label="Data")
axs[1].plot(x, line, color="black", label="Fit")

# R² annotation
axs[1].text(0.05, 0.95, f"$R^2 = {r_squared:.3f}$", transform=axs[1].transAxes,
            fontsize=14, verticalalignment='top')

# Labels and formatting
axs[1].set_xlabel("DNA input (ng)", fontsize=12)
axs[1].set_ylabel("Sequencing Yield (Mb)", fontsize=12)
axs[1].set_title("Sequencing Yield versus DNA input", fontsize=14)
axs[1].tick_params(axis='both', labelsize=12)
#axs[1].legend(fontsize=12)

# Layout
plt.tight_layout()
plt.savefig("Ecoli_data.svg", format="svg")
plt.show()
