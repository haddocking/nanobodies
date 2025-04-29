###Representation of models fnat improvement after Flexible Refinement (Figure 4)

import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


##FIGURE 4
##Histogram shows the distribution of Δfnat for all acceptable models after flexible refinement stage.

fnat_freq_data = pd.read_csv(Path("..", "data", "fnat_diff_data.tsv"), sep="\t")

# fnat_freq_loose_bound, fnat_freq_loose_unbound = fnat_freq_data[(fnat_freq_data["scenario"] == "Loose")&(fnat_freq_data["antigen"] == "bound")], fnat_freq_data[(fnat_freq_data["scenario"] == "Loose")&(fnat_freq_data["antigen"] == "unbound")]
# fnat_freq_twohit_bound, fnat_freq_twohit_unbound = fnat_freq_data[(fnat_freq_data["scenario"] == "Two-hit")&(fnat_freq_data["antigen"] == "bound")], fnat_freq_data[(fnat_freq_data["scenario"] == "Two-hit")&(fnat_freq_data["antigen"] == "unbound")]

#we plot 6 plots, for each scenario and for bound and unbound
fig, axs = plt.subplots(2, 3, figsize=(15, 8))

colors = {
    "real b" : "lightgreen",
    "real u" : "lightgreen",
    "loose b" : "lightblue",
    "loose u" : "lightblue",
    "twohit b" : "darkorange",
    "twohit u" : "darkorange"
}

bins = np.arange(-0.6, 0.6, 0.05)

for j, scenario in enumerate(["real", "loose", "twohit"]):
    for i, anti in enumerate(["b", "u"]):
        print(scenario, anti)
        fnat_freq = fnat_freq_data[(fnat_freq_data["scenario"] == scenario)&(fnat_freq_data["struct"] == anti)]["fnat_diff"].values
        #we plot the frequency of the data
        if anti == "b":
            alpha = 0.9
        else:
            alpha = 0.6
        axs[i][j].hist(fnat_freq, bins, density = True, color=colors[f"{scenario} {anti}"], alpha=alpha)
        axs[i][j].set_xlim(-0.6, 0.6)
        axs[i][j].set_ylim(0, 6.2)

        # we plot the median in the plot
        median_ij = np.median(fnat_freq)
        print(f"Median for {scenario} {anti}: {median_ij}")
        axs[i][j].text(
            median_ij + 0.5, 0.05,
            f"M = {median_ij:.3f}",
            horizontalalignment='center',
            verticalalignment='center',
            transform=axs[i][j].transAxes,
            fontsize=15,
            color='black',
            font="Arial")
        #axs[i][j].text(data["mean"].values[0] + 0.5, 0.33, f"{np.round(data['mean'].values[0], 3)}", horizontalalignment='center', verticalalignment='center', transform=axs[i][j].transAxes, fontsize=12, color='black')
        
        #and we plot the minimum and the maximum fnat improvement
        min_fnat_diff = np.min(fnat_freq)
        max_fnat_diff = np.max(fnat_freq)
        axs[i][j].axvline(min_fnat_diff, color='black', linestyle='dashed', linewidth=1)
        axs[i][j].text(min_fnat_diff-0.05, 0.25, f"Min: {np.round(min_fnat_diff, 4)}", rotation=90, fontsize=12)
        axs[i][j].axvline(max_fnat_diff, color='black', linestyle='dashed', linewidth=1)
        axs[i][j].text(max_fnat_diff+0.015, 0.25, f"Max: {np.round(max_fnat_diff, 4)}", rotation=90, fontsize=12)
# 
        #we print the amount of data points at the top right in a white box
        axs[i][j].text(0.5, 0.95, f"n = {fnat_freq.shape[0]}", horizontalalignment='center', verticalalignment='center', transform=axs[i][j].transAxes, fontsize=15, font="Arial")
# 
        axs[0][j].set_xticks([])
        axs[i][1].set_yticks([])
        axs[i][2].set_yticks([])
        axs[1][j].set_xlabel("$\Delta F_{nat}$", fontsize=18)

axs[0][0].set_title("True Interface", fontsize=18)
axs[0][1].set_title("Loose Interface", fontsize=18)
axs[0][2].set_title("Two-hit Interface", fontsize=18)
axs[0][0].set_ylabel("Density", fontsize=16)
axs[1][0].set_ylabel("Density", fontsize=16)

axs[0][0].text(-1.1, 3.1, "Bound\nAntigen", fontsize=18, horizontalalignment='left')
axs[1][0].text(-1.1, 3.1, "Unbound\nAntigen", fontsize=18, horizontalalignment='left')

#we change the size so that the last two texts fit the image
plt.subplots_adjust(#top=0.9,
                   #bottom=0.1,
                   #left=0.01,
                   #right=0.9,
                   hspace=0.05,
                   wspace=0.05
                   ) #we leave space between the plots
# do not have space between the subplots
plt.subplots_adjust(hspace=0, wspace=0)
plt.tight_layout()
plt.savefig(Path("figures", "figure3.png"), dpi=400)