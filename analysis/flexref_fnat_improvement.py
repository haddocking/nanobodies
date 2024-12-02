###Representation of models fnat improvement after Flexible Refinement (Figure 4)

import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


##FIGURE 4
##Histogram shows the distribution of Δfnat for all acceptable models after flexible refinement stage.

fnat_freq_data = pd.read_csv(Path("..", "data", "fnat_diff_flexref.tsv"), sep="\t")

fnat_freq_true_bound, fnat_freq_true_unbound = fnat_freq_data[(fnat_freq_data["scenario"] == "True")&(fnat_freq_data["antigen"] == "bound")], fnat_freq_data[(fnat_freq_data["scenario"] == "True")&(fnat_freq_data["antigen"] == "unbound")]
fnat_freq_loose_bound, fnat_freq_loose_unbound = fnat_freq_data[(fnat_freq_data["scenario"] == "Loose")&(fnat_freq_data["antigen"] == "bound")], fnat_freq_data[(fnat_freq_data["scenario"] == "Loose")&(fnat_freq_data["antigen"] == "unbound")]
fnat_freq_twohit_bound, fnat_freq_twohit_unbound = fnat_freq_data[(fnat_freq_data["scenario"] == "Two-hit")&(fnat_freq_data["antigen"] == "bound")], fnat_freq_data[(fnat_freq_data["scenario"] == "Two-hit")&(fnat_freq_data["antigen"] == "unbound")]


#we plot 6 plots, for each scenario and for bound and unbound
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

for j, scenario in enumerate(["True", "Loose", "Two-hit"]):
    for i, anti in enumerate(["bound", "unbound"]):
        data = fnat_freq_data[(fnat_freq_data["scenario"] == scenario)&(fnat_freq_data["antigen"] == anti)]
        fnat_freq = data["freq"].values.tolist()
        fnat_bin = data["bin"].values.tolist()
        
        #we plot the frequency of the data
        axs[i][j].hist(fnat_bin, fnat_bin, weights = fnat_freq)
        axs[i][j].set_xlim(-0.59, 0.59)
        axs[i][j].set_ylim(0, 1)

        #we plot the median in the plot
        axs[i][j].text(data["median"].values[0] + 0.5, 0.33, f"{np.round(data['median'].values[0], 4)}", horizontalalignment='center', verticalalignment='center', transform=axs[i][j].transAxes, fontsize=12, color='black')
        
        #and we plot the minimum and the maximum fnat improvement
        axs[i][j].axvline(data["min"].values[0], color='r', linestyle='dashed', linewidth=1)
        axs[i][j].text(data["min"].values[0]-0.04, 0.01, f"Min: {np.round(data['min'].values[0], 4)}", rotation=90)
        axs[i][j].axvline(data["max"].values[0], color='g', linestyle='dashed', linewidth=1)
        axs[i][j].text(data["max"].values[0]+0.01, 0.01, f"Max: {np.round(data['max'].values[0], 4)}", rotation=90)

        #we print the amount of data points at the top right in a white box
        axs[i][j].text(0.5, 0.95, f"n = {data['n'].values[0]}", horizontalalignment='center', verticalalignment='center', transform=axs[i][j].transAxes)

        axs[0][j].set_xticks([])
        axs[i][1].set_yticks([])
        axs[i][2].set_yticks([])
        axs[1][j].set_xlabel(r'$\Delta$' + "fnat", fontsize=12)

axs[0][0].set_title("True Interface", fontsize=15)
axs[0][1].set_title("Loose Interface", fontsize=15)
axs[0][2].set_title("Two-hit Interface", fontsize=15)
axs[0][0].set_ylabel("Frequency", fontsize=12)
axs[1][0].set_ylabel("Frequency", fontsize=12)

axs[0][0].text(-1.2, 700, "Bound\nAntigen", fontsize=15, horizontalalignment='left')
axs[1][0].text(-1.2, 700, "Unbound\nAntigen", fontsize=15, horizontalalignment='left')

#we change the size so that the last two texts fit the image
plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.9, hspace=0.05, wspace=0.05) #we leave space between the plots

plt.savefig(Path("figures", "fnat_diff_flexref.png"))