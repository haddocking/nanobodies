###Representation of the paratope differences between kinked and extended nanobodies (Supplementary Figure 5), residues involved in the two different clustered paratope representations (Figure 6A) and docking succes rates with the new mixed paratope restraints (Figure 6B)###

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import re

def int_incode_sort(el):
    #Gets the residue number, even with an insertion code
    if len(el) == 1:
        return int(el)
    elif re.findall(r"[A-Z]", el):
        return int(el[:-1])
    else:
        return int(el)


#SUPPLEMENTARY FIGURE 5
#Percent of ‘kinked’ and ‘extended’ CDR3 conformation structures that involve each nanobody region in antigen binding.

probability_regions = pd.read_csv(Path("..", "data", "probability_regions_interaction.tsv"), sep="\t")
probability_regions_kinked = probability_regions[probability_regions["angle_class"] == "kinked"]
probability_regions_extended = probability_regions[probability_regions["angle_class"] == "extended"]

#Now we plot a violin plot with all angle_classs
fig, ax = plt.subplots(1, 1, figsize=(12, 6))

fr1_probability_kinked, fr1_probability_extended = len([ fr1 for fr1 in probability_regions_kinked["fr1"] if fr1 > 0]), len([ fr1 for fr1 in probability_regions_extended["fr1"] if fr1 > 0])
cdr1_probability_kinked, cdr1_probability_extended = len([ cdr1 for cdr1 in probability_regions_kinked["cdr1"] if cdr1 > 0]), len([ cdr1 for cdr1 in probability_regions_extended["cdr1"] if cdr1 > 0])
fr2_probability_kinked, fr2_probability_extended = len([ fr2 for fr2 in probability_regions_kinked["fr2"] if fr2 > 0]), len([ fr2 for fr2 in probability_regions_extended["fr2"] if fr2 > 0])
cdr2_probability_kinked, cdr2_probability_extended = len([ cdr2 for cdr2 in probability_regions_kinked["cdr2"] if cdr2 > 0]), len([ cdr2 for cdr2 in probability_regions_extended["cdr2"] if cdr2 > 0])
fr3_probability_kinked, fr3_probability_extended = len([ fr3 for fr3 in probability_regions_kinked["fr3"] if fr3 > 0]), len([ fr3 for fr3 in probability_regions_extended["fr3"] if fr3 > 0])
cdr3_probability_kinked, cdr3_probability_extended = len([ cdr3 for cdr3 in probability_regions_kinked["cdr3"] if cdr3 > 0]), len([ cdr3 for cdr3 in probability_regions_extended["cdr3"] if cdr3 > 0])
fr4_probability_kinked, fr4_probability_extended = len([ fr4 for fr4 in probability_regions_kinked["fr4"] if fr4 > 0]), len([ fr4 for fr4 in probability_regions_extended["fr4"] if fr4 > 0])

fr1_probability_kinked, fr1_probability_extended = fr1_probability_kinked/probability_regions_kinked.shape[0]*100, fr1_probability_extended/probability_regions_extended.shape[0]*100
cdr1_probability_kinked, cdr1_probability_extended = cdr1_probability_kinked/probability_regions_kinked.shape[0]*100, cdr1_probability_extended/probability_regions_extended.shape[0]*100
fr2_probability_kinked, fr2_probability_extended = fr2_probability_kinked/probability_regions_kinked.shape[0]*100, fr2_probability_extended/probability_regions_extended.shape[0]*100
cdr2_probability_kinked, cdr2_probability_extended = cdr2_probability_kinked/probability_regions_kinked.shape[0]*100, cdr2_probability_extended/probability_regions_extended.shape[0]*100
fr3_probability_kinked, fr3_probability_extended = fr3_probability_kinked/probability_regions_kinked.shape[0]*100, fr3_probability_extended/probability_regions_extended.shape[0]*100
cdr3_probability_kinked, cdr3_probability_extended = cdr3_probability_kinked/probability_regions_kinked.shape[0]*100, cdr3_probability_extended/probability_regions_extended.shape[0]*100
fr4_probability_kinked, fr4_probability_extended = fr4_probability_kinked/probability_regions_kinked.shape[0]*100, fr4_probability_extended/probability_regions_extended.shape[0]*100

barwidth = 0.4
#We plot points for the probability of the regions interacting
ax.bar([1.3, 3.3, 5.3, 7.3, 9.3, 11.3, 13.3],
         [fr1_probability_kinked, cdr1_probability_kinked, fr2_probability_kinked, cdr2_probability_kinked, fr3_probability_kinked, cdr3_probability_kinked, fr4_probability_kinked],
            color = "blue", label="Kinked", edgecolor = "black", linewidth = 2, width = barwidth)

ax.bar([1.7, 3.7, 5.7, 7.7, 9.7, 11.7, 13.7],
        [fr1_probability_extended, cdr1_probability_extended, fr2_probability_extended, cdr2_probability_extended, fr3_probability_extended, cdr3_probability_extended, fr4_probability_extended],
            color = "red", label="Extended", edgecolor = "black", linewidth = 2, width = barwidth)


ax.set_ylabel("% of structures", size = 15)
ax.set_xticks([1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5])
ax.axvline(2.5, color="black", linestyle="--")
ax.axvline(4.5, color="black", linestyle="--")
ax.axvline(6.5, color="black", linestyle="--")
ax.axvline(8.5, color="black", linestyle="--")
ax.axvline(10.5, color="black", linestyle="--")
ax.axvline(12.5, color="black", linestyle="--")
ax.set_xticklabels(["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"], size = 15)
#we set the y axis from 0 to 1
ax.set_ylim(0, 105)
# yticks with 12 fontsize
ax.tick_params(axis='y', labelsize=12)

#And we add a legend with the angle_classs
leg = ax.legend(loc = "upper left", fontsize=15)

###Now we will perform a chi-squared test to see if the probability of the regions interacting is different between the angle_classes
#We will perform a chi-squared test for each region
fr1_obs = np.array([[len([ fr1 for fr1 in probability_regions_kinked["fr1"] if fr1 > 0]), probability_regions_kinked.shape[0] - len([ fr1 for fr1 in probability_regions_kinked["fr1"] if fr1 > 0])],[len([ fr1 for fr1 in probability_regions_extended["fr1"] if fr1 > 0]), probability_regions_extended.shape[0] - len([ fr1 for fr1 in probability_regions_extended["fr1"] if fr1 > 0])]])
cdr1_obs = np.array([[len([ cdr1 for cdr1 in probability_regions_kinked["cdr1"] if cdr1 > 0]), probability_regions_kinked.shape[0] - len([ cdr1 for cdr1 in probability_regions_kinked["cdr1"] if cdr1 > 0])],[len([ cdr1 for cdr1 in probability_regions_extended["cdr1"] if cdr1 > 0]), probability_regions_extended.shape[0] - len([ cdr1 for cdr1 in probability_regions_extended["cdr1"] if cdr1 > 0])]])
fr2_obs = np.array([[len([ fr2 for fr2 in probability_regions_kinked["fr2"] if fr2 > 0]), probability_regions_kinked.shape[0] - len([ fr2 for fr2 in probability_regions_kinked["fr2"] if fr2 > 0])],[len([ fr2 for fr2 in probability_regions_extended["fr2"] if fr2 > 0]), probability_regions_extended.shape[0] - len([ fr2 for fr2 in probability_regions_extended["fr2"] if fr2 > 0])]])
cdr2_obs = np.array([[len([ cdr2 for cdr2 in probability_regions_kinked["cdr2"] if cdr2 > 0]), probability_regions_kinked.shape[0] - len([ cdr2 for cdr2 in probability_regions_kinked["cdr2"] if cdr2 > 0])],[len([ cdr2 for cdr2 in probability_regions_extended["cdr2"] if cdr2 > 0]), probability_regions_extended.shape[0] - len([ cdr2 for cdr2 in probability_regions_extended["cdr2"] if cdr2 > 0])]])
fr3_obs = np.array([[len([ fr3 for fr3 in probability_regions_kinked["fr3"] if fr3 > 0]), probability_regions_kinked.shape[0] - len([ fr3 for fr3 in probability_regions_kinked["fr3"] if fr3 > 0])],[len([ fr3 for fr3 in probability_regions_extended["fr3"] if fr3 > 0]), probability_regions_extended.shape[0] - len([ fr3 for fr3 in probability_regions_extended["fr3"] if fr3 > 0])]])
cdr3_obs = np.array([[len([ cdr3 for cdr3 in probability_regions_kinked["cdr3"] if cdr3 > 0]), probability_regions_kinked.shape[0] - len([ cdr3 for cdr3 in probability_regions_kinked["cdr3"] if cdr3 > 0])],[len([ cdr3 for cdr3 in probability_regions_extended["cdr3"] if cdr3 > 0]), probability_regions_extended.shape[0] - len([ cdr3 for cdr3 in probability_regions_extended["cdr3"] if cdr3 > 0])]])
fr4_obs = np.array([[len([ fr4 for fr4 in probability_regions_kinked["fr4"] if fr4 > 0]), probability_regions_kinked.shape[0] - len([ fr4 for fr4 in probability_regions_kinked["fr4"] if fr4 > 0])],[len([ fr4 for fr4 in probability_regions_extended["fr4"] if fr4 > 0]), probability_regions_extended.shape[0] - len([ fr4 for fr4 in probability_regions_extended["fr4"] if fr4 > 0])]])

fr1_chi2, fr1_p, fr1_dof, fr1_exp = chi2_contingency(fr1_obs)
cdr1_chi2, cdr1_p, cdr1_dof, cdr1_exp = chi2_contingency(cdr1_obs)
fr2_chi2, fr2_p, fr2_dof, fr2_exp = chi2_contingency(fr2_obs)
cdr2_chi2, cdr2_p, cdr2_dof, cdr2_exp = chi2_contingency(cdr2_obs)
fr3_chi2, fr3_p, fr3_dof, fr3_exp = chi2_contingency(fr3_obs)
cdr3_chi2, cdr3_p, cdr3_dof, cdr3_exp = chi2_contingency(cdr3_obs)
fr4_chi2, fr4_p, fr4_dof, fr4_exp = chi2_contingency(fr4_obs)


###Now we will plot * if p<0.05, ** if p<0.01 and *** if p<0.001
###We will plot it above the highest point of each region and in the middle of both classes
if 0.01<fr1_p < 0.05:
    ax.text(1.5, max(fr1_probability_kinked, fr1_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<fr1_p < 0.01:
    ax.text(1.5, max(fr1_probability_kinked, fr1_probability_extended) + 0.05, "**", ha="center", size=15)
if fr1_p < 0.001:
    ax.text(1.5, max(fr1_probability_kinked, fr1_probability_extended) + 0.05, "***", ha="center", size=15)

if 0.01<cdr1_p < 0.05:
    ax.text(3.5, max(cdr1_probability_kinked, cdr1_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<cdr1_p < 0.01:
    ax.text(3.5, max(cdr1_probability_kinked, cdr1_probability_extended) + 0.05, "**", ha="center", size=15)
if cdr1_p < 0.001:
    ax.text(3.5, max(cdr1_probability_kinked, cdr1_probability_extended) + 0.05, "***", ha="center", size=15)

if 0.01<fr2_p < 0.05:
    ax.text(5.5, max(fr2_probability_kinked, fr2_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<fr2_p < 0.01:
    ax.text(5.5, max(fr2_probability_kinked, fr2_probability_extended) + 0.05, "**", ha="center", size=15)
if fr2_p < 0.001:
    ax.text(5.5, max(fr2_probability_kinked, fr2_probability_extended) + 0.05, "***", ha="center", size=15)

if 0.01<cdr2_p < 0.05:
    ax.text(7.5, max(cdr2_probability_kinked, cdr2_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<cdr2_p < 0.01:
    ax.text(7.5, max(cdr2_probability_kinked, cdr2_probability_extended) + 0.05, "**", ha="center", size=15)
if cdr2_p < 0.001:
    ax.text(7.5, max(cdr2_probability_kinked, cdr2_probability_extended) + 0.05, "***", ha="center", size=15)

if 0.01<fr3_p < 0.05:
    ax.text(9.5, max(fr3_probability_kinked, fr3_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<fr3_p < 0.01:
    ax.text(9.5, max(fr3_probability_kinked, fr3_probability_extended) + 0.05, "**", ha="center", size=15)
if fr3_p < 0.001:
    ax.text(9.5, max(fr3_probability_kinked, fr3_probability_extended) + 0.05, "***", ha="center", size=15)

if 0.01<cdr3_p < 0.05:
    ax.text(11.5, max(cdr3_probability_kinked, cdr3_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<cdr3_p < 0.01:
    ax.text(11.5, max(cdr3_probability_kinked, cdr3_probability_extended) + 0.05, "**", ha="center", size=15)
if cdr3_p < 0.001:
    ax.text(11.5, max(cdr3_probability_kinked, cdr3_probability_extended) + 0.05, "***", ha="center", size=15)

if 0.01<fr4_p < 0.05:
    ax.text(13.5, max(fr4_probability_kinked, fr4_probability_extended) + 0.05, "*", ha="center", size=15)
if 0.001<fr4_p < 0.01:
    ax.text(13.5, max(fr4_probability_kinked, fr4_probability_extended) + 0.05, "**", ha="center", size=15)
if fr4_p < 0.001:
    ax.text(13.5, max(fr4_probability_kinked, fr4_probability_extended) + 0.05, "***", ha="center", size=15)

plt.tight_layout()
plt.savefig(Path(".", "figures", "SI_figure5.png"), dpi=400)



##FIGURE 6A
##Bar plot showing the frequency of each residue in two cluster’s paratopes obtained after hierarchical clustering based on Jaccard distances.

treshold = 0.15
colors = ["green","purple"]
cluster_names = ["CDR", "CDR + FR"]
interacting_regions_grouped = pd.read_csv(Path("..", "data", "clustered_paratope_interacting_residues.tsv"), sep="\t")
#we order by cluster
interacting_regions_grouped = interacting_regions_grouped.sort_values(by="cluster")

fig, ax = plt.subplots (2, 1, figsize=(20, 10))
#the plot above will have the same x axis and the y axis goes from 0 to 1
#we will plot the frequency of each residue
#we count in how many pdbs each residue appears and divide by the total number of pdbs

for c in range(2):
    total_residues = []
    cluster_df = interacting_regions_grouped[interacting_regions_grouped["cluster"] == c+1]
    for n in range(cluster_df.shape[0]):
        regions = cluster_df.iloc[n]["unique_interaction_residues"]
        regions = regions.replace("[", "").replace("]", "").replace("'", "").replace(" ", "").split(",")
        regions = [int_incode_sort(el) for el in regions]
        regions = list(set(regions))
        total_residues += regions

    #we count the frequency of each residue
    residue_freq = {res: total_residues.count(res)/cluster_df.shape[0] for res in set(total_residues)}

    #we color the background in grey between the vertical lines
    ax[c].axvspan(25.6, 32.5, color='grey', alpha=0.3)
    ax[c].axvspan(51.5, 56.5, color='grey', alpha=0.3)
    ax[c].axvspan(94.5, 102.5, color='grey', alpha=0.3)
    ax[c].bar(residue_freq.keys(), residue_freq.values(), color=colors[c])
    ax[c].set_xlim(0, 110)
    ax[c].set_ylim(0, 1)
    ax[c].set_ylabel("Frequency", size = 18)
    ax[c].tick_params(axis='y', labelsize=15)
    ax[c].axhline(y=treshold, color='k', linestyle='--', alpha=0.5)

    #we plot on the top left corner the percentage of pdbs in each cluster
    ax[c].text(0.05, 0.90, f"{cluster_names[c]}: {cluster_df.shape[0]/interacting_regions_grouped.shape[0]*100:.1f}%", transform=ax[c].transAxes, size=20, ha='left', color=colors[c])

    #if a residue frequency is above a certain treshold, we plot it
    for res, freq in residue_freq.items():
        if freq > treshold:
            if res == 100 and c == 0:
                ax[c].text(res+0.1, freq, res, ha='center', va='bottom', color='black', size=9)
            else:
                ax[c].text(res, freq, res, ha='center', va='bottom', color='black', size=9)
            #     ax[c].text(res, freq + 0.01, res, ha='center', va='bottom', color='black', size=8)
            # else:
            #     ax[c].text(res, freq + 0.01, res, ha='center', va='bottom', color='black')
        
        if c == 1:
            # #we plot a star on the top of the bar if the residue is in the list (present in the restraints file)
            if res in [1,2,3,26,27,28,29,30,31,32,33,35,37,39,44,45,47,50,52,53,54,55,56,57,58,59,60,61,103,105,108]:
                ax[c].plot(res, freq-0.02, marker="*", color="white", markersize=10, markeredgecolor="black", markeredgewidth=0.2)

ax[1].set_xlabel("Nanobody regions", size=25)
ax[1].set_xticks([13, 29, 42, 54, 75.5, 98.5, 106])
ax[1].set_xticklabels(["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"], size=20)
ax[0].axes.get_xaxis().set_visible(False)

# add bars to clarify the RX regions
RX_regions = {
    "R1": [[1,3]],
    "R2": [[26,29]],
    "R3": [[30,33]],
    "R4": [[35,39]],
    "R5": [[44,47]],
    "R6": [[50,52]],
    "R7": [[53,55]],
    "R8": [[56,58]],
    "R9": [[59,61]],
    "R10": [[103,108]],
}

# now add to the list the max of the frequencies across that interval
print(residue_freq)
for rx in RX_regions:
    RX_regions[rx].append(max([residue_freq[res] for res in range(RX_regions[rx][0][0], RX_regions[rx][0][1]+1)]))

# plt.savefig(Path(".", "figures", "paratope_representations_presentations.png"), dpi=400)

for rx in RX_regions:
    print("hline for ", rx, RX_regions[rx])
    ax[1].hlines(RX_regions[rx][1] + 0.1, RX_regions[rx][0][0] - 0.5, RX_regions[rx][0][1] + 0.5, color='purple', linewidth=4)
    ax[1].text(RX_regions[rx][0][0] + (RX_regions[rx][0][1] - RX_regions[rx][0][0])/2, RX_regions[rx][1] + 0.1, rx, ha='center', va='bottom', color='purple', size=15)
plt.tight_layout()
plt.savefig(Path(".", "figures", "figure5.png"), dpi=400)



##FIGURE not shown
##Bar plots showing the docking success rates for the IBMu ensemble applying the new paratope restraints.

mixed_df = pd.read_csv(Path("..", "data", "haddock_docking_mixed_paratope_restraints_sr.tsv"), sep = "\t")

#we divide per pipeline stage and antigen
mix_loose_bound_capri_02, mix_loose_bound_capri_06, mix_loose_bound_capri_08 = mixed_df[(mixed_df["scenario"]=="mix_loose") & (mixed_df["antigen"]=="bound") & (mixed_df["capri"]==2)], mixed_df[(mixed_df["scenario"]=="mix_loose") & (mixed_df["antigen"]=="bound") & (mixed_df["capri"]==6)], mixed_df[(mixed_df["scenario"]=="mix_loose") & (mixed_df["antigen"]=="bound") & (mixed_df["capri"]==8)]
mix_loose_unbound_capri_02, mix_loose_unbound_capri_06, mix_loose_unbound_capri_08 = mixed_df[(mixed_df["scenario"]=="mix_loose") & (mixed_df["antigen"]=="unbound") & (mixed_df["capri"]==2)], mixed_df[(mixed_df["scenario"]=="mix_loose") & (mixed_df["antigen"]=="unbound") & (mixed_df["capri"]==6)], mixed_df[(mixed_df["scenario"]=="mix_loose") & (mixed_df["antigen"]=="unbound") & (mixed_df["capri"]==8)]
mix_twohit_bound_capri_02, mix_twohit_bound_capri_06, mix_twohit_bound_capri_08 = mixed_df[(mixed_df["scenario"]=="mix_twohit") & (mixed_df["antigen"]=="bound") & (mixed_df["capri"]==2)], mixed_df[(mixed_df["scenario"]=="mix_twohit") & (mixed_df["antigen"]=="bound") & (mixed_df["capri"]==6)], mixed_df[(mixed_df["scenario"]=="mix_twohit") & (mixed_df["antigen"]=="bound") & (mixed_df["capri"]==8)]
mix_twohit_unbound_capri_02, mix_twohit_unbound_capri_06, mix_twohit_unbound_capri_08 = mixed_df[(mixed_df["scenario"]=="mix_twohit") & (mixed_df["antigen"]=="unbound") & (mixed_df["capri"]==2)], mixed_df[(mixed_df["scenario"]=="mix_twohit") & (mixed_df["antigen"]=="unbound") & (mixed_df["capri"]==6)], mixed_df[(mixed_df["scenario"]=="mix_twohit") & (mixed_df["antigen"]=="unbound") & (mixed_df["capri"]==8)]

#we plot the results
fig,axs = plt.subplots(3, figsize = (8, 15))
barwidth = 0.4

ticks = [0.8, 1.2, 1.8, 2.2, 2.8, 3.2,  4.8,5.2, 5.8, 6.2, 6.8, 7.2]
tick_labels = ["T1 B", "T1 U", "T10 B", "T10 U", "T200 B", "T200 U", "T1 B", "T1 U", "T10 B", "T10 U", "T200 B", "T200 U"]

#Rigid Body stage
axs[0].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_acceptable'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_acceptable'].values[0], mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_acceptable'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_acceptable'].values[0], mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_acceptable'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_acceptable'].values[0]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
axs[0].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_high'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_high'].values[0], mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_high'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_high'].values[0], mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_high'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_high'].values[0]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
axs[0].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_medium'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_medium'].values[0], mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_medium'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_medium'].values[0], mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_medium'].values[0], mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_medium'].values[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
axs[0].axvline(4, color = "black")
axs[0].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_acceptable'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_acceptable'].values[0], mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_acceptable'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_acceptable'].values[0], mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_acceptable'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_acceptable'].values[0]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
axs[0].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_high'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_high'].values[0], mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_high'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_high'].values[0], mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_high'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_high'].values[0]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
axs[0].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_medium'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_medium'].values[0], mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_medium'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_medium'].values[0], mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_medium'].values[0], mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_medium'].values[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
axs[0].set_ylim(0, 100)
axs[0].set_xticks([])

#we plot just above the bars the SR values
axs[0].text(0.8, mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_acceptable'].values[0]+1, f"{mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(1.2, mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(1.8, mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_acceptable'].values[0]+1, f"{mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(2.2, mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(2.8, mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_acceptable'].values[0]-3, f"{mix_loose_bound_capri_02[mix_loose_bound_capri_02['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8, color = "white")
axs[0].text(3.2, mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_02[mix_loose_unbound_capri_02['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(4.8, mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(5.2, mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(5.8, mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(6.2, mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(6.8, mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_02[mix_twohit_bound_capri_02['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].text(7.2, mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_02[mix_twohit_unbound_capri_02['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[0].set_ylabel("SR (%)", size = 12)

#Flex Ref stage
axs[1].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_acceptable'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_acceptable'].values[0], mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_acceptable'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_acceptable'].values[0], mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_acceptable'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_acceptable'].values[0]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
axs[1].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_high'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_high'].values[0], mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_high'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_high'].values[0], mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_high'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_high'].values[0]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
axs[1].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_medium'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_medium'].values[0], mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_medium'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_medium'].values[0], mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_medium'].values[0], mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_medium'].values[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
axs[1].axvline(4, color = "black")
axs[1].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_acceptable'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_acceptable'].values[0], mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_acceptable'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_acceptable'].values[0], mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_acceptable'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_acceptable'].values[0]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
axs[1].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_high'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_high'].values[0], mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_high'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_high'].values[0], mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_high'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_high'].values[0]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
axs[1].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_medium'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_medium'].values[0], mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_medium'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_medium'].values[0], mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_medium'].values[0], mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_medium'].values[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
axs[1].set_ylim(0, 100)
axs[1].set_xticks([])

#we plot just above the bars the SR values
axs[1].text(0.8, mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_acceptable'].values[0]+1, f"{mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(1.2, mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(1.8, mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_acceptable'].values[0]+1, f"{mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(2.2, mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(2.8, mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_acceptable'].values[0]-3, f"{mix_loose_bound_capri_06[mix_loose_bound_capri_06['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8, color = "white")
axs[1].text(3.2, mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_06[mix_loose_unbound_capri_06['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(4.8, mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(5.2, mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(5.8, mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(6.2, mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(6.8, mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_06[mix_twohit_bound_capri_06['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].text(7.2, mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_06[mix_twohit_unbound_capri_06['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[1].set_ylabel("SR (%)", size = 12)

#EM Ref stage
axs[2].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_acceptable'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_acceptable'].values[0], mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_acceptable'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_acceptable'].values[0], mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_acceptable'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_acceptable'].values[0]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
axs[2].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_high'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_high'].values[0], mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_high'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_high'].values[0], mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_high'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_high'].values[0]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
axs[2].bar([0.8, 1.2, 1.8, 2.2, 2.8, 3.2], [mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_medium'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_medium'].values[0], mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_medium'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_medium'].values[0], mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_medium'].values[0], mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_medium'].values[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
axs[2].axvline(4, color = "black")
axs[2].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_acceptable'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_acceptable'].values[0], mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_acceptable'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_acceptable'].values[0], mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_acceptable'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_acceptable'].values[0]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
axs[2].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_high'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_high'].values[0], mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_high'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_high'].values[0], mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_high'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_high'].values[0]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
axs[2].bar([4.8,5.2, 5.8, 6.2, 6.8, 7.2], [mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_medium'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_medium'].values[0], mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_medium'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_medium'].values[0], mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_medium'].values[0], mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_acceptable'].values[0] - mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_medium'].values[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
axs[2].set_ylim(0, 100)
axs[2].set_xticks([])

#we plot just above the bars the SR values
axs[2].text(0.8, mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_acceptable'].values[0]+1, f"{mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(1.2, mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(1.8, mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_acceptable'].values[0]+1, f"{mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(2.2, mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(2.8, mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_acceptable'].values[0]-3, f"{mix_loose_bound_capri_08[mix_loose_bound_capri_08['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8, color = "white")
axs[2].text(3.2, mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_acceptable'].values[0]+1, f"{mix_loose_unbound_capri_08[mix_loose_unbound_capri_08['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(4.8, mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(5.2, mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==1]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(5.8, mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(6.2, mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==10]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(6.8, mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_acceptable'].values[0]+1, f"{mix_twohit_bound_capri_08[mix_twohit_bound_capri_08['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].text(7.2, mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_acceptable'].values[0]+1, f"{mix_twohit_unbound_capri_08[mix_twohit_unbound_capri_08['top']==200]['sr_acceptable'].values[0]:.2f}", ha = "center", va = "bottom", size = 8)
axs[2].set_ylabel("SR (%)", size = 12)

#we plot the titles for each scenario
axs[0].text(2, 106, "Mix Loose Interface", ha = "center", va = "center", rotation = 0, fontsize = 15)
axs[0].text(6, 106, "Mix Two-Hit Interface", ha = "center", va = "center", rotation = 0, fontsize = 15)

leg = axs[2].legend(loc = "upper right", fontsize = 10)
handles = leg.legend_handles[:3]
labels = ["Acceptable", "Medium", "High"]
leg = axs[2].legend(handles, labels, loc = "lower center", fontsize = 12, ncol = 3, bbox_to_anchor = (0.5, -0.21))

axs[0].text(-1, 50, "RigidBody\nstage", ha = "center", va = "center", rotation = 0, fontsize = 12)
axs[1].text(-1, 50, "FlexRef\nstage", ha = "center", va = "center", rotation = 0, fontsize = 12)
axs[2].text(-1, 50, "EMRef\nstage", ha = "center", va = "center", rotation = 0, fontsize = 12)

plt.tight_layout()

plt.savefig(Path(".", "figures", "mixed_paratope_restraints_sr.png"))