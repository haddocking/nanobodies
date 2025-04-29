###Representation of docking success rate (SR) for different ensmembles and scenarios (Fig.3, Supplementary Fig. 2 & 3)

import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

##FIGURE 3A 
##Bar plot showing individual models HADDOCK performance on nanobody-antigen docking after Flexible Refinement stage.

##We read the data from the Flexible Refinement stage
sr_data = pd.read_csv(Path("..", "data", "haddock_docking_sr_per_scenario.tsv"), sep = "\t")
#and from the AlphaFold2-Multimer and AlphaFold3 runs
alphafold_multimer = pd.read_csv(Path("..", "data/alphafold2multimer_25_predictions_sr.tsv"), sep = "\t")
# alphafold3
# alphafold3_unseen = pd.read_csv(Path("..", "data/alphafold3_unseen_25_predictions_sr.tsv"), sep = "\t")
alphafold3 = pd.read_csv(Path("..", "data/alphafold3_25_predictions_sr.tsv"), sep = "\t")
print(f"number of alphafold3 nanobodies {alphafold3['pdb'].unique().shape[0]}")

# print average dockq for af2m when rank == 1
print(f"AF2M dockq rank 1 {alphafold_multimer[alphafold_multimer['rank'] == 1]['dockq'].mean():.3f}")
# same for alphafold3
print(f"AF3 dockq rank 1 {alphafold3[alphafold3['rank'] == 1]['dockq'].mean():.3f}")

##We separate the data by scenario and antigen
true_bound = sr_data[(sr_data["scenario"]=="real")&(sr_data["struct"]=="b")]
true_unbound = sr_data[(sr_data["scenario"]=="real")&(sr_data["struct"]=="u")]
loose_bound = sr_data[(sr_data["scenario"]=="loose")&(sr_data["struct"]=="b")]
loose_unbound = sr_data[(sr_data["scenario"]=="loose")&(sr_data["struct"]=="u")]
twohit_bound = sr_data[(sr_data["scenario"]=="twohit")&(sr_data["struct"]=="b")]
twohit_unbound = sr_data[(sr_data["scenario"]=="twohit")&(sr_data["struct"]=="u")]

###We represent in three different plots the values for top1, top10 and top200 for the four ensembles and bound/unbound antigen runs, showing acceptable/medium/high

fig,axs = plt.subplots(3,5, figsize = (15, 12), width_ratios=[3, 3, 3, 3, 2])
# the fifth plot should be smaller

#xticks = [0.3, 0.7, 1.3, 1.7, 2.3, 2.7, 3.3, 3.7, 4.3, 4.7, 5.3, 5.7, 6.3, 6.7, 7.3, 7.7, 8.3, 8.7, 9.3, 9.7, 10.3, 10.7, 11.3, 11.7, 12.4, 12.8, 13.4, 13.8]
xlabels = ["B T1", "B T10", "B T200", "U T1", "U T10", "U T200"]
af3_xlabels = ["AFMu T1", "AFMu T10", "AF3 T1", "AF3 T10"]
barwidth = 0.4
ensembles = ["IB", "IBMo", "IBMu", "IBMoMu"]
scenarios = ["real", "loose", "twohit"]

bars_b = [0.3, 1.3, 2.3]
bars_u = [0.7, 1.7, 2.7]
text_fontsize = 9

# This full plot will go to the SI
sr_data_flex = sr_data[sr_data["stage"]=="flex"]
for i, ens in enumerate(ensembles):
    for j, scen in enumerate(scenarios):
        axs[j,i].set_ylim(0, 105)
        if j != 2:
            axs[j,i].set_xticks([])
        if i != 0:
            axs[j,i].set_yticks([])
        else:
            axs[j,i].set_ylabel("SR (%)", size = 15)
        ranks = [1, 10, 200]
        axs[j,i].bar(bars_b, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="b")]["acc_sr"].values[0]*100 for rank in ranks], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
        axs[j,i].bar(bars_b, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="b")]["med_sr"].values[0]*100 for rank in ranks], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
        axs[j,i].bar(bars_b, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="b")]["high_sr"].values[0]*100 for rank in ranks], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
        axs[j,i].bar(bars_u, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="u")]["acc_sr"].values[0]*100 for rank in ranks], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
        axs[j,i].bar(bars_u, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="u")]["med_sr"].values[0]*100 for rank in ranks], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
        axs[j,i].bar(bars_u, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="u")]["high_sr"].values[0]*100 for rank in ranks], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
        # let's label the acceptable SR
        for b_idx, b in enumerate(bars_b):
            b_value = sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==ranks[b_idx])&(sr_data_flex["struct"]=="b")]["acc_sr"].values[0]*100
            u_value = sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==ranks[b_idx])&(sr_data_flex["struct"]=="u")]["acc_sr"].values[0]*100
            axs[j,i].text(b, b_value , f"{b_value:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize)
            axs[j,i].text(bars_u[b_idx], u_value , f"{u_value:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize)
        

#we create lists with the values of the bars
lpdbs = alphafold_multimer["pdb"].unique().shape[0]
af2m_t1_acc = alphafold_multimer[(alphafold_multimer["rank"]==1)&(alphafold_multimer["capri"]!="incorrect")]
af2m_t1_med = alphafold_multimer[(alphafold_multimer["rank"]==1)&(alphafold_multimer["capri"]!="incorrect")&(alphafold_multimer["capri"]!="acceptable")]
af2m_t1_high = alphafold_multimer[(alphafold_multimer["rank"]==1)&(alphafold_multimer["capri"]=="high")]

af2m_t1_acc_sr = af2m_t1_acc["pdb"].unique().shape[0]/lpdbs*100
af2m_t1_med_sr = af2m_t1_med["pdb"].unique().shape[0]/lpdbs*100
af2m_t1_high_sr = af2m_t1_high["pdb"].unique().shape[0]/lpdbs*100
af_multimer_rank1_list = [af2m_t1_acc_sr, af2m_t1_med_sr, af2m_t1_high_sr]

af2m_t10_acc = alphafold_multimer[(alphafold_multimer["rank"]<11)&(alphafold_multimer["capri"]!="incorrect")]
af2m_t10_med = alphafold_multimer[(alphafold_multimer["rank"]<11)&(alphafold_multimer["capri"]!="incorrect")&(alphafold_multimer["capri"]!="acceptable")]
af2m_t10_high = alphafold_multimer[(alphafold_multimer["rank"]<11)&(alphafold_multimer["capri"]=="high")]

af2m_t10_acc_sr = af2m_t10_acc["pdb"].unique().shape[0]/lpdbs*100
af2m_t10_med_sr = af2m_t10_med["pdb"].unique().shape[0]/lpdbs*100
af2m_t10_high_sr = af2m_t10_high["pdb"].unique().shape[0]/lpdbs*100
af_multimer_rank10_list = [af2m_t10_acc_sr, af2m_t10_med_sr, af2m_t10_high_sr]

# let's print all the af_multimer bars
print(f"AF2M acc {af_multimer_rank1_list[0]:.2f} t10 {af_multimer_rank10_list[0]:.2f}")
print(f"AF2M med {af_multimer_rank1_list[1]:.2f} t10 {af_multimer_rank10_list[1]:.2f}")
print(f"AF2M high {af_multimer_rank1_list[2]:.2f} t10 {af_multimer_rank10_list[2]:.2f}")
# now alphafold3
lpdbs = alphafold3["pdb"].unique().shape[0]
af3_t1_acc = alphafold3[(alphafold3["rank"]==1)&(alphafold3["capri"]!="incorrect")]
af3_t1_med = alphafold3[(alphafold3["rank"]==1)&(alphafold3["capri"]!="incorrect")&(alphafold3["capri"]!="acceptable")]
af3_t1_high = alphafold3[(alphafold3["rank"]==1)&(alphafold3["capri"]=="high")]


af3_t1_acc_sr = round(af3_t1_acc["pdb"].unique().shape[0]/lpdbs*100, 2)
af3_t1_med_sr = round(af3_t1_med["pdb"].unique().shape[0]/lpdbs*100, 2)
af3_t1_high_sr = round(af3_t1_high["pdb"].unique().shape[0]/lpdbs*100, 2)
af3_rank1_list = [af3_t1_acc_sr, af3_t1_med_sr, af3_t1_high_sr]

af3_t10_acc = alphafold3[(alphafold3["rank"]<11)&(alphafold3["capri"]!="incorrect")]
af3_t10_med = alphafold3[(alphafold3["rank"]<11)&(alphafold3["capri"]!="incorrect")&(alphafold3["capri"]!="acceptable")]
af3_t10_high = alphafold3[(alphafold3["rank"]<11)&(alphafold3["capri"]=="high")]

af3_t10_acc_sr = round(af3_t10_acc["pdb"].unique().shape[0]/lpdbs*100, 2)
af3_t10_med_sr = round(af3_t10_med["pdb"].unique().shape[0]/lpdbs*100, 2)
af3_t10_high_sr = round(af3_t10_high["pdb"].unique().shape[0]/lpdbs*100, 2)
af3_rank10_list = [af3_t10_acc_sr, af3_t10_med_sr, af3_t10_high_sr]
# same as before, let's print them
print(f"AF3 acc {af3_rank1_list[0]:.1f} t10 {af3_rank10_list[0]:.1f}")
print(f"AF3 med {af3_rank1_list[1]:.1f} t10 {af3_rank10_list[1]:.1f}")
print(f"AF3 high {af3_rank1_list[2]:.1f} t10 {af3_rank10_list[2]:.1f}")

# we add for all three plots the AFMultimer and AF3 results for top1 and top10
bars_af2 = [0.4, 1.4]
bars_af3 = [0.8, 1.8]
for i in range(3):
    axs[i,4].set_ylim(0, 105)
    axs[i,4].set_yticks([])
    if i != 2:
        axs[i,4].set_xticks([])
    else:
        axs[i,4].set_xticks(bars_af2 + bars_af3)
        axs[i,4].set_xticklabels(af3_xlabels, rotation = 45, ha = "right", size = 12)
    axs[i,4].bar(bars_af2, [af_multimer_rank1_list[0], af_multimer_rank10_list[0]], color = "khaki", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Acceptable")
    axs[i,4].bar(bars_af2, [af_multimer_rank1_list[1], af_multimer_rank10_list[1]], color = "lightseagreen", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Medium")
    axs[i,4].bar(bars_af2, [af_multimer_rank1_list[2], af_multimer_rank10_list[2]], color = "darkblue", alpha = 1, edgecolor = "black", width=barwidth, label = "AF High")

    axs[i,4].bar(bars_af3, [af3_rank1_list[0], af3_rank10_list[0]], color = "khaki", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Acceptable")
    axs[i,4].bar(bars_af3, [af3_rank1_list[1], af3_rank10_list[1]], color = "lightseagreen", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Medium")
    axs[i,4].bar(bars_af3, [af3_rank1_list[2], af3_rank10_list[2]], color = "darkblue", alpha = 1, edgecolor = "black", width=barwidth, label = "AF High")

    #and we plot the sum above the bars
    axs[i,4].text(0.4, af_multimer_rank1_list[0], f"{af_multimer_rank1_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
    axs[i,4].text(0.8, af3_rank1_list[0], f"{af3_rank1_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
    axs[i,4].text(1.4, af_multimer_rank10_list[0], f"{af_multimer_rank10_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
    axs[i,4].text(1.8, af3_rank10_list[0], f"{af3_rank10_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
# 
# 
# remove spaces between subplots
plt.subplots_adjust(wspace=0.1, hspace=0.1)
# add a legend that takes labels from axs[2,2] and axs[2,4]
handles, labels = axs[2,2].get_legend_handles_labels()
handles_af, labels_af = axs[2,4].get_legend_handles_labels()
# now use these handles to create a legend for the whole figure
#
axs[0,0].set_title("IB\nensemble", fontsize = 15)
axs[0,1].set_title("IBMo\nensemble", fontsize = 15)
axs[0,2].set_title("IBMu\nensemble", fontsize = 15)
axs[0,3].set_title("IBMuMo\nensemble", fontsize = 15)
axs[0,4].set_title("Alphafold\nend-to-end", fontsize = 15)

axs[0,0].text(-1.5, 55, "True\nInterface", ha = "center", va = "center", fontsize = 15)
axs[1,0].text(-1.5, 55, "Loose\nInterface", ha = "center", va = "center", fontsize = 15)
axs[2,0].text(-1.5, 55, "Two-Hit\nInterface", ha = "center", va = "center", fontsize = 15)


# add label b) to ax[0]
axs[0,0].text(-1.5, 105, "a)", fontsize = 25)

# do tight layout but leave white space at the bottom of the plot
plt.tight_layout(h_pad=0.5, w_pad=-1.5)
plt.subplots_adjust(bottom=0.11)
plt.legend(handles[:3]+handles_af[:3], labels[:3]+labels_af[:3], loc = "lower center", fontsize = 15, ncol = 6, bbox_to_anchor=(-2.75, -0.43))
plt.savefig(Path("figures", "haddock_docking_flexref_sr_per_scenario.png"), dpi=400)

# now let's do the same but without the second and fourth ensemble
fig,axs = plt.subplots(3,4, figsize = (12, 12), width_ratios=[3, 3, 4, 2])
red_ensembles = ["IB", "IBMu"]
for i, ens in enumerate(red_ensembles):
    for j, scen in enumerate(scenarios):
        axs[j,i].set_xlim(0, 3)
        axs[j,i].set_ylim(0, 105)
        if j != 2:
            axs[j,i].set_xticks([])
        else:
            axs[j,i].set_xticks(bars_b + bars_u)
            axs[j,i].set_xticklabels(xlabels, rotation = 45, ha = "right", size = 12)
        if i != 0:
            axs[j,i].set_yticks([])
        else:
            axs[j,i].set_ylabel("SR (%)", size = 15)
        ranks = [1, 10, 200]
        axs[j,i].bar(bars_b, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="b")]["acc_sr"].values[0]*100 for rank in ranks], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
        axs[j,i].bar(bars_b, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="b")]["med_sr"].values[0]*100 for rank in ranks], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
        axs[j,i].bar(bars_b, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="b")]["high_sr"].values[0]*100 for rank in ranks], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
        axs[j,i].bar(bars_u, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="u")]["acc_sr"].values[0]*100 for rank in ranks], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "Acceptable")
        axs[j,i].bar(bars_u, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="u")]["med_sr"].values[0]*100 for rank in ranks], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "Medium")
        axs[j,i].bar(bars_u, [sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==rank)&(sr_data_flex["struct"]=="u")]["high_sr"].values[0]*100 for rank in ranks], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "High")
        # let's label the acceptable SR
        for b_idx, b in enumerate(bars_b):
            b_value = sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==ranks[b_idx])&(sr_data_flex["struct"]=="b")]["acc_sr"].values[0]*100
            u_value = sr_data_flex[(sr_data_flex["ensemble"]==ens)&(sr_data_flex["scenario"]==scen)&(sr_data_flex["rank"]==ranks[b_idx])&(sr_data_flex["struct"]=="u")]["acc_sr"].values[0]*100
            axs[j,i].text(b, b_value, f"{b_value:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize)
            axs[j,i].text(bars_u[b_idx], u_value, f"{u_value:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize)
        
# we add for all three plots the AFMultimer and AF3 results for top1 and top10
bars_af2 = [0.3, 1.3]
bars_af3 = [0.7, 1.7]
for i in range(3):
    axs[i,3].set_xlim(0, 2)
    axs[i,3].set_ylim(0, 105)
    axs[i,3].set_yticks([])
    if i != 2:
        axs[i,3].set_xticks([])
    else:
        axs[i,3].set_xticks(bars_af2 + bars_af3)
        axs[i,3].set_xticklabels(af3_xlabels, rotation = 45, ha = "right", size = 12)
    axs[i,3].bar(bars_af2, [af_multimer_rank1_list[0], af_multimer_rank10_list[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Acceptable")
    axs[i,3].bar(bars_af2, [af_multimer_rank1_list[1], af_multimer_rank10_list[1]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Medium")
    axs[i,3].bar(bars_af2, [af_multimer_rank1_list[2], af_multimer_rank10_list[2]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "AF High")

    axs[i,3].bar(bars_af3, [af3_rank1_list[0], af3_rank10_list[0]], color = "lightblue", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Acceptable")
    axs[i,3].bar(bars_af3, [af3_rank1_list[1], af3_rank10_list[1]], color = "lightgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "AF Medium")
    axs[i,3].bar(bars_af3, [af3_rank1_list[2], af3_rank10_list[2]], color = "darkgreen", alpha = 1, edgecolor = "black", width=barwidth, label = "AF High")

    #and we plot the sum above the bars
    axs[i,3].text(0.3, af_multimer_rank1_list[0], f"{af_multimer_rank1_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
    axs[i,3].text(0.7, af3_rank1_list[0], f"{af3_rank1_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
    axs[i,3].text(1.3, af_multimer_rank10_list[0], f"{af_multimer_rank10_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')
    axs[i,3].text(1.7, af3_rank10_list[0], f"{af3_rank10_list[0]:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize, horizontalalignment='center', verticalalignment='center')

# remove spaces between subplots
plt.subplots_adjust(wspace=0.1, hspace=0.1)
# add a legend that takes labels from axs[2,2] and axs[2,4]
handles, labels = axs[2,0].get_legend_handles_labels()
handles_af, labels_af = axs[2,2].get_legend_handles_labels()
# now use these handles to create a legend for the whole figure
#
axs[0,0].set_title("IB\nensemble", fontsize = 15)
axs[0,1].set_title("IBMu\nensemble", fontsize = 15)
axs[0,2].set_title("IBMu ensemble\n(cluster-based)", fontsize = 15)
axs[0,3].set_title("Alphafold\nend-to-end", fontsize = 15)

axs[0,0].text(-1.5, 55, "True\nInterface", ha = "center", va = "center", fontsize = 15)
axs[1,0].text(-1.5, 55, "Loose\nInterface", ha = "center", va = "center", fontsize = 15)
axs[2,0].text(-1.5, 55, "Two-Hit\nInterface", ha = "center", va = "center", fontsize = 15)

cluster_sr_data = pd.read_csv(Path("..", "data", "haddock_docking_clustered_sr_per_scenario.tsv"), sep = "\t")
# 
cluster_sr_data_ibmu = cluster_sr_data[cluster_sr_data["ensemble"]=="IBMu"]
#
scenarios = ["real", "loose", "twohit"]
# n_array = [1, 2, 3, 4, 5, 10]
n_array = [1,2,5,10]
acc_labels = ["Acceptable", "Medium", "High"]
acc_colors = ["lightblue", "lightgreen", "darkgreen"]
scenarios_titles = ["True Interface", "Loose Interface", "Two-Hit Interface"]
bound_bars = [0.3, 1.3, 2.3, 3.3]
unb_bars = [el + 0.4 for el in bound_bars]
ticks = [el for el in bound_bars] + [el + 0.4 for el in bound_bars]
tick_labels = ["B T1", "B T2", "B T5", "B T10", "U T1", "U T2", "U T5", "U T10"]
for i, scen in enumerate(scenarios):
    # features of the plot
    axs[i,2].set_ylim(0, 105)
    axs[i,2].set_xlim(0, 4)
    axs[i,2].set_yticks([])
    #axs[i,2].set_xticklabels(tick_labels, size = 10, rotation = 45, ha = "right")
    if i != 2:
        axs[i,2].set_xticks([])
    else:
        axs[i,2].set_xticks(ticks)
        axs[i,2].set_xticklabels(tick_labels, size = 12, rotation = 45, ha = "right")

    # we plot the bars
    df_scen = cluster_sr_data_ibmu[cluster_sr_data_ibmu["scenario"]==scen]
    df_scen_u = df_scen[df_scen["struct"]=="u"]
    df_scen_b = df_scen[df_scen["struct"]=="b"]
    for j, sr in enumerate(["acc_sr", "med_sr", "high_sr"]):
        vals_bound = [df_scen_b[df_scen_b["rank"]==n][sr].values[0]*100 for n in n_array]
        vals_unbound = [df_scen_u[df_scen_u["rank"]==n][sr].values[0]*100 for n in n_array]
        axs[i,2].bar(bound_bars, vals_bound, color = acc_colors[j], alpha = 1, edgecolor = "black", width=barwidth, label = acc_labels[j])
        axs[i,2].bar(unb_bars, vals_unbound, color = acc_colors[j], alpha = 1, edgecolor = "black", width=barwidth, label = acc_labels[j])
        if j == 0: # acc
            for bar_idx, n in enumerate(n_array):
                value_bound = vals_bound[bar_idx]
                value_unbound = vals_unbound[bar_idx]
                axs[i,2].text(bound_bars[bar_idx], vals_bound[bar_idx], f"{value_bound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize) 
                axs[i,2].text(unb_bars[bar_idx], vals_unbound[bar_idx], f"{value_unbound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize)

## do tight layout but leave white space at the bottom of the plot
plt.tight_layout(h_pad=0.5, w_pad=-1.5)
plt.subplots_adjust(bottom=0.12)
plt.legend(handles[:3], labels[:3], loc = "lower center", fontsize = 15, ncol = 3, bbox_to_anchor=(-2.2, -0.425))
plt.savefig(Path("figures", "figure2.png"), dpi=400)

# given sr_data, extract the top10 SR for the mix-loose scenario
mix_loose_sr_u = sr_data_flex[(sr_data_flex["scenario"]=="mix-loose")&(sr_data_flex["struct"]=="u")&(sr_data_flex["rank"]==10)]
mix_loose_sr_b = sr_data_flex[(sr_data_flex["scenario"]=="mix-loose")&(sr_data_flex["struct"]=="b")&(sr_data_flex["rank"]==10)]
print(f"Mix-loose T10 SR unbound: {mix_loose_sr_u['acc_sr'].values[0]*100:.2f}")
print(f"Mix-loose T10 SR bound: {mix_loose_sr_b['acc_sr'].values[0]*100:.2f}")
mix_twohit_sr_u = sr_data_flex[(sr_data_flex["scenario"]=="mix-twohit")&(sr_data_flex["struct"]=="u")&(sr_data_flex["rank"]==10)]
mix_twohit_sr_b = sr_data_flex[(sr_data_flex["scenario"]=="mix-twohit")&(sr_data_flex["struct"]=="b")&(sr_data_flex["rank"]==10)]
print(f"Mix-twohit T10 SR unbound: {mix_twohit_sr_u['acc_sr'].values[0]*100:.2f}")
print(f"Mix-twohit T10 SR bound: {mix_twohit_sr_b['acc_sr'].values[0]*100:.2f}")

#SUPPLEMENTARY FIGURE 2
#Bar plots showing docking success rate (SR) for All Surface scenario as a function of the TopN ranked structures.
fig,axs = plt.subplots(1,3, figsize = (12, 6), width_ratios=[3, 3, 3])
all_surf_data = sr_data[sr_data["scenario"]=="surf"]
# 
barwidth = 0.7
# 
axs[0].set_xlim(0, 3)
bound_bars = [0.5, 1.5, 2.5]
n_array = [1,10,200]
surf_rigid = all_surf_data[all_surf_data["stage"]=="rigid"]
print(surf_rigid)
for j, sr in enumerate(["acc_sr", "med_sr", "high_sr"]):
    vals_bound = [surf_rigid[surf_rigid["rank"]==n][sr].values[0]*100 for n in n_array]
    axs[0].bar(bound_bars, vals_bound, color = acc_colors[j], alpha = 1, edgecolor = "black", width=barwidth, label = acc_labels[j])
    if j == 0: # acc
        for bar_idx, n in enumerate(n_array):
            value_bound = vals_bound[bar_idx]
            axs[0].text(bound_bars[bar_idx], vals_bound[bar_idx], f"{value_bound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = 12) 
surf_flex = all_surf_data[all_surf_data["stage"]=="flex"]
for j, sr in enumerate(["acc_sr", "med_sr", "high_sr"]):
    vals_bound = [surf_flex[surf_flex["rank"]==n][sr].values[0]*100 for n in n_array]
    axs[1].bar(bound_bars, vals_bound, color = acc_colors[j], alpha = 1, edgecolor = "black", width=barwidth, label = acc_labels[j])
    if j == 0: # acc
        for bar_idx, n in enumerate(n_array):
            value_bound = vals_bound[bar_idx]
            axs[1].text(bound_bars[bar_idx], vals_bound[bar_idx], f"{value_bound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = 12) 
surf_em = all_surf_data[all_surf_data["stage"]=="emref"]
for j, sr in enumerate(["acc_sr", "med_sr", "high_sr"]):
    vals_bound = [surf_em[surf_em["rank"]==n][sr].values[0]*100 for n in n_array]
    axs[2].bar(bound_bars, vals_bound, color = acc_colors[j], alpha = 1, edgecolor = "black", width=barwidth, label = acc_labels[j])
    if j == 0: # acc
        for bar_idx, n in enumerate(n_array):
            value_bound = vals_bound[bar_idx]
            axs[2].text(bound_bars[bar_idx], vals_bound[bar_idx], f"{value_bound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = 12) 
for i in range(3):
    axs[i].set_ylim(0, 100)
    axs[i].set_xticks(bound_bars)
    axs[i].set_xticklabels(["T1", "T10", "T200"], size = 15)

#we plot the name of the stages
axs[0].set_title("Rigidbody stage", fontsize = 15)
axs[1].set_title("Flexref stage", fontsize = 15)
axs[2].set_title("Emref stage", fontsize = 15)
axs[0].set_ylabel("SR (%)", size = 15)
# empty y tick labels on axs[1] and axs[2]
axs[1].set_yticks([])
axs[2].set_yticks([])
axs[0].set_yticks([0, 20, 40, 60, 80, 100])
axs[0].set_yticklabels([0, 20, 40, 60, 80, 100], size = 12)

handles, labels = axs[0].get_legend_handles_labels()
plt.tight_layout()
plt.subplots_adjust(bottom=0.14)
plt.legend(handles, labels, loc = "lower center", fontsize = 15, ncol = 3, bbox_to_anchor=(-0.57, -0.18))
plt.savefig(Path("figures", "SI_figure2.png"), dpi=400)

# ##SUPPLEMENTARY FIGURE 3A
# ##Bar plot comparing model docking success rates (SRs) for TopN structures after HADDOCK scoring (HS) and Voronoi scoring (VS) for different information scenarios.

fig,axs = plt.subplots(1,3, figsize = (12, 6), width_ratios=[3, 3, 3])
barwidth = 0.4
xlabels = ["VS T1", "VS T10", "VS T200", "HS T1", "HS T10", "HS T200"]
sr_data_voro = sr_data[sr_data["stage"]=="emref-voro"]
sr_data_haddock_emref = sr_data[(sr_data["stage"]=="emref")&(sr_data["ensemble"]=="IBMu")&(sr_data["struct"]=="u")]

for j, scen in enumerate(scenarios):
    sr_data_voro_scen = sr_data_voro[sr_data_voro["scenario"]==scen]
    sr_data_haddock_emref_scen = sr_data_haddock_emref[sr_data_haddock_emref["scenario"]==scen]
    axs[j].set_title(scenarios_titles[j], fontsize = 15)
    axs[j].set_xlim(0, 3)
    axs[j].set_ylim(0, 105)
    axs[j].set_xticks([])
    axs[j].set_xticks(bars_b + bars_u)
    axs[j].set_xticklabels(xlabels, rotation = 45, ha = "right", size = 12)
    if j != 0:
        axs[j].set_yticks([])
    else:
        axs[j].set_ylabel("SR (%)", size = 15)
    ranks = [1, 10, 200]
    for i, sr in enumerate(["acc_sr", "med_sr", "high_sr"]):
        voro_acc = [sr_data_voro_scen[sr_data_voro_scen["rank"]==rank][sr].values[0]*100 for rank in ranks]
        haddock_acc = [sr_data_haddock_emref_scen[sr_data_haddock_emref_scen["rank"]==rank][sr].values[0]*100 for rank in ranks]
        print(f"voro {sr} ({scen}): {voro_acc}")
        axs[j].bar(bars_b, voro_acc, color = acc_colors[i], alpha = 1, edgecolor = "black", width=barwidth, label = acc_labels[i])
        axs[j].bar(bars_u, haddock_acc, color = acc_colors[i], alpha = 1, edgecolor = "black", width=barwidth)
        # let's label the acceptable SR
        if i == 0: # acc
            for bar_idx, rank in enumerate(ranks):
                value_bound = voro_acc[bar_idx]
                value_unbound = haddock_acc[bar_idx]
                axs[j].text(bars_b[bar_idx], voro_acc[bar_idx], f"{value_bound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize) 
                axs[j].text(bars_u[bar_idx], haddock_acc[bar_idx], f"{value_unbound:.1f}", ha = "center", va = "bottom", rotation = 0, fontsize = text_fontsize)
    
handles, labels = axs[0].get_legend_handles_labels()
plt.tight_layout()
plt.subplots_adjust(bottom=0.20)
plt.legend(handles, labels, loc = "lower center", fontsize = 15, ncol = 3, bbox_to_anchor=(-0.57, -0.28)) 
plt.savefig(Path("figures", "SI_figure3.png"), dpi=400)