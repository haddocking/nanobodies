###Representation of AlphaFold2-Multimer and AlphaFold3 docking success on Nanobody-Antigen complex prediction (Figure 1)
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Patch



##FIGURE 1A
##Bar plot showing acceptable success rates (SRs) for AlphaFold2-Multimer (AF Mu) and AlphaFold3 (AF3) as a function of the Top1,5,10 and 25 ranked models

#we read the data
top_data = pd.read_csv(Path("..", "data/alphafold3_complex_predictions_sr.tsv"), sep = "\t")
top_data_af2 = pd.read_csv(Path("..", "data/alphafold2multimer_complex_predictions_sr.tsv"), sep = "\t")
#get the percentage 
top_data = top_data*100
top_data_af2 = top_data_af2*100

top_data_plot_acceptable_3 = top_data.iloc[:]["acceptable"] - top_data.iloc[:]["medium"]
top_data_plot_medium_3 = top_data.iloc[:]["medium"] - top_data.iloc[:]["high"]
top_data_plot_acceptable_2 = top_data_af2.iloc[:]["acceptable"] - top_data_af2.iloc[:]["medium"]
top_data_plot_medium_2= top_data_af2.iloc[:]["medium"] - top_data_af2.iloc[:]["high"]


fig,ax = plt.subplots(figsize = (8,5))
barWidth = 0.5

#we first plot acceptable, medium and high for top1 Alphafold2
ax.bar(1, top_data_plot_acceptable_2[0], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(1, top_data_plot_medium_2[0], bottom = top_data_plot_acceptable_2[0], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(1, top_data_af2.iloc[:]["high"][0], bottom = top_data_plot_medium_2[0] + top_data_plot_acceptable_2[0], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top1 Alphafold3
ax.bar(1.75, top_data_plot_acceptable_3[0], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(1.75, top_data_plot_medium_3[0], bottom = top_data_plot_acceptable_3[0], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(1.75, top_data.iloc[:]["high"][0], bottom = top_data_plot_medium_3[0] + top_data_plot_acceptable_3[0], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top5 Alphafold2
ax.bar(3, top_data_plot_acceptable_2[4], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(3, top_data_plot_medium_2[4], bottom = top_data_plot_acceptable_2[4], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(3, top_data_af2.iloc[:]["high"][4], bottom = top_data_plot_medium_2[4] + top_data_plot_acceptable_2[4], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top5 Alphafold3
ax.bar(3.75, top_data_plot_acceptable_3[1], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(3.75, top_data_plot_medium_3[1], bottom = top_data_plot_acceptable_3[1], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(3.75, top_data.iloc[:]["high"][1], bottom = top_data_plot_medium_3[1] + top_data_plot_acceptable_3[1], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top10 Alphafold2
ax.bar(5, top_data_plot_acceptable_2[5], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(5, top_data_plot_medium_2[5], bottom = top_data_plot_acceptable_2[5], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(5, top_data_af2.iloc[:]["high"][5], bottom = top_data_plot_medium_2[5] + top_data_plot_acceptable_2[5], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top10 Alphafold3
ax.bar(5.75, top_data_plot_acceptable_3[2], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(5.75, top_data_plot_medium_3[2], bottom = top_data_plot_acceptable_3[2], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(5.75, top_data.iloc[:]["high"][2], bottom = top_data_plot_medium_3[2] + top_data_plot_acceptable_3[2], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top25 Alphafold2
ax.bar(7, top_data_plot_acceptable_2[7], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(7, top_data_plot_medium_2[7], bottom = top_data_plot_acceptable_2[7], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(7, top_data_af2.iloc[:]["high"][7], bottom = top_data_plot_medium_2[7] + top_data_plot_acceptable_2[7], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

#we plot acceptable, medium and high for top25 Alphafold3
ax.bar(7.75, top_data_plot_acceptable_3[3], color = 'lightblue', edgecolor = 'grey', width = barWidth, label = "Acceptable")
ax.bar(7.75, top_data_plot_medium_3[3], bottom = top_data_plot_acceptable_3[3], color = 'palegreen', edgecolor = 'grey', width = barWidth, label = "Medium")
ax.bar(7.75, top_data.iloc[:]["high"][3], bottom = top_data_plot_medium_3[3] + top_data_plot_acceptable_3[3], color = 'darkgreen', edgecolor = 'grey', width = barWidth, label = "High")

ax.set_xticks([1, 1.75, 3, 3.75, 5, 5.75, 7, 7.75])
ax.set_xticklabels(["Top1\n AF Mu", "Top1\n AF3", "Top5\n AF Mu", "Top5\n AF3", "Top10\n AF Mu", "Top10\n AF3", "Top25\n AF Mu", "Top25\n AF3"], rotation = 0, fontsize = "medium")
ax.set_ylabel("SR (%)", fontsize = 12)
ax.set_ylim(0,100)

#we plot the acceptable SR% on top of the bars
ax.text(1, top_data_plot_acceptable_2[0] + top_data_plot_medium_2[0] + top_data_af2.iloc[:]["high"][0] + 1, round(top_data_plot_acceptable_2[0] + top_data_plot_medium_2[0] + top_data_af2.iloc[:]["high"][0], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(1.75, top_data_plot_acceptable_3[0] + top_data_plot_medium_3[0] + top_data.iloc[:]["high"][0] + 1, round(top_data_plot_acceptable_3[0] + top_data_plot_medium_3[0] + top_data.iloc[:]["high"][0], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(3, top_data_plot_acceptable_2[4] + top_data_plot_medium_2[4] + top_data_af2.iloc[:]["high"][4] + 1, round(top_data_plot_acceptable_2[4] + top_data_plot_medium_2[4] + top_data_af2.iloc[:]["high"][4], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(3.75, top_data_plot_acceptable_3[1] + top_data_plot_medium_3[1] + top_data.iloc[:]["high"][1] + 1, round(top_data_plot_acceptable_3[1] + top_data_plot_medium_3[1] + top_data.iloc[:]["high"][1], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(5, top_data_plot_acceptable_2[5] + top_data_plot_medium_2[5] + top_data_af2.iloc[:]["high"][5] + 1, round(top_data_plot_acceptable_2[5] + top_data_plot_medium_2[5] + top_data_af2.iloc[:]["high"][5], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(5.75, top_data_plot_acceptable_3[2] + top_data_plot_medium_3[2] + top_data.iloc[:]["high"][2] + 1, round(top_data_plot_acceptable_3[2] + top_data_plot_medium_3[2] + top_data.iloc[:]["high"][2], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(7, top_data_plot_acceptable_2[7] + top_data_plot_medium_2[7] + top_data_af2.iloc[:]["high"][7] + 1, round(top_data_plot_acceptable_2[7] + top_data_plot_medium_2[7] + top_data_af2.iloc[:]["high"][7], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)
ax.text(7.75, top_data_plot_acceptable_3[3] + top_data_plot_medium_3[3] + top_data.iloc[:]["high"][3] + 1, round(top_data_plot_acceptable_3[3] + top_data_plot_medium_3[3] + top_data.iloc[:]["high"][3], 2), horizontalalignment='center', style = 'oblique', color='black',size = 10)

#We plot a legend for the colors
legend_elements = [Patch(facecolor='lightblue', edgecolor='grey', label='Acceptable'), Patch(facecolor='palegreen', edgecolor='grey', label='Medium'), Patch(facecolor='darkgreen', edgecolor='grey', label='High')]
ax.legend(handles=legend_elements, loc = "upper left", fontsize = "medium")

plt.tight_layout()
plt.savefig(Path(".","figures", "alphafold2multimer_af3_sr_comparison.png"))


#FIGURE 1B
#Violin plot showing the DockQ distribution of AlphaFold2-Multimer (AF Multimer) and AlphaFold3 (AF3) Top1-ranked models.

#we read the dockQ data
afmultimer = pd.read_csv(Path("..", "data","alphafold2multimer_complex_predictions_dockq_top1.tsv"), sep = "\t")
af3 = pd.read_csv(Path("..", "data","alphafold3_complex_predictions_dockq_top1.tsv"), sep = "\t")

fig, ax = plt.subplots(figsize = (6,5))
violin = plt.violinplot([afmultimer["dockq"], af3["dockq"]], showmeans = True, showmedians = False, showextrema = False)

#we plot the value of the median just above the line
ax.text(1, afmultimer["dockq"].mean() + 0.01, round(afmultimer["dockq"].mean(),3), horizontalalignment='center', style = 'oblique', color='black', size = 'small')
ax.text(2, af3["dockq"].mean() + 0.01, round(af3["dockq"].mean(),3), horizontalalignment='center', style = 'oblique', color='black', size = 'small')

ax.set_xticks([1,2])
ax.set_xticklabels(["AF Multimer (Top1)", "AF3 (Top1)"], fontsize = 12)
ax.set_ylabel("DockQ", fontsize = 12)
ax.set_ylim(-0.03,1)

#we draw a --- line at 0.23
plt.axhline(y=0.23, color='r', linestyle='--')
plt.axhline(y=0.6, color='green', linestyle='--')

#we plot the amount of structures with dockQ < 0.23 in red, and the amount of structures with dockQ > 0.6 in green
ax.text(1, 0.1, f"{afmultimer[afmultimer['dockq']<0.23].shape[0]} structures", horizontalalignment='center', style = 'oblique', color='red', size = 'medium')
ax.text(2, 0.1, f"{af3[af3['dockq']<0.23].shape[0]} structures", horizontalalignment='center', style = 'oblique', color='red', size = 'medium')
ax.text(1, 0.65, f"{afmultimer[afmultimer['dockq']>0.6].shape[0]} structures", horizontalalignment='center', style = 'oblique', color='green', size = 'medium')
ax.text(2, 0.65, f"{af3[af3['dockq']>0.6].shape[0]} structures", horizontalalignment='center', style = 'oblique', color='green', size = 'medium')

plt.tight_layout()
plt.savefig(Path(".", "figures", "alphafold2multimer_af3_dockq_top1_comparison.png"))

