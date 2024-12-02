###Representation of the accuracy of the unbound nanobody and antigen structure ensembles for the docking (Table1, Figure 2 and Supplementary Figure 1)

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


#we load the RMSD calculations
af_multimer = pd.read_csv(Path("..","data","alphafold2multimer_rmsd.tsv"), sep="\t")
af_monomer = pd.read_csv(Path("..","data","alphafold2_rmsd.tsv"), sep="\t")
immunebuilder = pd.read_csv(Path("..","data","immunebuilder_rmsd.tsv"), sep="\t")
raptorx = pd.read_csv(Path("..","data","raptorxsingle_rmsd.tsv"), sep="\t")
nanonet = pd.read_csv(Path("..","data","nanonet_rmsd.tsv"), sep="\t")
af3 = pd.read_csv(Path("..","data","alphafold3_rmsd.tsv"), sep="\t")

##TABLE 1
##Backbone mean RMSD (± standard deviation) of the different Nb regions by the different prediction methods

##We get the top1 ranked for all but RaptorX-Single and NanoNet
af_multimer_top1 = af_multimer[af_multimer['rank']==1]
af_monomer_top1 = af_monomer[af_monomer['rank']==1]
immunebuilder_top1 = immunebuilder[immunebuilder['rank']==1]
af3_top1 = af3[af3['rank']==1]

##We create a .csv with all 6 methods (top1 for AlphaFold2, AlphaFold2-Multimer and ImmuneBuilder)
##Columns are FR_RMS mean, FR_RMS std, FR_RMS median, CDR1_RMS mean, CDR1_RMS std, CDR1_RMS median, CDR2_RMS mean, CDR2_RMS std, CDR2_RMS median, CDR3_RMS mean, CDR3_RMS std, CDR3_RMS median
rms_csv = []
for methods in [af_multimer_top1, af_monomer_top1, immunebuilder_top1, raptorx, nanonet, af3]:
    FR_mean = round(methods['FR_RMS'].mean(), 2)
    FR_std = round(methods['FR_RMS'].std(), 2)
    FR_median = round(methods['FR_RMS'].median(), 2)
    CDR1_mean = round(methods['CDR1_RMS'].mean(), 2)
    CDR1_std = round(methods['CDR1_RMS'].std(), 2)
    CDR1_median = round(methods['CDR1_RMS'].median(), 2)
    CDR2_mean = round(methods['CDR2_RMS'].mean(), 2)
    CDR2_std = round(methods['CDR2_RMS'].std(), 2)
    CDR2_median = round(methods['CDR2_RMS'].median(), 2)
    CDR3_mean = round(methods['CDR3_RMS'].mean(), 2)
    CDR3_std = round(methods['CDR3_RMS'].std(), 2)
    CDR3_median = round(methods['CDR3_RMS'].median(), 2)

    rms_csv.append([FR_mean, FR_std, FR_median, CDR1_mean, CDR1_std, CDR1_median, CDR2_mean, CDR2_std, CDR2_median, CDR3_mean, CDR3_std, CDR3_median])

rms_df = pd.DataFrame(rms_csv, columns=['FR_RMS mean', 'FR_RMS std', 'FR_RMS median', 'CDR1_RMS mean', 'CDR1_RMS std', 'CDR1_RMS median', 'CDR2_RMS mean', 'CDR2_RMS std', 'CDR2_RMS median', 'CDR3_RMS mean', 'CDR3_RMS std', 'CDR3_RMS median'],
                      index=['af_multimer top1', 'af_monomer top1', 'immunebuilder top1', 'raptorx', 'nanonet', 'af3'])

rms_df.to_csv(Path(".", "figures", "nanobody_rmsd_summary.tsv"), sep="\t")


##SUPPLEMENTARY FIGURE 1
##Violin plot showing the distribution of the backbone CDR3 RMSD from the best ranked nanobody prediction for the different methods.

##We get the top1 ranked for all but RaptorX-Single and NanoNet
af_multimer_top1 = af_multimer[af_multimer['rank']==1]
af_monomer_top1 = af_monomer[af_monomer['rank']==1]
immunebuilder_top1 = immunebuilder[immunebuilder['rank']==1]
af3_top1 = af3[af3['rank']==1]

fig, ax = plt.subplots(figsize=(10, 5))

violin = plt.violinplot([af_multimer_top1['CDR3_RMS'], af_monomer_top1['CDR3_RMS'], immunebuilder_top1['CDR3_RMS'], raptorx['CDR3_RMS'], nanonet['CDR3_RMS'], af3_top1['CDR3_RMS']],
 showmeans=False, showmedians=True, showextrema=False)

#we set a different color for each method
colors = ['darkblue', 'lightblue', 'red', 'pink', 'orange', 'purple']
for i in range(len(violin['bodies'])):
    violin['bodies'][i].set_facecolor(colors[i])
    violin['bodies'][i].set_edgecolor('black')
    violin['bodies'][i].set_linewidth(1)

#we plot the median value
ax.text(1, af_multimer_top1['CDR3_RMS'].median() + 0.1, str(round(af_multimer_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(2, af_monomer_top1['CDR3_RMS'].median() + 0.1, str(round(af_monomer_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(3, immunebuilder_top1['CDR3_RMS'].median() + 0.1, str(round(immunebuilder_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(4, raptorx['CDR3_RMS'].median() + 0.1, str(round(raptorx['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(5, nanonet['CDR3_RMS'].median() + 0.1, str(round(nanonet['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(6, af3_top1['CDR3_RMS'].median() + 0.1, str(round(af3_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)

plt.xticks([1, 2, 3, 4, 5, 6], ['AF Multimer\n(Top1 models)', 'AF Monomer\n(Top1 models)', 'IB\n(Top1 models)', 'RXS', 'NN', 'AF3\n(Top1 models)'], size=10)
plt.ylabel('CDR3 RMSD (Å)', size = "medium")
plt.tight_layout()
plt.savefig(Path(".", "figures", "cdr3_rmsd_top1_prediction_comparison.png"))


#FIGURE 2A
#Violin plot showing the distribution of the best backbone CDR3 RMSD in the ensembles.

#we get the best available prediction for each pdb in the ensembles
best_ib = immunebuilder.groupby('pdb').min().reset_index()
best_afm = af_monomer.groupby('pdb').min().reset_index()
best_afmm = af_multimer.groupby('pdb').min().reset_index()
af3_ranked_df = af3.groupby('pdb').min().reset_index()

#now we combine the different ensembles
all_rmsd_df = pd.concat([immunebuilder, af_multimer, af_monomer], axis=0)
ib_mu_df = pd.concat([immunebuilder, af_multimer], axis=0)
ib_mo_df = pd.concat([immunebuilder, af_monomer], axis=0)
#and we also get the best available prediction for each pdb in the combined ensembles
all_rmsd_df = all_rmsd_df.groupby('pdb').min().reset_index()
ib_mu_df = ib_mu_df.groupby('pdb').min().reset_index()
ib_mo_df = ib_mo_df.groupby('pdb').min().reset_index()


#we read the data for the clustered ensembles
ib_mo_clustered = pd.read_csv(Path("..","data","alphahold2_monomer_immunebuilder_centroids_rmsd.tsv"), sep="\t")
ib_mu_clustered = pd.read_csv(Path("..","data","alphafold2multimer_immunebuilder_centroids_rmsd.tsv"), sep="\t")
all_clustered = pd.read_csv(Path("..","data","alphafold2multimer_monomer_immunebuilder_centroids_rmsd.tsv"), sep="\t")
#and we also group them by pdb and get the minimum CDR3 RMSD
all_clustered = all_clustered.groupby('pdb').min().reset_index()
ib_mo_clustered = ib_mo_clustered.groupby('pdb').min().reset_index()
ib_mu_clustered = ib_mu_clustered.groupby('pdb').min().reset_index()

##we print mean ± std and median of CDR3 RMSD for each combination
print("IB ensemble:")
print("Mean ± std: ", best_ib['CDR3_RMS'].mean(), " ± ", best_ib['CDR3_RMS'].std())
print("Median: ", best_ib['CDR3_RMS'].median())
print("\nAF Monomer ensemble:")
print("Mean ± std: ", best_afm['CDR3_RMS'].mean(), " ± ", best_afm['CDR3_RMS'].std())
print("Median: ", best_afm['CDR3_RMS'].median())
print("\nAF Multimer ensemble:")
print("Mean ± std: ", best_afmm['CDR3_RMS'].mean(), " ± ", best_afmm['CDR3_RMS'].std())
print("Median: ", best_afmm['CDR3_RMS'].median())
print("\nAF3 ensemble:")
print("Mean ± std: ", af3_ranked_df['CDR3_RMS'].mean(), " ± ", af3_ranked_df['CDR3_RMS'].std())
print("Median: ", af3_ranked_df['CDR3_RMS'].median())

print("\n\nIB + AF Monomer + AF Multimer:")
print("Mean ± std: ", all_rmsd_df['CDR3_RMS'].mean(), " ± ", all_rmsd_df['CDR3_RMS'].std())
print("Median: ", all_rmsd_df['CDR3_RMS'].median())
print("\nIB + AF Monomer:")
print("Mean ± std: ", ib_mo_df['CDR3_RMS'].mean(), " ± ", ib_mo_df['CDR3_RMS'].std())
print("Median: ", ib_mo_df['CDR3_RMS'].median())
print("\nIB + AF Multimer:")
print("Mean ± std: ", ib_mu_df['CDR3_RMS'].mean(), " ± ", ib_mu_df['CDR3_RMS'].std())
print("Median: ", ib_mu_df['CDR3_RMS'].median())

print("\n\nIB + AF Monomer + AF Multimer clustered:")
print("Mean ± std: ", all_clustered['CDR3_RMS'].mean(), " ± ", all_clustered['CDR3_RMS'].std())
print("Median: ", all_clustered['CDR3_RMS'].median())
print("\nIB + AF Monomer clustered:")
print("Mean ± std: ", ib_mo_clustered['CDR3_RMS'].mean(), " ± ", ib_mo_clustered['CDR3_RMS'].std())
print("Median: ", ib_mo_clustered['CDR3_RMS'].median())
print("\nIB + AF Multimer clustered:")
print("Mean ± std: ", ib_mu_clustered['CDR3_RMS'].mean(), " ± ", ib_mu_clustered['CDR3_RMS'].std())
print("Median: ", ib_mu_clustered['CDR3_RMS'].median())

#we plot the violin plot
fig, ax = plt.subplots(figsize=(10, 5))

violin = plt.violinplot([best_ib['CDR3_RMS'], best_afm['CDR3_RMS'], best_afmm['CDR3_RMS'],
                            ib_mo_df['CDR3_RMS'], ib_mu_df['CDR3_RMS'], all_rmsd_df['CDR3_RMS'],
                            ib_mo_clustered['CDR3_RMS'], ib_mu_clustered['CDR3_RMS'], all_clustered['CDR3_RMS'],
                            af3_ranked_df['CDR3_RMS']],
    showmeans=False, showmedians=True, showextrema=False)

#we set a different color for each method
colors = ['red', 'lightblue', 'darkblue', 'yellow', 'orange', 'brown', 'yellow', 'orange', 'brown', "purple"]
for i in range(len(violin['bodies'])):
    violin['bodies'][i].set_facecolor(colors[i])
    violin['bodies'][i].set_edgecolor('black')
    violin['bodies'][i].set_linewidth(1)

#we plot the median value
ax.text(1, best_ib['CDR3_RMS'].median() + 0.1, str(round(best_ib['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(2, best_afm['CDR3_RMS'].median() + 0.1, str(round(best_afm['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(3, best_afmm['CDR3_RMS'].median() + 0.1, str(round(best_afmm['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(4, ib_mo_df['CDR3_RMS'].median() + 0.1, str(round(ib_mo_df['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(5, ib_mu_df['CDR3_RMS'].median() + 0.1, str(round(ib_mu_df['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(6, all_rmsd_df['CDR3_RMS'].median() + 0.1, str(round(all_rmsd_df['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(7, ib_mo_clustered['CDR3_RMS'].median() + 0.1, str(round(ib_mo_clustered['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(8, ib_mu_clustered['CDR3_RMS'].median() + 0.1, str(round(ib_mu_clustered['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(9, all_clustered['CDR3_RMS'].median() + 0.1, str(round(all_clustered['CDR3_RMS'].median(), 2)), color='black', ha='center')
ax.text(10, af3_ranked_df['CDR3_RMS'].median() + 0.1, str(round(af3_ranked_df['CDR3_RMS'].median(), 2)), color='black', ha='center')

plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ['IB', 'AFMo', 'AFMu', 'IBMo', 'IBMu', 'IBMM', 'IBMo-Cl', 'IBMu-Cl', 'IBMM-Cl', 'AF3'], size=10)
plt.ylabel('CDR3 RMSD (Å)')

plt.tight_layout()
plt.savefig(Path(".", "figures", "cdr3_rmsd_best_ensemble_best_prediction_comparison.png"))


##FIGURE 2B
##Violin plot showing the distribution of the antigen and epitope backbone RMSD

#we read the unbound antigen RMSD data
antigen_rmsd = pd.read_csv(Path("..","data","unbound_antigens_rmsd.tsv"), sep="\t")


fig, ax = plt.subplots(figsize=(7, 5))
medians = [antigen_rmsd["af_multimer_whole_rmsd"].median(), antigen_rmsd["af_multimer_epitope_rmsd"].median(), antigen_rmsd["af_multimer_pdb_whole_rmsd"].median(), antigen_rmsd["af_multimer_pdb_epitope_rmsd"].median()]

violin = plt.violinplot([antigen_rmsd["af_multimer_whole_rmsd"], antigen_rmsd["af_multimer_epitope_rmsd"], antigen_rmsd["af_multimer_pdb_whole_rmsd"], antigen_rmsd["af_multimer_pdb_epitope_rmsd"]], showmeans = False, showmedians = True, showextrema=False)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(["Antigen RMSD", "Epitope RMSD", "Antigen RMSD", "Epitope RMSD"])
plt.ylabel("Backbone RMSD (Å)")
ax.set_ylim(0, 30)

for tick in [1,2,3,4]:
   ax.text(tick, medians[tick-1] + 0.2 , round(medians[tick -1], 2),
            horizontalalignment='center',
            style = 'oblique',
            color='black')
ax.axhline(y=5, color='black', linestyle='--')
ax.text(1, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_whole_rmsd'] > 5])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.text(2, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_epitope_rmsd'] > 5])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.text(3, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_whole_rmsd'] > 5])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.text(4, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_epitope_rmsd'] > 5])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.axhline(y=10, color='black', linestyle='--')
ax.text(1, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_whole_rmsd'] > 10])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.text(2, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_epitope_rmsd'] > 10])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.text(3, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_whole_rmsd'] > 10])} structures", horizontalalignment='center', style = 'oblique', color='black')
ax.text(4, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_epitope_rmsd'] > 10])} structures", horizontalalignment='center', style = 'oblique', color='black')
#we print the max value of antigen and epitope RMSD
ax.text(1, antigen_rmsd["af_multimer_whole_rmsd"].max()+1, f"{round(antigen_rmsd['af_multimer_whole_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black')
ax.text(2, antigen_rmsd["af_multimer_epitope_rmsd"].max()+1, f"{round(antigen_rmsd['af_multimer_epitope_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black')
ax.text(3, antigen_rmsd["af_multimer_pdb_whole_rmsd"].max()+1, f"{round(antigen_rmsd['af_multimer_pdb_whole_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black')
ax.text(4, antigen_rmsd["af_multimer_pdb_epitope_rmsd"].max()+1, f"{round(antigen_rmsd['af_multimer_pdb_epitope_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black')

#we change the color for the two semibound violins to orange
for pc in violin['bodies']:
    pc.set_facecolor('blue')
    pc.set_edgecolor('black')
    pc.set_alpha(0.5)
    pc.set_label('AF Multimer')
violin['bodies'][2].set_facecolor('orange')
violin['bodies'][3].set_facecolor('orange')
violin['bodies'][2].set_label('AF Multimer + PDB')
violin['bodies'][3].set_label('AF Multimer + PDB')

#we plot a legend for the colors
leg = ax.legend(loc='upper right', labels=['AF Multimer', 'AF Multimer + PDB'], fontsize = 12)
leg.legend_handles[0].set_color('blue')
leg.legend_handles[1].set_color('orange')
leg.legend_handles[0].set_edgecolor('black')
leg.legend_handles[1].set_edgecolor('black')
leg.legend_handles[0].set_alpha(0.5)
leg.legend_handles[1].set_alpha(0.5)

plt.tight_layout()
plt.savefig(Path(".", "figures", "rmsd_unbound_antigen_comparison.png"))