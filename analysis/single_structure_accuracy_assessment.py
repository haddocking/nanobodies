###Representation of the accuracy of the unbound nanobody and antigen structure ensembles for the docking (Table1, Figure 2 and Supplementary Figure 1)

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


#we load the RMSD calculations
af_multimer = pd.read_csv(Path("..","data","alphafold2_multimer_plddt_rmsd.tsv"), sep="\t")
af_monomer = pd.read_csv(Path("..","data","alphafold2_monomer_plddt_rmsd.tsv"), sep="\t")
immunebuilder = pd.read_csv(Path("..","data","immunebuilder_rmsd.tsv"), sep="\t")
raptorx = pd.read_csv(Path("..","data","raptorxsingle_rmsd.tsv"), sep="\t")
nanonet = pd.read_csv(Path("..","data","nanonet_rmsd.tsv"), sep="\t")
af3 = pd.read_csv(Path("..","data","alphafold3_rmsd.tsv"), sep="\t")

print(f"af3 shape: {af3.shape} npdbs {af3['pdb'].nunique()}")

# for the other dataset, let's make sure they share the same pdbs
# and that the total number is 40
prev_pdbs = []
for dataset in [af_multimer, af_monomer, immunebuilder, raptorx, nanonet]:
    print(f"dataset {dataset.shape}")
    pdbs = dataset['pdb'].unique()
    list_pdbs = sorted(list(pdbs))
    if prev_pdbs != []:
        assert prev_pdbs == list_pdbs
    assert len(pdbs) == 40
    prev_pdbs = sorted(list(pdbs))

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
methods = [af_multimer_top1, af_monomer_top1, immunebuilder_top1, raptorx, nanonet, af3_top1]
methods_strings = ["af_multimer", "af_monomer", "immunebuilder", "raptorx", "nanonet", "af3"]
for i, method in enumerate(methods):

    FR_mean = round(method['FR_RMS'].mean(), 2)
    FR_std = round(method['FR_RMS'].std(), 2)
    FR_median = round(method['FR_RMS'].median(), 2)
    CDR1_mean = round(method['CDR1_RMS'].mean(), 2)
    CDR1_std = round(method['CDR1_RMS'].std(), 2)
    CDR1_median = round(method['CDR1_RMS'].median(), 2)
    CDR2_mean = round(method['CDR2_RMS'].mean(), 2)
    CDR2_std = round(method['CDR2_RMS'].std(), 2)
    CDR2_median = round(method['CDR2_RMS'].median(), 2)
    CDR3_mean = round(method['CDR3_RMS'].mean(), 2)
    CDR3_std = round(method['CDR3_RMS'].std(), 2)
    CDR3_median = round(method['CDR3_RMS'].median(), 2)
    print(f"method {methods_strings[i]} FR RMSD mean {FR_mean} std {FR_std} median {FR_median}")
    print(f"method {methods_strings[i]} CDR1 RMSD mean {CDR1_mean} std {CDR1_std} median {CDR1_median}")
    print(f"method {methods_strings[i]} CDR2 RMSD mean {CDR2_mean} std {CDR2_std} median {CDR2_median}")
    print(f"method {methods_strings[i]} CDR3 RMSD mean {CDR3_mean} std {CDR3_std} median {CDR3_median}")
    rms_csv.append([FR_mean, FR_std, FR_median, CDR1_mean, CDR1_std, CDR1_median, CDR2_mean, CDR2_std, CDR2_median, CDR3_mean, CDR3_std, CDR3_median])

rms_df = pd.DataFrame(rms_csv, columns=['FR_RMS mean', 'FR_RMS std', 'FR_RMS median', 'CDR1_RMS mean', 'CDR1_RMS std', 'CDR1_RMS median', 'CDR2_RMS mean', 'CDR2_RMS std', 'CDR2_RMS median', 'CDR3_RMS mean', 'CDR3_RMS std', 'CDR3_RMS median'],
                      index=['af_multimer top1', 'af_monomer top1', 'immunebuilder top1', 'raptorx', 'nanonet', 'af3'])

rms_df.to_csv(Path(".", "figures", "nanobody_rmsd_summary.tsv"), sep="\t")

##SUPPLEMENTARY FIGURE 1
##Violin plot showing the distribution of the backbone CDR3 RMSD from the best ranked nanobody prediction for the different methods.

fig, ax = plt.subplots(figsize=(10, 5))

violin = plt.violinplot([immunebuilder_top1['CDR3_RMS'], af_monomer_top1['CDR3_RMS'],
                         af_multimer_top1['CDR3_RMS'], raptorx['CDR3_RMS'],
                         nanonet['CDR3_RMS'], af3_top1['CDR3_RMS']],
                        showmeans=False, showmedians=True, showextrema=False)

violin['cmedians'].set_colors("black")

#we set a different color for each method
colors = plt.cm.tab10.colors
for i in range(len(violin['bodies'])):
    violin['bodies'][i].set_facecolor(colors[i])
    violin['bodies'][i].set_edgecolor('black')
    violin['bodies'][i].set_linewidth(1)

# color the AF3 violin plot differently
violin['bodies'][-1].set_facecolor(plt.cm.tab10.colors[-1])

#we plot the median value
ax.text(1, immunebuilder_top1['CDR3_RMS'].median() + 0.1, str(round(immunebuilder_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(2, af_monomer_top1['CDR3_RMS'].median() + 0.1, str(round(af_monomer_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(3, af_multimer_top1['CDR3_RMS'].median() + 0.1, str(round(af_multimer_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(4, raptorx['CDR3_RMS'].median() + 0.1, str(round(raptorx['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(5, nanonet['CDR3_RMS'].median() + 0.1, str(round(nanonet['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)
ax.text(6, af3_top1['CDR3_RMS'].median() + 0.1, str(round(af3_top1['CDR3_RMS'].median(), 2)), color='black', ha='center', size=12)

# plt.xticks([1, 2, 3, 4, 5, 6], ['IB\n(Top1 models)', 'AF Monomer\n(Top1 models)', 'AF Multimer\n(Top1 models)', 'RXS', 'NN', 'AF3\n(Top1 models)'], size=12)
plt.xticks([1, 2, 3, 4, 5, 6], ['IB', 'AF', 'AF2M', 'RXS', 'NN', 'AF3'], size=12)
plt.ylabel('CDR3 RMSD (Å)', fontsize=15)
plt.title('Top ranked model accuracy', fontsize=15)
plt.tight_layout()
plt.savefig(Path(".", "figures", "SI_figure1.png"), dpi=400)

#FIGURE 1
#Violin plot showing the distribution of the best backbone CDR3 RMSD in the ensembles.

#we get the best available prediction for each pdb in the ensembles
best_ib = immunebuilder.groupby('pdb').min().reset_index()
best_afm = af_monomer.groupby('pdb').min().reset_index()
best_afmm = af_multimer.groupby('pdb').min().reset_index()
best_af3 = af3.groupby('pdb').min().reset_index()
print(f"best_af3 : {best_af3}")
#now we combine the different ensembles
all_rmsd_df = pd.concat([immunebuilder, af_multimer, af_monomer], axis=0)
ib_mu_df = pd.concat([immunebuilder, af_multimer], axis=0)
ib_mo_df = pd.concat([immunebuilder, af_monomer], axis=0)
#and we also get the best available prediction for each pdb in the combined ensembles
all_rmsd_df = all_rmsd_df.groupby('pdb').min().reset_index()
ib_mu_df = ib_mu_df.groupby('pdb').min().reset_index()
ib_mo_df = ib_mo_df.groupby('pdb').min().reset_index()

ib_mo_clustered_centroids = pd.read_csv(Path("..","data","af_monomer_imbuilder_centroids.tsv"), sep="\t")
ib_mu_clustered_centroids = pd.read_csv(Path("..","data","af_multi_imbuilder_centroids.tsv"), sep="\t")
all_clustered_centroids = pd.read_csv(Path("..","data","af_multi_monom_imbuilder_centroids.tsv"), sep="\t")

# for each pdb, we get the models from ib_mo_clustered_centroids and then we get the minimum CDR3 RMSD from the ib_mo_df dataset
dataframes = []
for dataset in [ib_mo_clustered_centroids, ib_mu_clustered_centroids, all_clustered_centroids]:
    dataset_df = pd.DataFrame()
    for pdb in dataset['pdb'].unique():
        models = dataset[dataset['pdb'] == pdb]
        # ib_models are the ones containing the string "imbuilder"
        ib_models = models[models['centroid_file'].str.contains("imbuilder")]
        # alphafold2 monomer models are the ones containing the string "alphafold2_ptm"
        af2_models = models[models['centroid_file'].str.contains("alphafold2_ptm")]
        # alphafold2_multimer models are the ones containing the string "alphafold2_multimer"
        af2m_models = models[models['centroid_file'].str.contains("alphafold2_multimer")]
        # suppress SettingWithCopyWarning
        ib_models = ib_models.copy()
        af2_models = af2_models.copy()
        af2m_models = af2m_models.copy()
        if ib_models.shape[0] > 0:
            # extract the rank of these models: it's the fifth element of centroid_file.split("_")
            ib_models.loc[:, 'rank'] = ib_models['centroid_file'].apply(lambda x: 1 + int(x.split("_")[2].lstrip("rank")))
            # now thanks to the rank we can extract the CDR3_RMS from the corresponding model in immunebuilder df
            ib_models = ib_models.merge(immunebuilder, on=['rank', "pdb"], how='inner')
            # print(ib_models)
        if af2_models.shape[0] > 0:
            af2_models.loc[:, 'rank'] = af2_models['centroid_file'].apply(lambda x: int(x.split("_")[4]))
            af2_models = af2_models.merge(af_monomer, on=['rank', "pdb"], how='inner')
            # print(af2_models)
        if af2m_models.shape[0] > 0:
            af2m_models['rank'] = af2m_models['centroid_file'].apply(lambda x: int(x.split("rank_")[1].split("_")[0]))
            af2m_models = af2m_models.merge(af_multimer, on=['rank', "pdb"], how='inner')
            # print(af2m_models)
        pdb_dataset = pd.concat([ib_models, af2_models, af2m_models], axis=0)
        dataset_df = pd.concat([dataset_df, pdb_dataset], axis=0)
    dataframes.append(dataset_df)
ib_mo_clustered = dataframes[0].groupby('pdb').min().reset_index()
ib_mu_clustered = dataframes[1].groupby('pdb').min().reset_index()
all_clustered = dataframes[2].groupby('pdb').min().reset_index()

##we print mean ± std and median of CDR3 RMSD for each combination
for df, name in zip([best_ib, best_afm, best_afmm, ib_mo_df, ib_mu_df, all_rmsd_df, ib_mo_clustered, ib_mu_clustered, all_clustered, best_af3],
                    ['IB ensemble', 'AF Monomer ensemble', 'AF Multimer ensemble', 'IB + AF Monomer', 'IB + AF Multimer', 'IB + AF Monomer + AF Multimer', 'IB + AF Monomer clustered', 'IB + AF Multimer clustered', 'IB + AF Monomer + AF Multimer clustered', 'AF3 ensemble']):
    print(f"{name}:")
    print(f"Mean ± std: {df['CDR3_RMS'].mean():.2f} ± {df['CDR3_RMS'].std():.2f}")
    print(f"Median: {df['CDR3_RMS'].median():.2f}")

#we plot the violin plot
fig, axs = plt.subplots(1,2, figsize=(15, 5), width_ratios=[10, 7])

violin = axs[0].violinplot([best_ib['CDR3_RMS'], best_afm['CDR3_RMS'], best_afmm['CDR3_RMS'],
                            ib_mo_df['CDR3_RMS'], ib_mu_df['CDR3_RMS'], all_rmsd_df['CDR3_RMS'],
                            ib_mo_clustered['CDR3_RMS'], ib_mu_clustered['CDR3_RMS'], all_clustered['CDR3_RMS'],
                            best_af3['CDR3_RMS']],
    showmeans=False, showmedians=True, showextrema=False)

violin['cmedians'].set_colors("black")

# we plot a dashed line at 3 and 4 angstroms
axs[0].axhline(y=3.5, color='black', linestyle='--')
for i, df in enumerate([best_ib, best_afm, best_afmm, ib_mo_df, ib_mu_df, all_rmsd_df, ib_mo_clustered, ib_mu_clustered, all_clustered, best_af3]):
    # how many structures exceed 3 and 4 angstroms?
    npdbs = df["pdb"].nunique()
    axs[0].text(i+1, 3.7, f"{(len(df[df['CDR3_RMS'] > 3.5])/npdbs*100):.1f}%", color='black', ha='center', size=12)
    # axs[0].text(i+1, 4.2, f"{(len(df[df['CDR3_RMS'] > 4])/npdbs*100):.1f}%", color='black', ha='center', size=12)

#we set a different color for each method
colors = plt.cm.tab10.colors
for i in range(len(violin['bodies'])):
    violin['bodies'][i].set_facecolor(colors[i])
    # violin['bodies'][i].set_edgecolor('black')
    # make the distribution slightly larger
    violin['bodies'][i].set_linewidth(1)

#we plot the median value
axs[0].text(1, best_ib['CDR3_RMS'].median() + 0.1, f"{best_ib['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(2, best_afm['CDR3_RMS'].median() + 0.1, f"{best_afm['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(3, best_afmm['CDR3_RMS'].median() + 0.1, f"{best_afmm['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(4, ib_mo_df['CDR3_RMS'].median() + 0.1, f"{ib_mo_df['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(5, ib_mu_df['CDR3_RMS'].median() + 0.1, f"{ib_mu_df['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(6, all_rmsd_df['CDR3_RMS'].median() + 0.1, f"{all_rmsd_df['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(7, ib_mo_clustered['CDR3_RMS'].median() + 0.1, f"{ib_mo_clustered['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(8, ib_mu_clustered['CDR3_RMS'].median() + 0.1, f"{ib_mu_clustered['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(9, all_clustered['CDR3_RMS'].median() + 0.1, f"{all_clustered['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(10, best_af3['CDR3_RMS'].median() + 0.1, f"{best_af3['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
axs[0].text(-0.1, 0.86, "a)", fontsize=25, va="bottom", transform = axs[0].transAxes)

axs[0].set_xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], ['IB', 'AF', 'AF2M', 'IBMo', 'IBMu', 'IBMM', 'IBMo-Cl', 'IBMu-Cl', 'IBMM-Cl', 'AF3'], size=12)
# set yticks font size
axs[0].tick_params(axis='y', labelsize=12)
axs[0].set_ylabel('CDR3 RMSD (Å)', fontsize=15)

##FIGURE 2B
##Violin plot showing the distribution of the antigen and epitope backbone RMSD

#we read the unbound antigen RMSD data
antigen_rmsd = pd.read_csv(Path("..","data","unbound_antigens_rmsd.tsv"), sep="\t")
# assert the pdbs are 40 and equal to the nanobody pdbs
npdbs = antigen_rmsd['pdb'].nunique()
assert npdbs == 40
assert sorted(list(antigen_rmsd['pdb'].unique())) == sorted(list(prev_pdbs))

# fig, ax = plt.subplots(figsize=(7, 5))
medians = [antigen_rmsd["af_multimer_whole_rmsd"].median(), antigen_rmsd["af_multimer_epitope_rmsd"].median(), antigen_rmsd["af_multimer_pdb_whole_rmsd"].median(), antigen_rmsd["af_multimer_pdb_epitope_rmsd"].median()]

violin = axs[1].violinplot([antigen_rmsd["af_multimer_whole_rmsd"], antigen_rmsd["af_multimer_epitope_rmsd"], antigen_rmsd["af_multimer_pdb_whole_rmsd"], antigen_rmsd["af_multimer_pdb_epitope_rmsd"]], showmeans = False, showmedians = True, showextrema=False)

violin['cmedians'].set_colors("black")

axs[1].set_xticks([1, 2, 3, 4])
axs[1].set_xticklabels(["Antigen RMSD", "Epitope RMSD", "Antigen RMSD", "Epitope RMSD"], size=12)
axs[1].set_ylabel("Backbone RMSD (Å)", fontsize=15)
# set the fontsize of the yticks to 12
axs[1].tick_params(axis='y', labelsize=12)
axs[1].set_ylim(0, 30)

for tick in [1,2,3,4]:
   axs[1].text(tick, medians[tick-1] + 0.2 , round(medians[tick -1], 2),
            horizontalalignment='center',
            style = 'oblique',
            color='black',
            size=12)
axs[1].axhline(y=5, color='black', linestyle='--')
axs[1].text(1, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_whole_rmsd'] > 5])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].text(2, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_epitope_rmsd'] > 5])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].text(3, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_whole_rmsd'] > 5])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].text(4, 5.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_epitope_rmsd'] > 5])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].axhline(y=10, color='black', linestyle='--')
axs[1].text(1, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_whole_rmsd'] > 10])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].text(2, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_epitope_rmsd'] > 10])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].text(3, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_whole_rmsd'] > 10])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
axs[1].text(4, 10.5, f"{len(antigen_rmsd[antigen_rmsd['af_multimer_pdb_epitope_rmsd'] > 10])*100/npdbs:.1f}%", horizontalalignment='center', style = 'oblique', color='black', size=12)
#we print the max value of antigen and epitope RMSD
axs[1].text(1, antigen_rmsd["af_multimer_whole_rmsd"].max(), f"{round(antigen_rmsd['af_multimer_whole_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black', va="bottom", size=12)
axs[1].text(2, antigen_rmsd["af_multimer_epitope_rmsd"].max(), f"{round(antigen_rmsd['af_multimer_epitope_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black', va="bottom", size=12)
axs[1].text(3, antigen_rmsd["af_multimer_pdb_whole_rmsd"].max(), f"{round(antigen_rmsd['af_multimer_pdb_whole_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black', va="bottom", size=12)
axs[1].text(4, antigen_rmsd["af_multimer_pdb_epitope_rmsd"].max(), f"{round(antigen_rmsd['af_multimer_pdb_epitope_rmsd'].max(),2)}Å", horizontalalignment='center', style = 'oblique', color='black', va="bottom", size=12)

axs[1].text(-0.1, 0.86, "b)", fontsize=25, va="bottom", transform = axs[1].transAxes)
#we change the color for the two semibound violins to orange
for pc in violin['bodies']:
    pc.set_facecolor('blue')
    # pc.set_edgecolor('black')
    pc.set_alpha(0.3)
    pc.set_label('AF Multimer')
violin['bodies'][2].set_facecolor('orange')
violin['bodies'][3].set_facecolor('orange')
violin['bodies'][2].set_label('AF Multimer + PDB')
violin['bodies'][3].set_label('AF Multimer + PDB')

#we plot a legend for the colors
leg = axs[1].legend(loc='upper right', labels=['AF Multimer', 'AF Multimer + PDB'], fontsize = 12)
leg.legend_handles[0].set_color('blue')
leg.legend_handles[1].set_color('orange')
leg.legend_handles[0].set_edgecolor('black')
leg.legend_handles[1].set_edgecolor('black')
leg.legend_handles[0].set_alpha(0.5)
leg.legend_handles[1].set_alpha(0.5)

plt.tight_layout()
plt.savefig(Path(".", "figures", "figure1.png"), dpi=400)
plt.close()

# now we do only the violinplot and only for the following entries best_ib, best_afm, best_afmm, ib_mu_df, ib_mu_clustered
fig, ax = plt.subplots(figsize=(8, 5))
violin = plt.violinplot([best_ib['CDR3_RMS'], best_afm['CDR3_RMS'], best_afmm['CDR3_RMS'], ib_mu_df['CDR3_RMS'], ib_mu_clustered['CDR3_RMS']],
                        showmeans=False, showmedians=True, showextrema=False)
violin['cmedians'].set_colors("black")
#we set a different color for each method
colors = plt.cm.tab10.colors
for i in range(len(violin['bodies'])):
    violin['bodies'][i].set_facecolor(colors[i])
    violin['bodies'][i].set_edgecolor('black')
    violin['bodies'][i].set_linewidth(1)
for i, df in enumerate([best_ib, best_afm, best_afmm, ib_mu_df, ib_mu_clustered]):
    # how many structures exceed 3 and 4 angstroms?
    npdbs = df["pdb"].nunique()
    ax.text(i+1, 3.7, f"{(len(df[df['CDR3_RMS'] > 3.5])/npdbs*100):.1f}%", color='black', ha='center', size=12)
ax.axhline(y=3.5, color='black', linestyle='--')
#we plot the median value
ax.text(1, best_ib['CDR3_RMS'].median() + 0.1, f"{best_ib['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
ax.text(2, best_afm['CDR3_RMS'].median() + 0.1, f"{best_afm['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
ax.text(3, best_afmm['CDR3_RMS'].median() + 0.1, f"{best_afmm['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
ax.text(4, ib_mu_df['CDR3_RMS'].median() + 0.1, f"{ib_mu_df['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)
ax.text(5, ib_mu_clustered['CDR3_RMS'].median() + 0.1, f"{ib_mu_clustered['CDR3_RMS'].median():.2f}", color='black', ha='center', size=12)

plt.xticks([1, 2, 3, 4, 5], ['IB', 'AFMo', 'AFMu', 'IBMu', 'IBMu-Cl'], size=12)
# set yticks font size
plt.tick_params(axis='y', labelsize=12)
plt.ylabel('CDR3 RMSD (Å)', fontsize=15)
plt.tight_layout()
plt.savefig(Path(".", "figures", "presentation_figure1.png"), dpi=400)
plt.close()

# a nice analysis of H3_plddt and h3_rmsd for af_monomer
# pick the rank == 1 in the dataframe
af_monomer_top1 = af_monomer[af_monomer['rank'] == 1]
# scatter plot between mean_cdr3_plddt and CDR3_RMS
fig, ax = plt.subplots(figsize=(8, 8))
plt.scatter(af_monomer_top1['mean_cdr3_plddt'], af_monomer_top1['CDR3_RMS'])
# calculate the pearson correlation
corr = af_monomer_top1['mean_cdr3_plddt'].corr(af_monomer_top1['CDR3_RMS'])
print(f"Pearson correlation between mean_cdr3_plddt and CDR3_RMS (af_monomer): {corr:.2f}")
plt.xlabel('mean_cdr3_plddt')
plt.ylabel('CDR3_RMS')
plt.title('Scatter plot between mean_cdr3_plddt and CDR3_RMS')
plt.savefig(Path(".", "figures", "h3-plddt-rmsd-af_monomer.png"), dpi=400)

# now we do the same for af_multimer
af_multimer_top1 = af_multimer[af_multimer['rank'] == 1]
# extract from alphafold2multimer_25_predictions_sr.tsv the dockq of the prediction
af_multimer_dockq = pd.read_csv(Path("..", "data", "alphafold2multimer_25_predictions_sr.tsv"), sep="\t")
af_multimer_top1 = af_multimer_top1.merge(af_multimer_dockq[['pdb', 'rank', 'dockq']], on=['pdb', 'rank'], how='inner')

# scatter plot between mean_cdr3_plddt and CDR3_RMS
fig, ax = plt.subplots(figsize=(8, 8))
plt.scatter(af_multimer_top1['mean_cdr3_plddt'], af_multimer_top1['CDR3_RMS'], c=af_multimer_top1['dockq'], cmap='viridis')
# calculate the pearson correlation
corr = af_multimer_top1['mean_cdr3_plddt'].corr(af_multimer_top1['CDR3_RMS'])
print(f"Pearson correlation between mean_cdr3_plddt and CDR3_RMS (af_multimer): {corr:.2f}")
plt.xlabel('mean_cdr3_plddt')
plt.ylabel('CDR3_RMS')
plt.title('Scatter plot between mean_cdr3_plddt and CDR3_RMS')
plt.savefig(Path(".", "figures", "h3-plddt-rmsd-af_multimer.png"), dpi=400)

alphaflow_pdbs = ["7tprE", "7uiaB", "7unzA", "7wkiB", "7xqvB", "7yz9B", "7zmvG", "8dtnE", "8dtuA"]
alphaflow_centroids = pd.read_csv(Path("..", "data", "centroids_alphaflow_rmsd.tsv"), sep="\t")
aflow_rmsd = pd.read_csv(Path("..", "data", "aflow_predictions_rmsd.tsv"), sep="\t")

all_pdbs = af_monomer['pdb'].unique()
# for all pdbs extract the lowest CDR3_RMS for all the datasets
diff_wrt_centroids = []
diff_wrt_overall = []
for el in all_pdbs:
    min_cdr3_afmo = af_monomer[af_monomer['pdb'] == el]['CDR3_RMS'].min()
    min_cdr3_afmm = af_multimer[af_multimer['pdb'] == el]['CDR3_RMS'].min()
    min_cdr3_ib = immunebuilder[immunebuilder['pdb'] == el]['CDR3_RMS'].min()
    print(f"{el} afmo {min_cdr3_afmo} afmm {min_cdr3_afmm} ib {min_cdr3_ib}")
    if el in alphaflow_pdbs:
        print("This is an AlphaFlow pdb")
        min_af_std = min(min_cdr3_afmo, min_cdr3_afmm)
        # over the alphaflow_centroids, are we doing any better?
        alphaflow_centr_pdb = alphaflow_centroids[alphaflow_centroids['pdb'] == el]
        # print minimum CDR3_RMS
        min_cdr3_aflow_centr = alphaflow_centr_pdb['CDR3_RMS'].min()
        print(f"{el} aflow min CDR3_RMS (20 centroids) {min_cdr3_aflow_centr}")
        # overall
        min_cdr3_aflow_overall = aflow_rmsd[aflow_rmsd['pdb'] == el]['CDR3_RMS'].min()
        print(f"{el} aflow min CDR3_RMS (overall) {min_cdr3_aflow_overall}")
        # do some differences
        diff_cdr3_aflow_centr = min_cdr3_aflow_centr - min_af_std
        diff_cdr3_aflow_overall = min_cdr3_aflow_overall - min_af_std
        print(f"{el} aflow diff CDR3_RMS (20 centroids) {diff_cdr3_aflow_centr}")
        print(f"{el} aflow diff CDR3_RMS (overall) {diff_cdr3_aflow_overall}")
        diff_wrt_centroids.append(diff_cdr3_aflow_centr)
        diff_wrt_overall.append(diff_cdr3_aflow_overall)
import numpy as np
print(f"avg Diff wrt centroids: {np.mean(diff_wrt_centroids)}")
print(f"avg Diff wrt overall: {np.mean(diff_wrt_overall)}")