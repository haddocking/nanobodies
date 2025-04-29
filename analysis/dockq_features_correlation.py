###Representation of the correlations between max DockQ values and CDR3RMSD (Figure 5A), Epitope RMSD (Supplementary Figure 4) and %Paratope reresented in the restraints (Figure 5B)

import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

##FIGURE 5A
##Scatter plot showing the relation between the minimum CDR3 RMSD in the input models and the maximum DockQ achieved among the top10 models for each ensemble
cdr3rmsd_dockq = pd.read_csv(Path("..", "data/cdr3_rmsd_dockq.tsv"), sep="\t")
cdr3rmsd_dockq_filtered = cdr3rmsd_dockq[(cdr3rmsd_dockq["scenario"] == "real")&(cdr3rmsd_dockq["struct"] == "b")]

#we get the values for each ensemble (only for the real scenario)
cdr3rmsd_dockq_ib = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IB"]
cdr3rmsd_dockq_ib_monom = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IBMo"]
cdr3rmsd_dockq_ib_multim = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IBMu"]
cdr3rmsd_dockq_ib_monom_multim = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IBMoMu"]

fig, axs = plt.subplots(1, 4, figsize=(16.8,5))
gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, 1, 1])
gs.update(wspace=0.5)
axs[0].scatter(cdr3rmsd_dockq_ib["min_cdr3_rmsd"], cdr3rmsd_dockq_ib["max_dockq"], color="purple", label="IB", marker="o", alpha=0.5, edgecolors="black")
axs[0].scatter(cdr3rmsd_dockq_ib_monom["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom["max_dockq"], color="purple", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
axs[0].scatter(cdr3rmsd_dockq_ib_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_multim["max_dockq"], color="purple", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
axs[0].scatter(cdr3rmsd_dockq_ib_monom_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom_multim["max_dockq"], color="purple", label="IBMM", marker="x", alpha=0.5)

axs[0].set_xlabel("Min CDR3 RMSD", size = 18)
axs[0].set_ylabel("Max DockQ in top10 models", size = 18)

#we plot a discontinuous line at dockq = 0.23 (docking success)
axs[0].axhline(y=0.23, color="black", linestyle="--", alpha=0.5)
##we calculate the correlation and plot it
corr_ib = np.corrcoef(cdr3rmsd_dockq_ib["min_cdr3_rmsd"], cdr3rmsd_dockq_ib["max_dockq"])[0,1]
corr_ib_monom = np.corrcoef(cdr3rmsd_dockq_ib_monom["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom["max_dockq"])[0,1]
corr_ib_multim = np.corrcoef(cdr3rmsd_dockq_ib_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_multim["max_dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(cdr3rmsd_dockq_ib_monom_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom_multim["max_dockq"])[0,1]
# print correlation values
print(f"IB correlation: {corr_ib}")
print(f"IBMo correlation: {corr_ib_monom}")
print(f"IBMu correlation: {corr_ib_multim}")
print(f"IBMM correlation: {corr_ib_monom_multim}")

#we plot a at the top right the correlations for each ensemble
axs[0].legend([f"IB $\\rho$ = {corr_ib:.2f}", f"IBMo $\\rho$ = {corr_ib_monom:.2f}", f"IBMu $\\rho$ = {corr_ib_multim:.2f}", f"IBMM $\\rho$ = {corr_ib_monom_multim:.2f}"],
 loc="upper right", fontsize=9)
axs[0].set_title(f"True Interface", size = 18)
plt.savefig(Path("figures", "dockq_analysis", "dockq_cdr3_rmsd_correlation_top10.png"),dpi=400)

# let's re-do the first plot in the loose scenario
cdr3rmsd_dockq_filtered = cdr3rmsd_dockq[(cdr3rmsd_dockq["scenario"] == "loose")&(cdr3rmsd_dockq["struct"] == "b")]
cdr3rmsd_dockq_ib = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IB"]
cdr3rmsd_dockq_ib_monom = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IBMo"]
cdr3rmsd_dockq_ib_multim = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IBMu"]
cdr3rmsd_dockq_ib_monom_multim = cdr3rmsd_dockq_filtered[cdr3rmsd_dockq_filtered["ensemble"] == "IBMoMu"]

axs[1].scatter(cdr3rmsd_dockq_ib["min_cdr3_rmsd"], cdr3rmsd_dockq_ib["max_dockq"], color="purple", label="IB", marker="o", alpha=0.5, edgecolors="black")
axs[1].scatter(cdr3rmsd_dockq_ib_monom["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom["max_dockq"], color="purple", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
axs[1].scatter(cdr3rmsd_dockq_ib_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_multim["max_dockq"], color="purple", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
axs[1].scatter(cdr3rmsd_dockq_ib_monom_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom_multim["max_dockq"], color="purple", label="IBMM", marker="x", alpha=0.5)
axs[1].set_xlabel("Min CDR3 RMSD", size = 18)
#axs[1].set_ylabel("Max DockQ in top10 models", size = 18)

# correlations
corr_ib = np.corrcoef(cdr3rmsd_dockq_ib["min_cdr3_rmsd"], cdr3rmsd_dockq_ib["max_dockq"])[0,1]
corr_ib_monom = np.corrcoef(cdr3rmsd_dockq_ib_monom["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom["max_dockq"])[0,1]
corr_ib_multim = np.corrcoef(cdr3rmsd_dockq_ib_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_multim["max_dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(cdr3rmsd_dockq_ib_monom_multim["min_cdr3_rmsd"], cdr3rmsd_dockq_ib_monom_multim["max_dockq"])[0,1]

axs[1].legend([f"IB $\\rho$ = {corr_ib:.2f}", f"IBMo $\\rho$ = {corr_ib_monom:.2f}", f"IBMu $\\rho$ = {corr_ib_multim:.2f}", f"IBMM $\\rho$ = {corr_ib_monom_multim:.2f}"],
 loc="upper right", fontsize=9)
axs[1].axhline(y=0.23, color="black", linestyle="--", alpha=0.5)
axs[1].set_title(f"Loose Interface", size = 18)

plt.savefig(Path("figures", "dockq_analysis", "dockq_cdr3_rmsd_correlation_top10_loose.png"),dpi=400)

##FIGURE 5B
##Scatter plots showing the relation between the % of paratope represented in the restraints and the maximum DockQ among the top10 models for each ensemble after final refinement

repr_df = pd.read_csv(Path("..", "data/paratope_epitope_repr.tsv"), sep="\t")
# divide the column paratope_represented by 100
repr_df["paratope_represented"] = repr_df["paratope_repr"]/100
loose_dockq = cdr3rmsd_dockq[(cdr3rmsd_dockq["scenario"] == "loose")&(cdr3rmsd_dockq["struct"] == "b")]
two_hit_dockq = cdr3rmsd_dockq[(cdr3rmsd_dockq["scenario"] == "twohit")&(cdr3rmsd_dockq["struct"] == "b")]
# sort all the three datasets by pdb
loose_dockq = loose_dockq.sort_values(by="pdb")
two_hit_dockq = two_hit_dockq.sort_values(by="pdb")
repr_df = repr_df.sort_values(by="pdb")

#Loose Interface
axs[2].scatter(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IB"]["max_dockq"], color="blue", label="IB", marker="o", alpha=0.5, edgecolors="black")
axs[2].scatter(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IBMo"]["max_dockq"], color="blue", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
axs[2].scatter(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IBMu"]["max_dockq"], color="blue", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
axs[2].scatter(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IBMoMu"]["max_dockq"], color="blue", label="IBMM", marker="x", alpha=0.5)
axs[2].set_title(f"Loose Interface", size = 18)
axs[2].set_xlabel("$F_{para}$", size = 18)
# axs[1].set_ylabel("Max DockQ in top10 models", size = 15)
# #we calculate the correlation and plot it at the top left corner
corr_ib = np.corrcoef(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IB"]["max_dockq"])[0,1]
corr_ib_monom = np.corrcoef(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IBMo"]["max_dockq"])[0,1]
corr_ib_multim = np.corrcoef(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IBMu"]["max_dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(repr_df["paratope_represented"], loose_dockq[loose_dockq["ensemble"] == "IBMoMu"]["max_dockq"])[0,1]
# print correlation values
print(f"IB correlation: {corr_ib}")
print(f"IBMo correlation: {corr_ib_monom}")
print(f"IBMu correlation: {corr_ib_multim}")
print(f"IBMM correlation: {corr_ib_monom_multim}")

#we plot a legend for the markers
axs[2].legend([f"IB $\\rho$ = {corr_ib:.2f}", f"IBMo $\\rho$ = {corr_ib_monom:.2f}", f"IBMu $\\rho$ = {corr_ib_multim:.2f}", f"IBMM $\\rho$ = {corr_ib_monom_multim:.2f}"],
 loc="upper left", ncol=1, fontsize=9)	

#we plot a discontinuous line at dockq = 0.23 (docking success)
axs[2].axhline(y=0.23, color="black", linestyle="--", alpha=0.5)

#Two-Hit Interface
print(f'shape repr_df: {repr_df.shape} shape two_hit_dockq: {two_hit_dockq[two_hit_dockq["ensemble"] == "IB"].shape}')
axs[3].scatter(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IB"]["max_dockq"], color="blue", label="IB", marker="o", alpha=0.5, edgecolors="black")
axs[3].scatter(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IBMo"]["max_dockq"], color="blue", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
axs[3].scatter(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IBMu"]["max_dockq"], color="blue", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
axs[3].scatter(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IBMoMu"]["max_dockq"], color="blue", label="IBMM", marker="x", alpha=0.5)
axs[3].set_xlabel("$F_{para}$", size = 18)
axs[3].set_title("Two-Hit Interface", size = 18)

# #we calculate the correlation and plot it at the top left corner

corr_ib = np.corrcoef(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IB"]["max_dockq"])[0,1]
corr_ib_monom = np.corrcoef(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IBMo"]["max_dockq"])[0,1]
corr_ib_multim = np.corrcoef(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IBMu"]["max_dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(repr_df["paratope_represented"], two_hit_dockq[two_hit_dockq["ensemble"] == "IBMoMu"]["max_dockq"])[0,1]

axs[3].legend([f"IB $\\rho$ = {corr_ib:.2f}", f"IBMo $\\rho$ = {corr_ib_monom:.2f}", f"IBMu $\\rho$ = {corr_ib_multim:.2f}", f"IBMM $\\rho$ = {corr_ib_monom_multim:.2f}"],
    loc="upper left", ncol=1, fontsize=9)
#we plot a discontinuous line at dockq = 0.23
axs[3].axhline(y=0.23, color="black", linestyle="--", alpha=0.5)

for i in range(4):
    axs[i].set_ylim(-0.02,1)
    if i != 0:
        # remove yticks
        axs[i].set_yticks([])
#plt.savefig(Path("figures", "dockq_paratope_represented_correlation_top10.png"))
# putting a) and b) text
axs[0].text(-0.12, 1.02, "a)", fontsize=25, va="bottom", transform = axs[0].transAxes)
axs[2].text(-0.08, 1.02, "b)", fontsize=25, va="bottom", transform = axs[2].transAxes)
plt.tight_layout()
pos1 = axs[0].get_position()
pos2 = axs[1].get_position()
pos3 = axs[2].get_position()
pos4 = axs[3].get_position()
axs[1].set_position([pos2.x0 - (pos2.x0 - pos1.x1), pos2.y0, pos2.width, pos2.height])
axs[3].set_position([pos4.x0 - (pos4.x0 - pos3.x1), pos4.y0, pos4.width, pos4.height])

plt.savefig(Path("figures", "figure4.png"), dpi=400)
plt.close()

#SUPPLEMENTARY FIGURE 4
##Scatter plot showing the relationship between maximum DockQ among Top10 nanobody-antigen models and epitope RMSD.
epitope_rmsd_data = pd.read_csv(Path("..", "data/unbound_antigens_rmsd.tsv"), sep="\t")
# order by pdb
epitope_rmsd_data = epitope_rmsd_data.sort_values(by="pdb")
# drop the 7qneG
epitope_rmsd_data = epitope_rmsd_data[epitope_rmsd_data["pdb"] != "7qneG"]
# sort all the cdr3rmsd_dockq_ib and similar by pdb
cdr3rmsd_dockq_filtered = cdr3rmsd_dockq[(cdr3rmsd_dockq["scenario"] == "real")&(cdr3rmsd_dockq["struct"] == "u")]
cdr3rmsd_dockq_ib = cdr3rmsd_dockq_filtered[(cdr3rmsd_dockq_filtered["ensemble"] == "IB")&(cdr3rmsd_dockq_filtered["pdb"] != "7qneG")].sort_values(by="pdb")
cdr3rmsd_dockq_ib_monom = cdr3rmsd_dockq_filtered[(cdr3rmsd_dockq_filtered["ensemble"] == "IBMo")&(cdr3rmsd_dockq_filtered["pdb"] != "7qneG")].sort_values(by="pdb")
cdr3rmsd_dockq_ib_multim = cdr3rmsd_dockq_filtered[(cdr3rmsd_dockq_filtered["ensemble"] == "IBMu")&(cdr3rmsd_dockq_filtered["pdb"] != "7qneG")].sort_values(by="pdb")
cdr3rmsd_dockq_ib_monom_multim = cdr3rmsd_dockq_filtered[(cdr3rmsd_dockq_filtered["ensemble"] == "IBMoMu")&(cdr3rmsd_dockq_filtered["pdb"] != "7qneG")].sort_values(by="pdb")
# assert the same number of rows
assert epitope_rmsd_data.shape[0] == cdr3rmsd_dockq_ib.shape[0] == 39, f"Expected 40 rows, found {epitope_rmsd_data.shape[0]} {cdr3rmsd_dockq_ib.shape[0]} {cdr3rmsd_dockq_ib_monom.shape[0]} {cdr3rmsd_dockq_ib_multim.shape[0]} {cdr3rmsd_dockq_ib_monom_multim.shape[0]}"
#we will plot the data
fig, ax = plt.subplots(1, 1, figsize=(5,5))

ax.scatter(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib["max_dockq"], color="green", label="IB", marker="o", alpha=0.5, edgecolors="black")
ax.scatter(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib_monom["max_dockq"], color="green", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
ax.scatter(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib_multim["max_dockq"], color="green", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
ax.scatter(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib_monom_multim["max_dockq"], color="green", label="IBMM", marker="x", alpha=0.5)
ax.set_xlabel("Unbound Epitope RMSD", size = 15)
ax.set_ylabel("Max DockQ in top10 models", size = 15)
#we plot a discontinuous line at dockq = 0.23
ax.axhline(y=0.23, color="black", linestyle="--", alpha=0.5)
#we calculate the correlation and plot it
corr_ib = np.corrcoef(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib["max_dockq"])[0,1]
corr_ib_monom = np.corrcoef(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib_monom["max_dockq"])[0,1]
corr_ib_multim = np.corrcoef(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib_multim["max_dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(epitope_rmsd_data["af_multimer_pdb_epitope_rmsd"], cdr3rmsd_dockq_ib_monom_multim["max_dockq"])[0,1]

#we plot a at the top right the correlations for each ensemble
ax.legend(["IB " + r"$\rho$ = " + str(round(corr_ib,4)), "IBMo " + r"$\rho$ = " + str(round(corr_ib_monom,4)), "IBMu " + r"$\rho$ = " + str(round(corr_ib_multim,4)), "IBMM " + r"$\rho$ = " + str(round(corr_ib_monom_multim,4))],
    loc="upper right", fontsize=10)

plt.savefig(Path( "figures", "SI_figure4.png"), dpi=400)
plt.close()

