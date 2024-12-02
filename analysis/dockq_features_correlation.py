###Representation of the correlations between max DockQ values and CDR3RMSD (Figure 5A), Epitope RMSD (Supplementary Figure 4) and %Paratope reresented in the restraints (Figure 5B)

import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

##FIGURE 5A
##Scatter plot showing the relation between the minimum CDR3 RMSD in the input models and the maximum DockQ achieved among the top10 models for each ensemble

cdr3rmsd_dockq = pd.read_csv(Path("..", "data/cdr3_rmsd_dockq.tsv"), sep="\t")

#we get the values for each ensemble
cdr3rmsd_dockq_ib = cdr3rmsd_dockq[cdr3rmsd_dockq["ensemble"] == "ib"]
cdr3rmsd_dockq_ib_monom = cdr3rmsd_dockq[cdr3rmsd_dockq["ensemble"] == "ib_monom"]
cdr3rmsd_dockq_ib_multim = cdr3rmsd_dockq[cdr3rmsd_dockq["ensemble"] == "ib_multim"]
cdr3rmsd_dockq_ib_monom_multim = cdr3rmsd_dockq[cdr3rmsd_dockq["ensemble"] == "ib_monom_multim"]

#we plot it only for true scenario and emref_08
fig, ax = plt.subplots(1, 1, figsize=(5,5))

ax.scatter(cdr3rmsd_dockq_ib["cdr3_rmsd"], cdr3rmsd_dockq_ib["dockq"], color="purple", label="IB", marker="o", alpha=0.5, edgecolors="black")
ax.scatter(cdr3rmsd_dockq_ib_monom["cdr3_rmsd"], cdr3rmsd_dockq_ib_monom["dockq"], color="purple", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
ax.scatter(cdr3rmsd_dockq_ib_multim["cdr3_rmsd"], cdr3rmsd_dockq_ib_multim["dockq"], color="purple", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
ax.scatter(cdr3rmsd_dockq_ib_monom_multim["cdr3_rmsd"], cdr3rmsd_dockq_ib_monom_multim["dockq"], color="purple", label="IBMM", marker="x", alpha=0.5, edgecolors="black")

ax.set_xlabel("Min CDR3 RMSD", size = 15)
ax.set_ylabel("Max DockQ in top10 models", size = 15)

#we plot a discontinuous line at dockq = 0.23 (docking success)
ax.axhline(y=0.23, color="black", linestyle="--", alpha=0.5)
##we calculate the correlation and plot it
corr_ib = np.corrcoef(cdr3rmsd_dockq_ib["cdr3_rmsd"], cdr3rmsd_dockq_ib["dockq"])[0,1]
corr_ib_monom = np.corrcoef(cdr3rmsd_dockq_ib_monom["cdr3_rmsd"], cdr3rmsd_dockq_ib_monom["dockq"])[0,1]
corr_ib_multim = np.corrcoef(cdr3rmsd_dockq_ib_multim["cdr3_rmsd"], cdr3rmsd_dockq_ib_multim["dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(cdr3rmsd_dockq_ib_monom_multim["cdr3_rmsd"], cdr3rmsd_dockq_ib_monom_multim["dockq"])[0,1]

#we plot a at the top right the correlations for each ensemble
ax.legend(["IB " + r"$\rho$ = " + str(round(corr_ib,4)), "IBMo " + r"$\rho$ = " + str(round(corr_ib_monom,4)), "IBMu " + r"$\rho$ = " + str(round(corr_ib_multim,4)), "IBMM " + r"$\rho$ = " + str(round(corr_ib_monom_multim,4))],
 loc="upper right", fontsize=8)

plt.savefig(Path("figures", "dockq_cdr3_rmsd_correlation_top10.png"))



##SUPPLEMENTARY FIGURE 4
##Scatter plot showing the relation between maximum DockQ among Top10 nanobody-antigen models and epitope RMSD.

epitope_rmsd_dockq = pd.read_csv(Path("..", "data/epitope_rmsd_dockq.tsv"), sep="\t")

#we will plot the data
fig, ax = plt.subplots(1, 1, figsize=(5,5))

ax.scatter(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib"]["dockq"], color="green", label="IB", marker="o", alpha=0.5, edgecolors="black")
ax.scatter(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom"]["dockq"], color="green", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
ax.scatter(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_multim"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_multim"]["dockq"], color="green", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
ax.scatter(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom_multim"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom_multim"]["dockq"], color="green", label="IBMM", marker="x", alpha=0.5, edgecolors="black")
ax.set_xlabel("Unbound Epitope RMSD", size = 15)
ax.set_ylabel("Max DockQ in top10 models", size = 15)
#we plot a discontinuous line at dockq = 0.23
ax.axhline(y=0.23, color="black", linestyle="--", alpha=0.5)
#we calculate the correlation and plot it
corr_ib = np.corrcoef(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib"]["dockq"])[0,1]
corr_ib_monom = np.corrcoef(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom"]["dockq"])[0,1]
corr_ib_multim = np.corrcoef(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_multim"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_multim"]["dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom_multim"]["epitope_rmsd"], epitope_rmsd_dockq[epitope_rmsd_dockq["ensemble"] == "ib_monom_multim"]["dockq"])[0,1]

#we plot a at the top right the correlations for each ensemble
ax.legend(["IB " + r"$\rho$ = " + str(round(corr_ib,4)), "IBMo " + r"$\rho$ = " + str(round(corr_ib_monom,4)), "IBMu " + r"$\rho$ = " + str(round(corr_ib_multim,4)), "IBMM " + r"$\rho$ = " + str(round(corr_ib_monom_multim,4))],
    loc="upper right", fontsize=10)


plt.savefig(Path( "figures", "dockq_epitope_rmsd_correlation_top10.png"))



##FIGURE 5B
##Scatter plots showing the relation between the % of paratope represented in the restraints and the maximum DockQ among the top10 models for each ensemble after final refinement

represented_paratope_dockq = pd.read_csv(Path("..", "data/paratope_represented_dockq.tsv"), sep="\t")
loose_represented_paratope_dockq = represented_paratope_dockq[represented_paratope_dockq["scenario"] == "loose"]
two_hit_represented_paratope_dockq = represented_paratope_dockq[represented_paratope_dockq["scenario"] == "two-hit"]

#we will plot the data
fig, axs = plt.subplots(1,2, figsize=(10,5))

#Loose Interface
axs[0].scatter(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib"]["dockq"], color="blue", label="IB", marker="o", alpha=0.5, edgecolors="black")
axs[0].scatter(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom"]["dockq"], color="blue", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
axs[0].scatter(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_multim"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_multim"]["dockq"], color="blue", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
axs[0].scatter(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["dockq"], color="blue", label="IBMM", marker="x", alpha=0.5, edgecolors="black")
axs[0].set_title(f"Loose Interface", size = 15)
axs[0].set_xlabel("% Paratope represented", size = 15)
axs[0].set_ylabel("Max DockQ in top10 models", size = 15)
# #we calculate the correlation and plot it at the top left corner
corr_ib = np.corrcoef(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib"]["dockq"])[0,1]
corr_ib_monom = np.corrcoef(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom"]["dockq"])[0,1]
corr_ib_multim = np.corrcoef(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_multim"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_multim"]["dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["paratope_represented"], loose_represented_paratope_dockq[loose_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["dockq"])[0,1]

#we plot a legend for the markers
axs[0].legend(["IB "+ r'$\rho$ = ' + f"{np.round(corr_ib,4)}", "IBMo "+r'$\rho$ = ' + f"{np.round(corr_ib_monom,4)}","IBMu "+ r'$\rho$ = ' + f"{np.round(corr_ib_multim,4)}","IBMM "+ r'$\rho$ = ' + f"{np.round(corr_ib_monom_multim,4)}"],
 loc="upper left", ncol=1, fontsize=8)	

#we plot a discontinuous line at dockq = 0.23 (docking success)
axs[0].axhline(y=0.23, color="black", linestyle="--", alpha=0.5)

#Two-Hit Interface
axs[1].scatter(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib"]["dockq"], color="blue", label="IB", marker="o", alpha=0.5, edgecolors="black")
axs[1].scatter(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom"]["dockq"], color="blue", label="IBMo", marker="^", alpha=0.5, edgecolors="black")
axs[1].scatter(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_multim"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_multim"]["dockq"], color="blue", label="IBMu", marker="s", alpha=0.5, edgecolors="black")
axs[1].scatter(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["dockq"], color="blue", label="IBMM", marker="x", alpha=0.5, edgecolors="black")
axs[1].set_xlabel("% Paratope represented", size = 15)
axs[1].set_title("Two-Hit Interface", size = 15)

# #we calculate the correlation and plot it at the top left corner
corr_ib = np.corrcoef(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib"]["dockq"])[0,1]
corr_ib_monom = np.corrcoef(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom"]["dockq"])[0,1]
corr_ib_multim = np.corrcoef(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_multim"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_multim"]["dockq"])[0,1]
corr_ib_monom_multim = np.corrcoef(two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["paratope_represented"], two_hit_represented_paratope_dockq[two_hit_represented_paratope_dockq["ensemble"] == "ib_monom_multim"]["dockq"])[0,1]

axs[1].legend(["IB "+r'$\rho$ = ' + f"{np.round(corr_ib,4)}","IBMo "+ r'$\rho$ = ' + f"{np.round(corr_ib_monom,4)}","IBMu "+ r'$\rho$ = ' + f"{np.round(corr_ib_multim,4)}","IBMM "+ r'$\rho$ = ' + f"{np.round(corr_ib_monom_multim,4)}"],
    loc="upper left", ncol=1, fontsize=8)
#we plot a discontinuous line at dockq = 0.23
axs[1].axhline(y=0.23, color="black", linestyle="--", alpha=0.5)

plt.savefig(Path("figures", "dockq_paratope_represented_correlation_top10.png"))
