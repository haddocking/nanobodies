import pandas as pd
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#CDR3 length distribution analysis

cdr3_df = pd.read_csv(Path("..","data", "benchmark_dataset.tsv"), sep="\t")

cdr3_df["cdr3_length"] = cdr3_df["CDR3"].apply(len)

#we print the mean, max, min and median of the CDR3 length distribution
print("Mean CDR3 length: ", cdr3_df["cdr3_length"].mean(), "residues")
print("Max CDR3 length: ", cdr3_df["cdr3_length"].max(), "residues")
print("Min CDR3 length: ", cdr3_df["cdr3_length"].min(), "residues")
print("Median CDR3 length: ", cdr3_df["cdr3_length"].median(), "residues")

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# plt.hist(cdr3_df["cdr3_length"], bins=range(0, 30, 1), edgecolor='black')
# plt.xlabel("CDR3 length (residues)", fontsize=12)
# plt.ylabel("n of elements in dataset", fontsize=12)
# plt.xlim(5, 26)
# plt.xticks(range(5, 27, 2))
# plt.tight_layout()
# plt.savefig(Path("figures", "cdr3_length_distribution.png"), dpi=400)
# plt.close()

ax[0].hist(cdr3_df["cdr3_length"], bins=range(0, 30, 1), edgecolor='black', color='blue')
ax[0].set_xlabel("CDR3 length (residues)", fontsize=12)
ax[0].set_ylabel("n of elements in dataset", fontsize=12)
ax[0].set_xlim(5, 26)
ax[0].set_xticks(range(5, 27, 2))


#Antigen size distribution analysis

antigen_size_df = pd.read_csv(Path("..","data", "antigen_size.tsv"), sep="\t")

print("Mean antigen size: ", antigen_size_df["size (kDa)"].mean(), "kDa.")
print("Max antigen size: ", antigen_size_df["size (kDa)"].max(), "kDa.")
print("Min antigen size: ", antigen_size_df["size (kDa)"].min(), "kDa.")
print("Median antigen size: ", antigen_size_df["size (kDa)"].median(), "kDa.")

ax[1].hist(antigen_size_df["size (kDa)"], bins=range(0, 61, 5), edgecolor='black', color='orange')
ax[1].set_xlabel("Antigen size (kDa)", fontsize=12)
ax[1].set_ylabel("n of elements in dataset", fontsize=12)
ax[1].set_xlim(0, 60)
ax[1].set_xticks(range(0, 61, 10))
#plt.tight_layout()
#plt.savefig(Path("figures", "antigen_size_distribution.png"), dpi=400)
#plt.close()

ax[0].text(-0.20, 1.02, "a)", fontsize=20, va="bottom", transform = ax[0].transAxes)
ax[1].text(-0.20, 1.02, "b)", fontsize=20, va="bottom", transform = ax[1].transAxes)
plt.subplots_adjust(wspace=0.4)


plt.savefig(Path("figures", "SI_figure1_trial.png"), dpi=400)