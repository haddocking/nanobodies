import pandas as pd
import os
from pathlib import Path
import matplotlib.pyplot as plt

#CDR3 length distribution analysis

cdr3_df = pd.read_csv(Path("..","data", "benchmark_dataset.tsv"), sep="\t")

cdr3_df["cdr3_length"] = cdr3_df["CDR3"].apply(len)

#we print the mean, max, min and median of the CDR3 length distribution
print("Mean CDR3 length: ", cdr3_df["cdr3_length"].mean(), "residues")
print("Max CDR3 length: ", cdr3_df["cdr3_length"].max(), "residues")
print("Min CDR3 length: ", cdr3_df["cdr3_length"].min(), "residues")
print("Median CDR3 length: ", cdr3_df["cdr3_length"].median(), "residues")

plt.hist(cdr3_df["cdr3_length"], bins=range(0, 30, 1), edgecolor='black')
plt.xlabel("CDR3 length (residues)", fontsize=12)
plt.ylabel("n of elements in dataset", fontsize=12)
plt.xlim(5, 26)
plt.xticks(range(5, 27, 2))
plt.tight_layout()
plt.savefig(Path("figures", "cdr3_length_distribution.png"), dpi=400)
plt.close()

#Antigen size distribution analysis

antigen_size_df = pd.read_csv(Path("..","data", "antigen_size.tsv"), sep="\t")

print("Mean antigen size: ", antigen_size_df["size (kDa)"].mean(), "kDa.")
print("Max antigen size: ", antigen_size_df["size (kDa)"].max(), "kDa.")
print("Min antigen size: ", antigen_size_df["size (kDa)"].min(), "kDa.")
print("Median antigen size: ", antigen_size_df["size (kDa)"].median(), "kDa.")

plt.hist(antigen_size_df["size (kDa)"], bins=range(0, 61, 5), edgecolor='black', color='orange')
plt.xlabel("Antigen size (kDa)", fontsize=12)
plt.ylabel("n of elements in dataset", fontsize=12)
plt.xlim(0, 60)
plt.xticks(range(0, 61, 10))
plt.tight_layout()
plt.savefig(Path("figures", "antigen_size_distribution.png"), dpi=400)
plt.close()