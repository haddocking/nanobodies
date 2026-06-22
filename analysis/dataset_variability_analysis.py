import pandas as pd
import os
from pathlib import Path
import matplotlib.pyplot as plt

cdr3_df = pd.read_csv(Path("..","data", "benchmark_dataset.tsv"), sep="\t")

cdr3_df["cdr3_length"] = cdr3_df["CDR3"].apply(len)

plt.hist(cdr3_df["cdr3_length"], bins=range(0, 30, 1), edgecolor='black')
plt.xlabel("CDR3 length (residues)", fontsize=12)
plt.ylabel("n of elements in dataset", fontsize=12)
plt.xlim(5, 26)
plt.xticks(range(5, 27, 2))
plt.tight_layout()
plt.savefig(Path("figures", "cdr3_length_distribution.png"), dpi=400)
