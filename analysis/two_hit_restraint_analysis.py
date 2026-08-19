import pandas as pd
from pathlib import Path
import os

def extract_avg_t10_dockq(pdb):
    capri_path = pd.read_csv(Path("..", "benchmark_pdb_files", pdb, f"{pdb}_capri.tsv"), sep="\t")
    # filter: keep only stage = emref data and scenario = two-hit
    capri_df = capri_path[(capri_path["stage"] == "emref") & (capri_path["scenario"] == "twohit") & (capri_path["ensemble"] == "IBMu")]
    
    dockq_bound = capri_df[(capri_path["struct"] == "b") & (capri_df["caprieval_rank"] < 11)]
    dockq_unbound = capri_df[(capri_df["struct"] == "u") & (capri_df["caprieval_rank"] < 11)]
    print(dockq_bound[["struct", "caprieval_rank", "dockq"]])
    avg_dockq_bound = dockq_bound["dockq"].mean()
    avg_dockq_unbound = dockq_unbound["dockq"].mean()
    max_dockq_bound = dockq_bound["dockq"].max()
    max_dockq_unbound = dockq_unbound["dockq"].max()
    print(f"Average dockq for bound: {avg_dockq_bound}")
    print(f"Average dockq for unbound: {avg_dockq_unbound}")
    print(f"Max dockq for bound: {max_dockq_bound}")
    print(f"Max dockq for unbound: {max_dockq_unbound}")
    return avg_dockq_bound, avg_dockq_unbound, max_dockq_bound, max_dockq_unbound

def read_ambig(ambig_path):
    first_ass = False
    epitope_residues = []
    with open(ambig_path, "r") as f:
        print(f"Reading ambig file: {ambig_path}")
        for line in f:
            if line.startswith("assign") and first_ass is False:
                first_ass = True
            elif line.startswith("assign") and first_ass is True:
                break
            elif "segid B)" in line:
                epi_res = int(line.split()[1])
                epitope_residues.append(epi_res)
    print(f"Extracted epitope residues: {epitope_residues}")
    return epitope_residues


def extract_restrain_perc(pdb):
    two_hit_ambig = Path("..", "benchmark_pdb_files", pdb, f"{pdb}_twohit_ambig.tbl")
    two_hit_residues = read_ambig(two_hit_ambig)
    real_ambig = Path("..", "benchmark_pdb_files", pdb, f"{pdb}_real_ambig.tbl")
    real_residues = read_ambig(real_ambig)
    # quantify how much of the real epitope residues are restrained in the two-hit scenario
    restrained_residues = set(real_residues).intersection(set(two_hit_residues))
    print(restrained_residues)
    perc = len(restrained_residues) / len(real_residues) * 100
    print(f"Percentage of real epitope residues restrained in two-hit scenario: {perc}")
    return perc

#extract_avg_t10_dockq("8en2D")

#extract_restrain_perc("8en2D")

analysis_filename = "../data/two_hit_restraint_analysis.csv"
if Path(analysis_filename).exists():
    print(f"{analysis_filename} already exists. Skipping analysis.")
    df = pd.read_csv(analysis_filename)
else:
    print(f"{analysis_filename} does not exist. Running analysis.")
    data = []
    for el in os.listdir(Path("..", "benchmark_pdb_files")):
        if os.path.isdir(Path("..", "benchmark_pdb_files", el)):
            dockq_b, dockq_u, max_b, max_u = extract_avg_t10_dockq(el)
            perc = extract_restrain_perc(el)
            data.append([el, dockq_b, dockq_u, max_b, max_u, perc])

    df = pd.DataFrame(data, columns=["pdb", "avg_dockq_bound", "avg_dockq_unbound", "max_dockq_bound", "max_dockq_unbound", "perc_restrained"])
    df.to_csv("../data/two_hit_restraint_analysis.csv", index=False)

# make a 2x2 plot with dock values against percentage of restrained residues
import matplotlib.pyplot as plt
ax, axs = plt.subplots(2, 2, figsize=(10, 10))
axs[0, 0].scatter(df["perc_restrained"], df["avg_dockq_bound"])
axs[0, 0].set_xlabel("Percentage of true epitope residues in two-hit restraints")
axs[0, 0].set_ylabel("Average DockQ (T10 bound)")
axs[0, 1].scatter(df["perc_restrained"], df["avg_dockq_unbound"])
axs[0, 1].set_xlabel("Percentage of true epitope residues in two-hit restraints")
axs[0, 1].set_ylabel("Average DockQ (T10 unbound)")
axs[1, 0].scatter(df["perc_restrained"], df["max_dockq_bound"])
axs[1, 0].set_xlabel("Percentage of true epitope residues in two-hit restraints")
axs[1, 0].set_ylabel("Max DockQ (T10 bound)")
axs[1, 1].scatter(df["perc_restrained"], df["max_dockq_unbound"])
axs[1, 1].set_xlabel("Percentage of true epitope residues in two-hit restraints")
axs[1, 1].set_ylabel("Max DockQ (T10 unbound)")
# color the point
plt.tight_layout()
plt.savefig("figures/SI_figure4.png")