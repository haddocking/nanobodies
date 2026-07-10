#Retrieves the epitope from every antigen and computes the TM-score between every pair of epitopes. 

import itertools
from itertools import combinations
import os
from pathlib import Path
import MDAnalysis as mda
from scipy.spatial.distance import cdist
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from tmtools import tm_align

nanobodies = pd.read_csv(Path("..", "data", "benchmark_dataset.tsv"), sep="\t")

# ##Calculate epitopes

# def keep_atoms(inp_pdb_f):
#     """
#     Keep atoms.

#     Parameters
#     ----------
#     inp_pdb_f : Path
#         Path to PDB file.

#     Returns
#     -------
#     out_pdb_fname : Path
#         Path to PDB file.
#     """
#     out_pdb_fname = Path(f"{inp_pdb_f.stem}-atoms.pdb")
#     with open(inp_pdb_f, "r") as pdb_fh:
#         with open(out_pdb_fname, "w") as f:
#             for line in pdb_fh:
#                 if line.startswith("ATOM"):
#                     f.write(line)
#     return out_pdb_fname


# def load_contacts_mda(pdb_f, cutoff=4.5): #distance cutoff of 4.5 Angstroms
#         """
#         Load residue-based contacts.
#         """
#         # get chains
#         mdu = mda.Universe(pdb_f)
#         con_list, con_num_list = [], []
#         unique_chains = [str(el)[-2] for el in mdu.segments]
#         unique_chain_combs = list(combinations(unique_chains, 2))
#         # calculating contacts
#         for pair in unique_chain_combs:
#                 chain_n = mdu.select_atoms(f"chainID {pair[0]} and not name H*")
#                 chain_s = mdu.select_atoms(f"chainID {pair[1]} and not name H*")
#                 dist = cdist(chain_n.positions, chain_s.positions)
#                 npw = np.where(dist < cutoff)
#                 del dist
#                 for k in range(npw[0].shape[0]):
#                         con = (pair[0], f"{chain_n.atoms[npw[0][k]].resid}{chain_n.atoms[npw[0][k]].icode}", pair[1], f"{chain_s.atoms[npw[1][k]].resid}{chain_s.atoms[npw[1][k]].icode}")
#                         con_num = (pair[0], chain_n.atoms[npw[0][k]].resid, pair[1], chain_s.atoms[npw[1][k]].resid)
#                         #print(f"con {con}")
#                         con_list.append(con)
#                         con_num_list.append(con_num)
#         return [set(con_list), set(con_num_list)]

# per_pdb_epitope = []

# for n in range(nanobodies.shape[0]):
#     pdb, chain = nanobodies.iloc[n]["PDB"], nanobodies.iloc[n]["Nb Chain"]

#     #First we get the complex pdb file
#     complex_pdb = Path("..", "benchmark_pdb_files", f"{pdb}{chain}", f"{pdb}{chain}_ref.pdb")
#     complex_atoms = keep_atoms(complex_pdb)
#     # #We calculate the real interface scenario
#     mda_contacts = load_contacts_mda(complex_atoms)
#     os.remove(complex_atoms)
#     contact = mda_contacts[0]
#     contact_nums = mda_contacts[1]
#     #print(f"{pdb}{chain} contacts: {contact_nums}")
#     inter_resnums = []
#     for con in contact_nums:
#         if con[0] == "A" and con[2] == "B":
#             inter_resnums.append(con[3])
#         elif con[2] == "A" and con[0] == "B":
#             inter_resnums.append(con[1])
    
#     unique_inter_resnums = list(set(inter_resnums))
#     #print(pdb, chain, unique_inter_resnums)
    
#     per_pdb_epitope.append([pdb, chain, ",".join([str(res) for res in sorted(unique_inter_resnums)])])

# per_pdb_epitope_df = pd.DataFrame(per_pdb_epitope, columns=["pdb", "chain", "epitope"])
# #print(per_pdb_epitope_df.head())
# per_pdb_epitope_df.to_csv(Path("..", "data", "benchmark_dataset_epitopes.tsv"), sep="\t", index=False)

#We extract the epitopes from the antigen PDB file and compute the TM-score between every pair of epitopes

def parse_epitope(epitope_str, default_chain="B"):
    """
    Turn "102,103,A150,A151" into {"A": {102, 103, 150, 151}}
    (assuming default_chain='A' for the bare numbers).
    Returns dict: chain_id -> set of residue numbers.
    """
    chain_residues = defaultdict(set)
    for token in str(epitope_str).split(","):
        token = token.strip()
        if not token:
            continue
        if token[0].isalpha():
            chain_id, resnum = token[0], int(token[1:])
        else:
            chain_id, resnum = default_chain, int(token)
        chain_residues[chain_id].add(resnum)
    return chain_residues
 
 
def pdb_path_for(pdb_id, pdb_chain):
    return os.path.join("../benchmark_pdb_files", pdb_id+pdb_chain, f"{pdb_id+pdb_chain}_l_b.pdb")

def get_epitope_coords_seq(structure, chain_residues):
    """
    Walk the parsed structure and pull out CA coordinates + 1-letter
    sequence only for residues listed in chain_residues, directly from the
    Biopython Structure object. Nothing is written to disk.
    """
    coords, seq = [], []
    model = structure[0]  # first model
    for chain in model:
        if chain.id not in chain_residues:
            continue
        wanted = chain_residues[chain.id]
        for residue in chain.get_residues():
            # residue.id[0] == " " excludes waters/heteroatoms
            if residue.id[0] == " " and residue.id[1] in wanted:
                if "CA" in residue.child_dict:
                    coords.append(residue.child_dict["CA"].coord)
                    seq.append(protein_letters_3to1.get(residue.resname, "X"))
 
    if not coords:
        raise ValueError(
            f"No epitope residues found (requested chains/residues: "
            f"{dict(chain_residues)}). Check chain IDs / residue numbering."
        )
 
    return np.vstack(coords), "".join(seq)
 
 
def load_all_epitope_data(df):
    """Read the CSV rows, parse each PDB once, extract epitope data. Returns
    {pdb_id+pdb_chain: (coords, seq)}."""
    parser = PDBParser(QUIET=True)
    data = {}
 
    for _, row in df.iterrows():
        pdb_id = str(row["pdb"]).strip()
        pdb_chain = str(row["chain"]).strip()
        chain_residues = parse_epitope(row["epitope"])
 
        path = pdb_path_for(pdb_id, pdb_chain)
        structure = parser.get_structure(pdb_id, path)
 
        coords, seq = get_epitope_coords_seq(structure, chain_residues)
        if len(seq) < 3:
            print(f"WARNING: {pdb_id}{pdb_chain} epitope has only {len(seq)} residues "
                  f"with CA atoms; TM-align needs at least 3.")
 
        data[pdb_id+pdb_chain] = (coords, seq)
        print(f"{pdb_id}{pdb_chain}: {len(seq)} epitope residues loaded")
 
    return data

def compute_tmscore_matrix(epitope_data):
    """
    Pairwise TM-score matrix. TM-align returns two scores per pair (each
    normalized by one structure's length); we average them for a symmetric
    matrix, the usual convention when comparing many structures of
    different lengths.
    """
    names = list(epitope_data.keys())
    n = len(names)
    tm_matrix = np.eye(n)
    #rmsd_matrix = np.zeros((n, n))
 
    for i, j in itertools.combinations(range(n), 2):
        coords_i, seq_i = epitope_data[names[i]]
        coords_j, seq_j = epitope_data[names[j]]
 
        result = tm_align(coords_i, coords_j, seq_i, seq_j)
 
        tm_score = (result.tm_norm_chain1 + result.tm_norm_chain2) / 2
        tm_matrix[i, j] = tm_matrix[j, i] = tm_score
        #rmsd_matrix[i, j] = rmsd_matrix[j, i] = result.rmsd
 
        print(f"{names[i]:15s} vs {names[j]:15s}  TM-score={tm_score:.3f}  ")
              #f"RMSD={result.rmsd:.2f} A")
 
    tm_df = pd.DataFrame(tm_matrix, index=names, columns=names)
    #rmsd_df = pd.DataFrame(rmsd_matrix, index=names, columns=names)
    return tm_df#, rmsd_df

def plot_heatmap(tm_df, out_path="tmscore_heatmap.png"):
    """
    Plot the pairwise TM-score matrix showing each pair only once (upper
    triangle, diagonal skipped since self-similarity is trivially 1.0), and
    explicitly highlight the minimum and maximum TM-score pairs.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
 
    # Mask the lower triangle + diagonal -> only the upper triangle is drawn
    mask = np.tril(np.ones_like(tm_df, dtype=bool))
 
    # Find min/max among the *visible* (upper-triangle) values only
    upper_vals = tm_df.where(~mask)
    min_val = upper_vals.min().min()
    max_val = upper_vals.max().max()
    min_pair = upper_vals.stack().idxmin()  # (row_name, col_name)
    max_pair = upper_vals.stack().idxmax()
 
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        tm_df,
        mask=mask,
        cmap="viridis",
        vmin=0,
        vmax=1,
        square=True,
        cbar_kws={"label": "TM-score"},
        ax=ax,
    )
 
    # Draw a colored box around the min and max cells so they're easy to spot
    for pair, color in [(min_pair, "red"), (max_pair, "lime")]:
        row_name, col_name = pair
        r = tm_df.index.get_loc(row_name)
        c = tm_df.columns.get_loc(col_name)
        ax.add_patch(
            plt.Rectangle((c, r), 1, 1, fill=False, edgecolor=color, lw=3)
        )
 
    # subtitle = (
    #     f"min = {min_val:.3f} ({min_pair[0]} vs {min_pair[1]}, red)   "
    #     f"max = {max_val:.3f} ({max_pair[0]} vs {max_pair[1]}, green)"
    # )
    #ax.set_title(f"Pairwise epitope TM-score\n{subtitle}", fontsize=10)

    ax.text(-0.15, 1.02, "c)", fontsize=20, va="bottom", transform = ax.transAxes)
 
    #plt.tight_layout()
    plt.savefig(out_path, dpi=400)
    print(f"Saved heatmap to {out_path}")
    print(f"Min TM-score: {min_val:.3f} ({min_pair[0]} vs {min_pair[1]})")
    print(f"Max TM-score: {max_val:.3f} ({max_pair[0]} vs {max_pair[1]})")
    print(f"Mean TM-score: {upper_vals.mean().mean():.3f}")

 

epitopes_df = pd.read_csv(Path("..", "data", "benchmark_dataset_epitopes.tsv"), sep="\t")
 
epitope_data = load_all_epitope_data(epitopes_df)
 
#tm_df, rmsd_df = compute_tmscore_matrix(epitope_data)
tm_df = compute_tmscore_matrix(epitope_data) 

tm_df.to_csv("figures/tmscore_matrix.tsv", sep="\t")
#rmsd_df.to_csv("figures/rmsd_matrix.tsv", sep="\t")
print("Saved tmscore_matrix.tsv")

plot_heatmap(tm_df, out_path="figures/tmscore_heatmap.png")

