import os
import pandas as pd
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import cdist
from itertools import combinations


def keep_atoms(inp_pdb_f):
    """
    Keep atoms.

    Parameters
    ----------
    inp_pdb_f : Path
        Path to PDB file.

    Returns
    -------
    out_pdb_fname : Path
        Path to PDB file.
    """
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-atoms.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in pdb_fh:
                if line.startswith("ATOM"):
                    f.write(line)
    return out_pdb_fname


def load_contacts_mda(pdb_f, cutoff=4.5): #distance cutoff of 4.5 Angstroms
        """
        Load residue-based contacts.
        """
        # get chains
        mdu = mda.Universe(pdb_f)
        con_list, con_num_list = [], []
        unique_chains = [str(el)[-2] for el in mdu.segments]
        unique_chain_combs = list(combinations(unique_chains, 2))
        # calculating contacts
        for pair in unique_chain_combs:
                chain_n = mdu.select_atoms(f"chainID {pair[0]} and not name H*")
                chain_s = mdu.select_atoms(f"chainID {pair[1]} and not name H*")
                dist = cdist(chain_n.positions, chain_s.positions)
                npw = np.where(dist < cutoff)
                del dist
                for k in range(npw[0].shape[0]):
                        con = (pair[0], f"{chain_n.atoms[npw[0][k]].resid}{chain_n.atoms[npw[0][k]].icode}", pair[1], f"{chain_s.atoms[npw[1][k]].resid}{chain_s.atoms[npw[1][k]].icode}")
                        con_num = (pair[0], chain_n.atoms[npw[0][k]].resid, pair[1], chain_s.atoms[npw[1][k]].resid)
                        #print(f"con {con}")
                        con_list.append(con)
                        con_num_list.append(con_num)
        return [set(con_list), set(con_num_list)]

nanobodies = pd.read_csv(Path("..", "data", "nanobodies_interface_analysis_dataset_defined_unbiased.tsv"), sep="\t")
summary_sabdab = pd.read_csv(Path("..", "data", "sabdab_nano_summary_all.tsv"), sep="\t")

#We will use the chothia numbering scheme to analyse the location of the interaction residues

per_pdb_interacting_residues = []
per_pdb_interacting_resnums = []

per_pdb_interacting_regions, per_pdb_loose_interaction_regions = [], []
per_pdb_interacting_frameworks, per_pdb_interacting_conserved_frameworks = [], []

for n in range(nanobodies.shape[0]):
    pdb, chain = nanobodies.iloc[n]["pdb"], nanobodies.iloc[n]["chain"]

    #First we get the chothia pdb file
    chothia_nanobody = Path("..", "paratope_analysis_chothia_pdb_files", f"{pdb}.pdb")
    chothia_atoms = keep_atoms(chothia_nanobody)
    # #We calculate the real interface scenario
    mda_contacts = load_contacts_mda(chothia_atoms)
    os.remove(chothia_atoms)
    contact = mda_contacts[0]
    contact_nums = mda_contacts[1]
    
    inter_resids = []
    for con in contact:
        if con[0] == chain:
            inter_resids.append(con[1])
        elif con[2] == chain:
            inter_resids.append(con[3])
    
    unique_inter_resids = list(set(inter_resids))
    print(unique_inter_resids)

    inter_resnums = []
    for con in contact_nums:
        if con[0] == chain:
            inter_resnums.append(con[1])
        elif con[2] == chain:
            inter_resnums.append(con[3])
    
    unique_inter_resnums = list(set(inter_resnums))
    print(unique_inter_resnums)
    
    per_pdb_interacting_residues.append([pdb, chain, unique_inter_resids, inter_resids])
    per_pdb_interacting_resnums.append([pdb, chain, "-".join([str(res) for res in sorted(unique_inter_resnums)]), "-".join([str(res) for res in sorted(inter_resnums)])])

#Now we convert the list into dataframes
per_pdb_interacting_residues_df = pd.DataFrame(per_pdb_interacting_residues, columns=["pdb", "chain", "unique_interaction_residues", "interacting_residues"])
per_pdb_interacting_resnums_df = pd.DataFrame(per_pdb_interacting_resnums, columns=["pdb", "chain", "unique_interaction_resnums","interacting_resnums"])

#We will save the dataframes
per_pdb_interacting_residues_df.to_csv(Path("..", "data", "per_pdb_interacting_residues_unbiased_new.tsv"), sep="\t", index=False)
per_pdb_interacting_resnums_df.to_csv(Path("..", "data", "per_pdb_interacting_resnums_unbiased_new.tsv"), sep="\t", index=False)

nanobodies_angles = pd.read_csv(Path("..", "data", "nanobodies_torsion_angles.tsv"), sep="\t")

probability_regions = []
for n in range(per_pdb_interacting_resnums_df.shape[0]):
    pdb,chain = per_pdb_interacting_resnums_df.iloc[n]["pdb"], per_pdb_interacting_resnums_df.iloc[n]["chain"]

    angle_pdb = nanobodies_angles[nanobodies_angles["pdb"] == pdb]
    angle_class = angle_pdb[angle_pdb["chain"]==chain]["angle_class"].values[0]

    species_chain = summary_sabdab[summary_sabdab["Hchain"] == chain]
    species = species_chain[species_chain["pdb"]==pdb]["heavy_species"].values[0]

    fr1, fr2, fr3, fr4, cdr1, cdr2, cdr3 = 0, 0, 0, 0, 0, 0, 0
    res_no_icode = per_pdb_interacting_resnums_df.iloc[n]["unique_interaction_resnums"]
    res_no_icode = [int(el) for el in res_no_icode.split("-")]

    for res in res_no_icode:
        res = int(res)
        if res > 135:
            continue
        if res in range(1, 26):
            fr1 += 1
        elif res in range(26, 33): #chothia definition from 26 to 32
            cdr1 += 1
        elif res in range(33, 52):
            fr2 += 1
        elif res in range(52, 57): #chothia definition from 52 to 56
            cdr2 += 1
        elif res in range(57, 95):
            fr3 += 1
        elif res in range(95, 103):#chothia definition from 95 to 102
            cdr3 += 1
        elif res > 102:
            fr4 += 1

    probability_regions.append([pdb, chain, angle_class, fr1/len(res_no_icode), cdr1/len(res_no_icode), fr2/len(res_no_icode), cdr2/len(res_no_icode), fr3/len(res_no_icode), cdr3/len(res_no_icode), fr4/len(res_no_icode), species])

probability_regions_df = pd.DataFrame(probability_regions, columns=["pdb", "chain", "angle_class", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "species"])
print(probability_regions_df)
probability_regions_df.drop(["species"], axis=1, inplace=True)
print(probability_regions_df["angle_class"].value_counts())

#we drop the ones with angle_class other
probability_regions_df = probability_regions_df[probability_regions_df["angle_class"] != "other"]
probability_regions_df.to_csv(Path("..", "data", "probability_regions_interaction.tsv"), sep="\t", index=False, float_format="%.3f")