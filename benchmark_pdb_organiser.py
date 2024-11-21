###We organise the pdbs in the Github folder

import os
import sys
import pandas as pd
from pathlib import Path
import re
import shutil
import numpy as np

DATADIR = "/trinity/login/msanchez/MiguelSM_2024"

nanobodies = pd.read_csv(f"{DATADIR}/nanobodies_dataset/nanobodies_sabdab_unseen_resolution_checkseen_allfilters_07-02-2024_cdr.tsv", sep="\t")


for n in range(nanobodies.shape[0]):
    pdb, chain = nanobodies.iloc[n]["pdb"], nanobodies.iloc[n]["chain"]

    # Create the folder for each pdb
    os.makedirs(Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}"), exist_ok=True)
    #we create a folder for each of the ensembles inside each pdb
    os.makedirs(Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IB_ensemble"), exist_ok=True)
    os.makedirs(Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IBMo_ensemble"), exist_ok=True)
    os.makedirs(Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IBMu_ensemble"), exist_ok=True)
    os.makedirs(Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IBMM_ensemble"), exist_ok=True)

    # Copy the pdb files to the folder
    for file in os.listdir(Path(DATADIR, "nanobodies_dataset/haddock_runner/imbuilder/bound")):
        if file.startswith(f"{pdb+chain}_r_u") and file.endswith(".pdb"):
            shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/imbuilder/bound/{file}"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IB_ensemble/{file}"))
    
    for file in os.listdir(Path(DATADIR, "nanobodies_dataset/haddock_runner/alphafold_monomer_imbuilder/bound")):
        if file.startswith(f"{pdb+chain}_r_u") and file.endswith(".pdb"):
            shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_monomer_imbuilder/bound/{file}"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IBMo_ensemble/{file}"))
    
    for file in os.listdir(Path(DATADIR, "nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound")):
        if file.startswith(f"{pdb+chain}_r_u") and file.endswith(".pdb"):
            shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{file}"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IBMu_ensemble/{file}"))
    
    for file in os.listdir(Path(DATADIR, "nanobodies_dataset/haddock_runner/alphafold_multimer_monomer_imbuilder/bound")):
        if file.startswith(f"{pdb+chain}_r_u") and file.endswith(".pdb"):
            shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_monomer_imbuilder/bound/{file}"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/IBMM_ensemble/{file}"))
    
    #we copy the bound receptor files
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/pdb_files_sequence/{pdb}-{chain}.pdb"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_r_b.pdb"))

    #now we copy the bound ligand files
    for file in os.listdir(Path(DATADIR, "nanobodies_dataset/haddock_runner/imbuilder/bound")):
        if file.startswith(f"{pdb+chain}_l_u") and file.endswith(".pdb"):
            shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/imbuilder/bound/{file}"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_l_b.pdb"))
    
    #we copy the unbound ligand files
    for file in os.listdir(Path(DATADIR, "nanobodies_dataset/haddock_runner/imbuilder/semibound")):
        if file.startswith(f"{pdb+chain}_l_u") and file.endswith(".pdb"):
            shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/imbuilder/semibound/{file}"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_l_u.pdb"))
    
    #we copy the reference files
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/imbuilder/bound/{pdb}{chain}_ref.pdb"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_ref.pdb"))

    #we copy the ambiguous restraint files
    #Real interface
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_real_ambig.tbl"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_real_ambig.tbl"))
    #Loose interface
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_loose_ambig.tbl"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_loose_ambig.tbl"))
    #Two-Hit interface
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_twohit_ambig.tbl"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_twohit_ambig.tbl"))
    #All surface
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_loose_all_ambig.tbl"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_allsurface_ambig.tbl"))
    #Loose Mix
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_mix.tbl.tgz"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_loose_mix_ambig.tbl.tgz"))
    #Two-Hit Mix
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_mix_2h.tbl.tgz"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_twohit_mix_ambig.tbl.tgz"))

    ###We copy the unambiguous restraint files for the bound antigen
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/bound/{pdb}{chain}_unambig.tbl"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_l_b_unambig.tbl"))
    #and for the unbound antigen
    shutil.copy(Path(DATADIR, f"nanobodies_dataset/haddock_runner/alphafold_multimer_imbuilder/semibound/{pdb}{chain}_unambig.tbl"), Path(DATADIR, f"github/Nanobody_antigen_benchmark_dataset/benchmark_pdb_files/{pdb}{chain}/{pdb+chain}_l_u_unambig.tbl"))

    
            