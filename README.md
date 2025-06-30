# Nanobody - Antigen HADDOCK docking

This repository contains all data and analysis required to reproduce the figures and supplementary figures in the paper 'Combining AI structure prediction and integrative
modelling for nanobody-antigen complexes' from Miguel Sánchez-Marín, Marco Giulini, Alexandre MJJ Bonvin

## Repository organisation
This repository contains the follwoing directories:
* **analysis**: Python scripts required to reproduce the paper figures.
  * **figures/**: directory with the report figures in png format.
  - **af_metrics_vs_sr.py** : reads the success rates and metrics for both Alphafold2-Multimer and AlphaFold3 and compares them to the docking success rates, highlighting the cases when HADDOCK docking was able to "rescue" a failed AlphaFold run.
  - **dockq_features_correlation.py** : reads the files `data/*_dockq.tsv` and performs several analysis with the correlation between the max dockQ in the top10 models and several quantities, such as the minimum CDR3 RMSD across the ensemble (see ).
  - **flexref_fnat_improvement.py** : reads the file `data/fnat_diff_flexref.tsv` and plots the Fnat improvement distribution after flexible refinement stage (`analysis/figures/figure4.png`).
  - **haddock_docking_sr_scenario_analysis.py** : reads the files `data/haddock_docking_sr_per_scenario.tsv`, and `data/haddock_docking_clustered_sr_per_scenario.tsv` and plots the docking success rates.
  - **paratope_analysis.py** : reads the file `data/probability_regions_interaction.tsv`and plots the regions frequency in the paratope of _kinked_ and _extended_ nanobodies (`figures/SI_figure5.png`). Reads the file `data/clustered_paratope_interacting_residues.tsv` and plots the frequency of each residue in the different CDR and CDR+FR paratope representations (`analysis/figures/figure5.png`). Finally, it calculates and plots the docking success rates for these new restraints (`analysis/figures/mixed_paratope_restraints_sr.png`).
  - **single_structure_assessment.py**: does all the analysis about the monomeric nanobody structure, including calculating the mean, standard deviation and median values for all regions RMSD (`analysis/figures/nanobody_rmsd_summary.tsv`), plotting the CDR3 RMSD distribution for top1 models (`analysis/figures/cdr3_rmsd_top1_prediction_comparison.png`) and for the best in each ensemble (`analysis/figures/cdr3_rmsd_best_ensemble_best_prediction_comparison.png`). The script calculates also the RMSD distribution for the unbound antigens (`analysis/figures/rmsd_unbound_antigen_comparison.png`).

* **benchmark_pdb_files**: directory with all the input nanobody structures, ligand structures and restraint files used in the report. Quality metrics of each obtained HADDOCK model are also reported. 
  * **{pdb}{chain}**: there is a directory for each PDB in the dataset containing:
   * *ensemble.pdb: the input nanobody ensemble from IB (ImmuneBuilder), IBMo (ImmuneBuilder + AlphaFold2), IBMu (ImmuneBuilder + AlphaFold2-Multimer) or IBMM (ImmuneBuilder + ALphaFold2 + AlphaFold2-Multimer) predictions.
    - {pdb}{chain}_allsurface_ambig.tbl: input ambiguous restraints for All Surface scenario.
    - {pdb}{chain}_l_u_unambig.tbl: input unambiguous restraints for unbound ligand structure.
    - {pdb}{chain}_real_ambig.tbl: input ambiguous restraints for Real Interface scenario.
    - {pdb}{chain}_l_b.pdb: input bound antigen structure.
    - {pdb}{chain}_loose_ambig.tbl: input ambiguous restraints for Loose Interface scenario.
    - {pdb}{chain}_ref.pdb: input reference structure.      
    - {pdb}{chain}_l_b_unambig.tbl: input unambiguous restraints for bound ligand structure.
    - {pdb}{chain}_loose_mix_ambig.tbl.tgz: input ambiguous restraints for Mixed paratope Loose Interface scenario.
    - {pdb}{chain}_twohit_ambig.tbl: input ambiguous restraints for Two-Hit Interface scenario.
    - {pdb}{chain}_l_u.pdb: input unbound antigen structure.
    - {pdb}{chain}_r_b.pdb: bound nanobody structure.
    - {pdb}{chain}_twohit_mix_ambig.tbl.tgz: input ambiguous restraints for Mixed paratope Two-Hit Interface scenario.
    - {pdb}{chain}_capri.tsv: CAPRI quality metrics for all the generated HADDOCK models.
    -  
* **cfg_files**: example .cfg files for docking with HADDOCK3 for all different information scenarios.
  - **nanobody_antigen_all_surface_haddock_docking.cfg**: example .cfg file for All Surface docking.
  - **nanobody_antigen_loose_interface_haddock_docking.cfg**: example .cfg file for Loose Interface docking.
  - **nanobody_antigen_mix_loose_interface_haddock_docking.cfg**: example .cfg file for Mixed Loose Interface docking.
  - **nanobody_antigen_mix_twohit_interface_haddock_docking.cfg**: example .cfg file for Mixed Two-Hit Interface docking.
  - **nanobody_antigen_real_interface_haddock_docking.cfg**: example .cfg file for Real Interface docking.
  - **nanobody_antigen_twohit_interface_haddock_docking.cfg**: example .cfg file for Two-Hit Interface docking.
  
* **data**: all data required to reproduce the figures in the report.
  - *_rmsd.tsv: RMSD of the different nanbody regions (FR, CDR1, CDR2, CDR3) for the different methods. If ending with _centroids_rmsd.tsv it only represents CDR3 RMSD from the centroids of the clustered ensembles of predictions. 
  - alphafold2* and alphafold3* : all the relevant statistics about AlphaFold2-Multimer and AlphaFold3 models.
  - cdr3_rmsd_dockq.tsv: minimum available CDR3 RMSD in the ensemble and maximum DockQ achieved in the top10 models
  - clustered_paratope_interacting_residues.tsv: paratope residues from each Paratope Dataset PDB and the cluster to which they belong.
  - epitope_rmsd_dockq.tsv minimum epitope RMSD and maximum DockQ achieved in the top10 models
  - fnat_diff_flexref.tsv: frequency of fnat improvements after flexible refinement classified in 0.05 width bins.
  - haddock_docking_*.tsv: docking success rates for different CAPRI classes and depending of topN ranked models/clusters, for different information scenarios and/or pipeline stages.
  - paratope_represented_dockq.tsv : % of represented paratope in the restraints and maximum DockQ achieved in the top10 models for Loose and Two-Hit Interface scenario, bound antigen EMRef stage.
  - probability_regions_interaction.tsv : percent of kinked/extended structures with each nanobody region being present in the paratope.
  - voro_scoring_*_sr.tsv: docking success rates for different information scenarios as a function of the topN models/clusters after VoroIF-GNN scoring the final HADDOCK-scored models.

