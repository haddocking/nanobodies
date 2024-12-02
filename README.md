# Nanobody - Antigen HADDOCK docking

This repository conatins all data and analysis required to reproduce the figures and supplementary figures in the 'Information-driven docking of machine learning structure ensembles achieves high success on nanobody-antigen complexes' Major Research Project Report from Miguel Sánchez Marín (MSc Bioinformatics and Biocomplexity, Utrecht University)

## Repository organisation
This repository contains the follwoing directories:
* **analysis**: Python scripts required to reproduce the report figures.
  * **figures/**: directory with the report figures in png format.
  - **alphafold_complex_assessment.py** : reads the files `data/alphafold*_complex_predictions.tsv` and `data/alphafold*_complex_predictions_dockq_top1.tsv` and plots the docking success rates (`analysis/figures/alphafold2multimer_af3_sr_comparison.png`) and DockQ distribution for Top1 models (`analysis/figures/alphafold2multimer_af3_dockq_top1_comparison.png`)
  - **dockq_features_correlation.py**
  - **flexref_fnat_improvement.py**
  - **haddock_docking_sr_scenario_analysis.py**
  - **paratope_analysis.py**
  - **single_structure_accuracy_assessment.py**

* **benchmark_pdb_files.tgz**: directory with all the input nanobody structures, ligand structures and restraint files used in the report. Due to size limitations the directory is compressed (uncompressed directory size 349 MB). In order to access, run `tar -xzvf benchmark_pdb_files.tgz`.
 
* **cfg_files**: example .cfg files for docking with HADDOCK3 for all different information scenarios.
  - **nanobody_antigen_all_surface_haddock_docking.cfg**
  - **nanobody_antigen_loose_interface_haddock_docking.cfg**
  - **nanobody_antigen_mix_loose_interface_haddock_docking.cfg**
  - **nanobody_antigen_mix_twohit_interface_haddock_docking.cfg**
  - **nanobody_antigen_real_interface_haddock_docking.cfg**
  - **nanobody_antigen_twohit_interface_haddock_docking.cfg**
  
* **data**: all data required to reproduce the figures in the report.
  - **alphafold2_monomer_immunebuilder_centroids_rmsd.tsv**
  - **alphafold2_rmsd.tsv**
  - **alphafold2multimer_complex_predictions_dockq_top1.tsv**
  - **alphafold2multimer_complex_predictions_sr.tsv**
  - **alphafold2multimer_immunebuilder_centroids_rmsd.tsv**
  - **alphafold2multimer_monomer_immunebuilder_centroids_rmsd.tsv**
  - **alphafold2multimer_rmsd.tsv**
  - **alphafold3_complex_predictions_dockq_top1.tsv**
  - **alphafold3_complex_predictions_sr.tsv**
  - **alphafold3_rmsd.tsv**
  - **alphahold2_monomer_immunebuilder_centroids_rmsd.tsv**
  - **cdr3_rmsd_dockq.tsv**
  - **clustered_paratope_interacting_residues.tsv**
  - **epitope_rmsd_dockq.tsv**
  - **fnat_diff_flexref.tsv**
  - **haddock_docking_clustered_sr_per_scenario.tsv**
  - **haddock_docking_mixed_paratope_restraints_sr.tsv**
  - **haddock_docking_sr_allsurface.tsv**
  - **haddock_docking_sr_per_scenario.tsv**
  - **immunebuilder_rmsd.tsv**
  - **nanonet_rmsd.tsv**
  - **paratope_represented_dockq.tsv**
  - **probability_regions_interaction.tsv**
  - **raptorxsingle_rmsd.tsv**
  - **unbound_antigens_rmsd.tsv**
  - **voro_scoring_clusters_sr.tsv**
  - **voro_scoring_models_sr.tsv**
