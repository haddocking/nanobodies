# Nanobody - Antigen HADDOCK docking

This repository contains all data and analysis required to reproduce the figures and supplementary figures in the 'Information-driven docking of machine learning structure ensembles achieves high success on nanobody-antigen complexes' Major Research Project Report from Miguel Sánchez Marín (MSc Bioinformatics and Biocomplexity, Utrecht University)

## Repository organisation
This repository contains the follwoing directories:
* **analysis**: Python scripts required to reproduce the report figures.
  * **figures/**: directory with the report figures in png format.
  - **alphafold_complex_assessment.py** : reads the files `data/alphafold*_complex_predictions.tsv` and `data/alphafold*_complex_predictions_dockq_top1.tsv` and plots the docking success rates (`analysis/figures/alphafold2multimer_af3_sr_comparison.png`) and DockQ distribution for Top1 models (`analysis/figures/alphafold2multimer_af3_dockq_top1_comparison.png`)
  - **dockq_features_correlation.py** : reads the files `data/*_dockq.tsv` and plots the correlation of top10 max dockQ and min CDR3 RMSD, min epitope RMSD and % of paratope represented in the restraints(`analysis/figures/dockq_*_correlation_top10.png`).
  - **flexref_fnat_improvement.py** : reads the file `data/fnat_diff_flexref.tsv` and plots the Fnat improvement distribution after flexible refinement stage (`analysis/figures/fnat_diff_flexref.png`).
  - **haddock_docking_sr_scenario_analysis.py** : reads the files `data/haddock_docking_sr_per_scenario.tsv`, `data/haddock_docking_clustered_sr_per_scenario.tsv` and `data/haddock_docking_sr_allsurface.tsv` and plots the docking success rates in `analysis/figures/haddock_docking_*.png`. It also reads `data/voro_scoring_*.tsv`and plots the success rates (`analysis/figures/voro_scoring_*.png`).
  - **paratope_analysis.py** : reads the file `data/probability_regions_interaction.tsv`and plots the regions frequency in the paratope of _kinked_ and _extended_ nanobodies (`analysis/figures/kinked_extended_paratope_region_freq.png`). Reads the file `data/clustered_paratope_interacting_residues.tsv` and plots the frequency of each residue in the different CDR and CDR+FR paratope representations (`analysis/figures/paratope_representations_residue_freq.png`). Finally, it reads `data/haddock_docking_mixed_paratope_restraints_sr.tsv`and plots the docking success rates for the new restraints (`analysis/figures/mixed_paratope_restraints_sr.png`).
  - **single_structure_accuracy_assessment.py**: reads all `data/*_rmsd.tsv` files, calculates the mean, standard deviation and median values for all regions RMSD (`analysis/figures/nanobody_rmsd_summary.tsv`), and plots the CDR3 RMSD distribution for top1 models (`analysis/figures/cdr3_rmsd_top1_prediction_comparison.png`) and for the best in each ensemble (`analysis/figures/cdr3_rmsd_best_ensemble_best_prediction_comparison.png`) and the RMSD distribution in the unbound antigens (`analysis/figures/rmsd_unbound_antigen_comparison.png`).

* **benchmark_pdb_files.tgz**: directory with all the input nanobody structures, ligand structures and restraint files used in the report. Due to size limitations the directory is compressed (uncompressed directory size 349 MB). In order to access, run `tar -xzvf benchmark_pdb_files.tgz`.
  * **{pdb}{chain}**: there is a directory for each PDB in the dataset containing:
   * IB*_ensemble/: directory containing the input nanobody ensemble from IB (ImmuneBuilder), IBMo (ImmuneBuilder + AlphaFold2), IBMu (ImmuneBuilder + AlphaFold2-Multimer) or IBMM (ImmuneBuilder + ALphaFold2 + AlphaFold2-Multimer) predictions.
      - {pdb}{chain}_r_u_ensemble.pdb: ensemble containing all nanobody structure predictions.
      - {pdb}{chain}_r_u_*.pdb: each individual predicted nanobody structure.
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

* **cfg_files**: example .cfg files for docking with HADDOCK3 for all different information scenarios.
  - **nanobody_antigen_all_surface_haddock_docking.cfg**: example .cfg file for All Surface docking.
  - **nanobody_antigen_loose_interface_haddock_docking.cfg**: example .cfg file for Loose Interface docking.
  - **nanobody_antigen_mix_loose_interface_haddock_docking.cfg**: example .cfg file for Mixed Loose Interface docking.
  - **nanobody_antigen_mix_twohit_interface_haddock_docking.cfg**: example .cfg file for Mixed Two-Hit Interface docking.
  - **nanobody_antigen_real_interface_haddock_docking.cfg**: example .cfg file for Real Interface docking.
  - **nanobody_antigen_twohit_interface_haddock_docking.cfg**: example .cfg file for Two-Hit Interface docking.
  
* **data**: all data required to reproduce the figures in the report.
  - *_rmsd.tsv: RMSD of the different nanbody regions (FR, CDR1, CDR2, CDR3) for the different methods. If ending with _centroids_rmsd.tsv it only represents CDR3 RMSD from the centroids of the clustered ensembles of predictions. 
  - *_complex_predictions_dockq_top1.tsv: DockQ values for top1 complex predictions from AlphaFold2-Multimer and AlphaFold3.
  - *_complex_predictions_sr.tsv: docking success rates for different CAPRI classes and depending of topN ranked structures from AlphaFold2-Multimer and AlphaFold3 predictions. 
  - cdr3_rmsd_dockq.tsv: minimum available CDR3 RMSD in the ensemble and maximum DockQ achieved in the top10 models for True Interface scenario, bound antigen EMRef stage.
  - clustered_paratope_interacting_residues.tsv: paratope residues from each Paratope Dataset PDB and the cluster to which they belong.
  - epitope_rmsd_dockq.tsv minimum epitope RMSD and maximum DockQ achieved in the top10 models for True Interface scenario, unbound antigen EMRef stage.
  - fnat_diff_flexref.tsv: frequency of fnat improvements after flexible refinement classified in 0.05 width bins.
  - haddock_docking_*.tsv: docking success rates for different CAPRI classes and depending of topN ranked models/clusters, for different information scenarios and/or pipeline stages.
  - paratope_represented_dockq.tsv : % of represented paratope in the restraints and maximum DockQ achieved in the top10 models for Loose and Two-Hit Interface scenario, bound antigen EMRef stage.
  - probability_regions_interaction.tsv : percent of kinked/extended structures with each nanobody region being present in the paratope.
  - voro_scoring_*_sr.tsv: docking success rates for different information scenarios as a function of the topN models/clusters after VoroIF-GNN scoring the final HADDOCK-scored models.

