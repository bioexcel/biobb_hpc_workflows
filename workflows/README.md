# BioExcel HPC Workflows

Welcome to the **BioExcel HPC Workflows** collection — a modular and scalable set of pipelines designed for high-performance molecular dynamics simulations and free energy calculations. These workflows are built using [BioBB (BioExcel Building Blocks)](https://biobb-docs.readthedocs.io/en/latest/) and optimized for execution in **HPC environments** like MareNostrum at BSC.

## Folder Structure

This `workflows` directory contains several pipelines templates organized by analysis type:

- MD:
    - `md_list`: Molecular dynamics simulations (minimization, equilibration, production)
    - `md_add_muts_wt, md_muts_sets`: MD workflows incorporating residue mutations

- PMX:
    - `pmx_prot_lig, pmx_cv_cufix_term`: Alchemical free energy calculation workflows

Each workflow consists of a pair of files:

- **Python file** with the **workflow logic**
- **Yaml file** with the **workflow parameters**

These templates are used by the **biobb_hpc_workflows** launchers to generate job scripts adapted to specific supercomputers. See global [README](../README.md) file of the repo.

## Copyright & Licensing

This software has been developed in the MMB group at the BSC & IRB for the European BioExcel project, funded by the European Commission (EU Horizon Europe 101093290, EU H2020 823830, EU H2020 675728).

© 2015–2025 Barcelona Supercomputing Center  
© 2015–2025 Institute for Research in Biomedicine

Licensed under the Apache License 2.0 — see the LICENSE file for details.
