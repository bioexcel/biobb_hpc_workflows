# BioExcel Building Blocks (BioBB) High Performance Computing (HPC) Workflow repository


We are working on some new workflows that will be added to the current list of workflows and launchers:

* Workflows for MD:
    * [md_list](workflows/MD/md_list.py): Performs a system setup and runs a molecular dynamics simulation on each one of the structures listed in the YAML properties file.
    * [md_muts_sets](workflows/MD/md_muts_sets.py): Performs a system setup and runs a molecular dynamics simulation for each one of the listed mutations in a given structure.
    * [md_add_muts_wt](workflows/MD/md_add_muts_wt.py): Applies a list of mutations over the initial structure obtaining a set of structures (initial structure + one mutation, initial structure + two mutations, initial structure + three mutations, ...). Finally performs a system setup and runs a molecular dynamics simulation for each of the structures in the set.

* Workflows for Free Energy Calculation (PMX):
    * [pmx_cv_cufix_term](workflows/PMX/pmx_cv_cufix_term.py): Performs a fast-growth mutation free energy calculation from two equilibrium trajectories.

* Launchers:
    * [md_launch](launchers/md_launch.py): Launcher for the [md_list](workflows/MD/md_list.py) workflow.
    * [mdmut_launch](launchers/mdmut_launch.py): Launcher for the [md_muts_sets](workflows/MD/md_muts_sets.py) and [md_add_muts_wt](workflows/MD/md_add_muts_wt.py) workflows.
    * [pmx_launch](launchers/pmx_launch.py): Launcher for the [pmx_cv_cufix_term](workflows/PMX/pmx_cv_cufix_term.py).

The launchers are available in:
    * Mare Nostrum IV (MN4) at the Barcelona Supercomputing Center (BSC) (`/apps/BIOBB/workflows`)

### Installation

* Download and install BioBB HPC Conda Pack (BioBB library + PyCOMPSs Workflow manager):<br>
https://mmb.irbbarcelona.org/biobb/availability/condapacks

* Clone GitHub repository:<br>
```git clone https://github.com/bioexcel/biobb_hpc_workflows.git```

* Generate and modify HPC cluster/supercomputer configuration file from template in [sc_conf](sc_conf/sc_conf.git.yml):<br>
```cp sc_conf/sc_conf.git.yml sc_conf/sc_conf.yaml```<br>
The template should be modified with specific HPC cluster/supercomputer paths. In particular, these two important fields should be added:<br>   
   > 
   > **workflows_path**: Path to this repo in the cluster file system.
   > 
   > **biobb_path**: Path to the BioBB Conda Pack in the cluster file system. 
   > 

### Testing

* Different tests can be found in the ([README.doc](README.doc)) file.
* Test data can be found in https://mmb.irbbarcelona.org/biobb/availability/condapacks

### Biobb modules used

* [biobb_pmx](https://github.com/bioexcel/biobb_pmx): Tools to setup and run Alchemical Free Energy calculations.
* [biobb_md](https://github.com/bioexcel/biobb_md): Tools to setup and run Molecular Dynamics simulations.
* [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.
* [biobb_structure_utils](https://github.com/bioexcel/biobb_structure_utils): Tools to  modify or extract information from a PDB structure file.

## Copyright & Licensing
This software has been developed in the [MMB group](http://mmb.irbbarcelona.org) at the [BSC](http://www.bsc.es/) & [IRB](https://www.irbbarcelona.org/) for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU H2020 [823830](http://cordis.europa.eu/projects/823830), EU H2020 [675728](http://cordis.europa.eu/projects/675728)).

* (c) 2015-2021 [Barcelona Supercomputing Center](https://www.bsc.es/)
* (c) 2015-2021 [Institute for Research in Biomedicine](https://www.irbbarcelona.org/)

Licensed under the
[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), see the file LICENSE for details.

![](https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png "Bioexcel")
