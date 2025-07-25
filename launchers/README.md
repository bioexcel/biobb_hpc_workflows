# BioExcel Building Blocks (BioBB) High Performance Computing (HPC) Workflow Repository

Collection of pre-exascale biomolecular simulation workflows built using the BioExcel Building Blocks (BioBB) library and controlled by the PyCOMPSs workflow manager.

- **BioExcel Building Blocks (BioBB)**: [http://mmb.irbbarcelona.org/biobb/](http://mmb.irbbarcelona.org/biobb/)
- **PyCOMPSs**: [https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar)
- **BioBB HPC Conda Pack**: [https://mmb.irbbarcelona.org/biobb/availability/condapacks](https://mmb.irbbarcelona.org/biobb/availability/condapacks)

---

## 📁 Index

1. [Workflows Library](#1-workflows-library)
2. [Workflow Launchers](#2-workflow-launchers)
3. [Examples of Use](#3-examples-of-use)
4. [Biobb Modules Used](#4-biobb-modules-used)
5. [Copyright & Licensing](#5-copyright--licensing)
---

## 1) Workflows Library

### 🧬 Workflows for MD

- `md_list`: Setup and run MD on listed structures.
- `md_muts_sets`: Run MD for each mutation in a structure.
- `md_add_muts_wt`: Apply multiple mutations to a structure and run MD per variant.

### ⚛️ Free Energy Calculation (PMX)

- `pmx_cv_cufix_term`: Fast-growth mutation free energy calc from equilibrium trajectories.

---

## 2) Workflow Launchers

- `md_launch`: Launcher for `md_list`
- `mdmut_launch`: Launcher for `md_muts_sets` and `md_add_muts_wt`
- `pmx_launch`: Launcher for `pmx_cv_cufix_term`

---

## 3) Examples of Use

### 🔬 Example 1: Launching MD simulation for a single structure

1.1 The Lysozyme protein structure is used as a test case. The execution uses 4 nodes in MPI. The length of the simulation is 5ns and the force-field used is Charmm27.

```
python md_launch.py -i test_data/LISTS/charmm.list -q default -nn 4 -mpi 4 -jn charmm -t 120 -l 5 -o Charmm -f charmm27
```

### 🧬 Example 2: Launching MD simulations for a list of PDB structures

2.1 The miniABC list (https://doi.org/10.1093/nar/gkz905) is used as a test dataset. The list is formed by 13 nucleic acids eighteen-mers containing all possible tetramers. 14 nodes are used, 1 of them reserved for PyCOMPSs. The system box type, size and ionic concentration are given as parameters.

```
python md_launch.py -i test_data/LISTS/miniABC.list -q default -nn 14 -mpi 1 -jn miniABC -t 120 -l 5 -o miniABC -s DNA -bs 0.7 -bt octahedron --concentration 0.07
```

2.2 A reduced subset of the "prot-dna" list (https://doi.org/10.1016/j.jmb.2019.07.021) is used as a test dataset. The list is formed by 2 protein-DNA complexes from the original list of 50 complexes studied in the published work. In this case, 16 nodes are used, using 8 nodes per MD simulation in MPI.

```
python md_launch.py -i test_data/LISTS/prot-dna.list -q default -nn 16 -mpi 8 -jn protDNA -t 120 -l 5 -o protDNA --system Protein-DNA
```

### 🧬 Example 3: Mutations on Protein Structures

3.1 Modeling a particular mutation (Alanine 369 to Serine) identified between the Angiotensin-Converting Enzime (ACE) protein in human versus the same protein in the Rhinolophus affinis species (bat). Two MD simulations are run in the same job, one for the WT protein (without the mutation) and another one with the modeled mutation (A369S). 16 nodes are used, using 8 nodes x MD simulation in MPI.

```
python mdmut_launch.py -wt test_data/hACE.pdb -m WT+A:Ala369Ser -l 5 -nn 16 -mpi 8 -o MDs-hACE -jn MDsACE -t 120 -q default -sc_conf sc_conf.yml -sc mn5
```

3.2 Same example as before, but now using a custom force-field, manually modified (GROMACS specbond.dat file) to keep distance restrictions on a couple of delta Histidines with the structural ZN ion of the ACE protein.

```
python mdmut_launch.py -wt test_data/hACE-ZN.pdb -m WT+A:Ala369Ser -l 5 -nn 16 -mpi 8 -o MDs-hACE-ZN -jn MDsACE-ZN -t 120 -q default -gl test_data/ff_mod -f amber99sb-ildn-mod -sc_conf sc_conf.yml -sc mn5 
```

3.3 Same example as 1.5, but now using the SARS-CoV-2 RBD complexed with the human ACE protein. The 4 simulations (ACE_WT, ACE_A369S, RBD-ACE_WT and RBD-ACE_A369S) can be used for a free energy calculation using non-equilibrium approximations (fast growth thermodynamic integration), see example 4.

```
python mdmut_launch.py -wt test_data/RBD-hACE.pdb -m WT+A:Ala369Ser -l 5 -nn 16 -mpi 8 -o MDs-RBD-hACE -jn MDsCOMPLEX -t 360 -q default -sc_conf sc_conf.yml -sc mn5 --w_in_m_cpus 112
```

### ⚛️ Example 4: Free Energy Calculations (PMX)

Computing the effect of a particular residue mutation on the protein binding free energy using non-equilibrium approximation (fast growth thermodynamic integration).

4.1 Computing the alchemical free energy for the Alanine369Serine mutation on the human ACE protein using the equilibrium MD simulations of the WT and the A369S mutation (see example 3.1). 32 nodes are used to compute the 100 x 2 (forward + reverse) short (50ps) thermodynamic integration simulations needed to extract the final deltaG. The total number of snapshots included in each of the equilibrium simulations (mut_start_end_num_frames, wt_start_end_num_frames), as well as the time of the last frame to read from these input simulations (wt_end, mut_end, in ps) are passed as parameters.

```
python pmx_launch.py -m Ala369Ser -wt_top test_data/CV/hACE/WT/gppmd.tpr -wt_trj test_data/CV/hACE/WT/md.xtc -mut_top test_data/CV/hACE/A369S/gppmd.tpr -mut_trj test_data/CV/hACE/A369S/md.xtc -i 2 -q default -t 350 -nn 32 -nf 100 -fe 50 --mut_start_end_num_frames 2500 --wt_start_end_num_frames 2500 --wt_end 5000 --mut_end 5000 -jn PMX_hACE -o PMX-hACE
```

4.2 Same example as the previous one, but now computing the alchemical free energy for the Alanine369Serine mutation on the SARS-CoV-2 Receptor Binding Domain (RBD) complexed with the human ACE (hACE) protein. The combination of the deltaG computed in the example 4.1 with the one computed in this example is giving the final binding free energy for the A369S mutation on the protein complex.

```
python pmx_launch.py -m Ala369Ser -wt_top test_data/CV/RBD-hACE/WT/gppmd.tpr -wt_trj test_data/CV/RBD-hACE/WT/md.xtc -mut_top test_data/CV/RBD-hACE/A369S/gppmd.tpr -mut_trj test_data/CV/RBD-hACE/A369S/md.xtc -i 2 -q default -t 350 -nn 32 -nf 100 -fe 50 --mut_start_end_num_frames 2500 --wt_start_end_num_frames 2500 --wt_end 5000 --mut_end 5000 -jn PMX_COMPLEX -o PMX-RBD-hACE
```

## 4) Biobb Modules Used

| Module                  | Description                                     | Link                                                                 |
|------------------------|-------------------------------------------------|----------------------------------------------------------------------|
| biobb_pmx              | Tools to set up and run alchemical free energy calculations | https://github.com/bioexcel/biobb_pmx            |
| biobb_md               | Tools to set up and run Molecular Dynamics simulations       | https://github.com/bioexcel/biobb_md             |
| biobb_analysis         | Tools to analyze Molecular Dynamics trajectories            | https://github.com/bioexcel/biobb_analysis       |
| biobb_structure_utils  | Tools to modify or extract info from PDB files               | https://github.com/bioexcel/biobb_structure_utils |

---

## 5) ⚖️ Copyright & Licensing

This software has been developed in the MMB group at the BSC & IRB for the European BioExcel project, funded by the European Commission (EU H2020 823830, EU H2020 675728).

© 2015–2021 Barcelona Supercomputing Center  
© 2015–2021 Institute for Research in Biomedicine

Licensed under the Apache License 2.0 — see the LICENSE file for details.
