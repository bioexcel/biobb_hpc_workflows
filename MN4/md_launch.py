#!/usr/bin/env python3

from pathlib import Path
import argparse
import re
import shutil
import subprocess
import oyaml as yaml

def get_template_config_dict(config_yaml_path):
    with open(config_yaml_path) as config_yaml_file:
        return yaml.safe_load(config_yaml_file)


def launch(input_structures, queue, num_nodes, compss_version, md_length,
           base_dir, compss_debug, time, output_dir, job_name, mpi_nodes,
           gmxlib, system, ff, concentration, box_type, box_size):

    base_dir = Path(base_dir)
    #pth_path = Path.home().joinpath('.local', 'lib', 'python3.6', 'site-packages', 'biobb.pth')

    template_py_path = base_dir.joinpath('workflows', 'MD', 'md_list.py')
    template_yaml_path = base_dir.joinpath('workflows', 'MD', 'md_list.yaml')

    # Create working dir path
    work_dir = Path('.')
    working_dir_path = work_dir.joinpath('MDs', 'md_set')

    if output_dir:
        if output_dir.startswith('/'):
            working_dir_path = Path(output_dir).resolve()
        else:
            working_dir_path = work_dir.joinpath('MDs', output_dir)
    working_dir_path.mkdir(parents=True, exist_ok=True)

    # Check if it's the first launch
    run_number = 0
    run_dir = working_dir_path.joinpath("wf_md")
    config_yaml_path = working_dir_path.joinpath(f"md.yaml")
    wf_py_path = working_dir_path.joinpath(f"md.py")
    launch_path = working_dir_path.joinpath(f"launch.sh")
    if not job_name:
        job_name = "mdlaunch_job"
    while run_dir.exists():
        run_number += 1
        run_dir = working_dir_path.joinpath(f"wf_md_{str(run_number)}")
        config_yaml_path = working_dir_path.joinpath(f"md_{str(run_number)}.yaml")
        wf_py_path = working_dir_path.joinpath(f"md_{str(run_number)}.py")
        launch_path = working_dir_path.joinpath(f"launch_{str(run_number)}.sh")
        job_name = f"{job_name}_{str(run_number)}"

    # Copy py file
    shutil.copyfile(template_py_path, wf_py_path)

    # Read yaml template file
    config_dict = get_template_config_dict(template_yaml_path)

    # Update config_dict
    config_dict['working_dir_path'] = str(run_dir)

    # Length of the simulations
    config_dict['step13_grompp_md']['properties']['mdp']['nsteps'] = int((md_length*1000)/0.002)

    # Custom (modified) Force Field (initializing GMXLIB env variable)
    config_dict['step1_pdb2gmx']['properties']['gmxlib'] = gmxlib
    config_dict['step4_grompp_genion']['properties']['gmxlib'] = gmxlib
    config_dict['step6_grompp_min']['properties']['gmxlib'] = gmxlib
    config_dict['step9_grompp_nvt']['properties']['gmxlib'] = gmxlib
    config_dict['step11_grompp_npt']['properties']['gmxlib'] = gmxlib
    config_dict['step13_grompp_md']['properties']['gmxlib'] = gmxlib

    # Force Field
    config_dict['step1_pdb2gmx']['properties']['force_field'] = ff

    # Ionic concentration
    config_dict['step5_genion']['properties']['concentration'] = concentration

    # Box Type and Size
    config_dict['step2_editconf']['properties']['box_type'] = box_type
    config_dict['step2_editconf']['properties']['distance_to_molecule'] = box_size

    if (system == "DNA"):
        config_dict['step9_grompp_nvt']['properties']['mdp']['tc_grps'] =  "DNA Water_and_ions"
        config_dict['step11_grompp_npt']['properties']['mdp']['tc_grps'] =  "DNA Water_and_ions"
        config_dict['step13_grompp_md']['properties']['mdp']['tc_grps'] =  "DNA Water_and_ions"
    elif (system == "Protein-DNA"):
        config_dict['step8_make_ndx']['properties']['selection'] =  "\\\"Protein\\\"|\\\"DNA\\\""
        config_dict['step9_grompp_nvt']['properties']['mdp']['tc_grps'] =  "Protein_DNA Water_and_ions"
        config_dict['step11_grompp_npt']['properties']['mdp']['tc_grps'] =  "Protein_DNA Water_and_ions"
        config_dict['step13_grompp_md']['properties']['mdp']['tc_grps'] =  "Protein_DNA Water_and_ions"

    # Getting input structures to be simulated
    structures = ''
    for line in open(input_structures, 'r'):
        if not structures:
            structures = line.rstrip()
        else:
            structures = structures + "," + line.rstrip()

    config_dict['input_structures'] = structures
    print(config_dict['input_structures'])

    with open(config_yaml_path, 'w') as config_yaml_file:
        config_yaml_file.write(yaml.dump(config_dict))

    # Create launch
    with open(launch_path, 'w') as launch_file:
        launch_file.write(f"#!/bin/bash\n")
        launch_file.write(f"\n")
        launch_file.write(f"module purge\n")
        launch_file.write(f"\n")
        launch_file.write(f"module load ANACONDA/2019.10\n")
        #launch_file.write(f"source activate biobb\n")
        launch_file.write(f"module load biobb/covid_dev\n")
        launch_file.write(f"\n")
        #launch_file.write(f"# COMPSs environment\n")
        #launch_file.write(f"export COMPSS_PYTHON_VERSION=none\n")
        #launch_file.write(f"# COMPSs release\n")
        #launch_file.write(f"module load COMPSs/{compss_version}\n")
        #launch_file.write(f"\n")
        #launch_file.write(f"# Singularity\n")
        #launch_file.write(f"module load singularity\n")
        #launch_file.write(f"\n")
        #launch_file.write(f"#GROMACS 2019\n")
        #launch_file.write(f"module load intel/2018.4 impi/2018.4 mkl/2018.4 gromacs/2019.1\n")
        #launch_file.write(f"\n")
        #launch_file.write(f"#TASK_TIME_OUT env var\n")
        #launch_file.write(f"#3600 seconds = 1h\n")
        #launch_file.write(f"export TASK_TIME_OUT=\"3600\"\n")
        launch_file.write(f"#MULTINODE MPI env var\n")
        launch_file.write(f"export TASK_COMPUTING_NODES=\"{mpi_nodes}\"\n")
        launch_file.write(f"\n")
        launch_file.write(f"#Permissions for everyone\n")
        launch_file.write(f"umask ugo+rwx\n")
        launch_file.write(f"\n")
        launch_file.write(f"enqueue_compss ")
        if compss_debug:
            launch_file.write(f"-d ")
        if num_nodes == 1:
            launch_file.write(f"--worker_in_master_cpus=48 ")
        launch_file.write(f"--job_name={job_name} \
                          --num_nodes={num_nodes} \
                          --exec_time={str(time)} \
                          --base_log_dir=$PWD \
                          --worker_working_dir=$PWD \
                          --master_working_dir=$PWD \
                          --network=ethernet \
                          --qos={queue}  \
                          {wf_py_path} \
                          --config {config_yaml_path} ")
        launch_file.write(f"\n")

    subprocess.call(f"bash {launch_path}", shell=True)


def main():
    parser = argparse.ArgumentParser(description="Workflow to setup and run MD simulations for a set of PDB structures.")
    parser.add_argument('-i', '--input_structures', required=True, default='input.lst', type=str, help="(input.lst) [File with list of structures to be simulated, one structure path per line]")
    parser.add_argument('-s', '--system', required=False, default='Protein', type=str, help="(Protein) Molecule type [Protein|DNA|Protein-DNA]")
    parser.add_argument('-f', '--force_field', required=False, default='amber99sb-ildn', type=str, help="(amber99sb-ildn) Force Field to be used in the simulation [amber99sb-ildn|charmm27|gromos54a7|oplsaa]")
    parser.add_argument('-bs', '--box_size', required=False, default=0.8, type=float, help="(0.8) System box size (nanometers)")
    parser.add_argument('-bt', '--box_type', required=False, default='octahedron', type=str, help="(octahedron) System box type (triclinic|cubic|dodecahedron|octahedron)")
    parser.add_argument('-c', '--concentration', required=False, default=0.05, type=float, help="(0.05) System ionic concentration (mol/liter)")
    parser.add_argument('-q', '--queue', required=False, default='bsc_ls', type=str, help="(bsc_ls) [bsc_ls|debug]")
    parser.add_argument('-t', '--time', required=False, default=120, type=int, help="(120) [integer] Time in minutes")
    parser.add_argument('-nn', '--num_nodes', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-cv', '--compss_version', required=False, default='2.6.1', type=str, help="(2.6.1) [version_name]")
    parser.add_argument('-d', '--compss_debug', required=False, help="Compss debug mode", action='store_true')
    parser.add_argument('-l', '--md_length', required=False, default=10, type=int, help="(10ns) [integer] MD length in nanoseconds")
    parser.add_argument('-mpi', '--mpi_nodes', required=False, default=1, type=int, help="(1) [integer] Number of MPI nodes to be used per MD simulation")
    parser.add_argument('--base_dir', required=False, default='/apps/BIOBB/workflows', type=str, help="('.') [path_to_base_dir]")
    parser.add_argument('-o', '--output_dir', required=False, default='', type=str, help="Output dir name: If output_dir is absolute it will be respected if it's a relative path: /base_dir/MDs/runs/output_dir', if output_dir not exists, the name is autogenerated.")
    parser.add_argument('-jn', '--job_name', required=False, default='', type=str, help="Job name if it not exists, the name is autogenerated.")
    parser.add_argument('-gl', '--gromacs_lib', required=False, default='.', type=str, help="Gromacs lib, path where to find the force field libraries.")
    args = parser.parse_args()

    # ALL possible force fields
    # 1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
    # 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
    # 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
    # 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
    # 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
    # 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
    # 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
    # 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
    # 9: GROMOS96 43a1 force field
    #10: GROMOS96 43a2 force field (improved alkane dihedrals)
    #11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
    #12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
    #13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
    #14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
    #15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

    # Specific call of each building block
    launch(input_structures=args.input_structures,
           queue=args.queue,
           time=args.time,
           num_nodes=args.num_nodes,
           compss_version=args.compss_version,
           compss_debug=args.compss_debug,
           md_length=args.md_length,
           mpi_nodes=args.mpi_nodes,
           output_dir=args.output_dir,
           job_name=args.job_name,
           gmxlib=args.gromacs_lib,
           system=args.system,
           ff=args.force_field,
           box_type=args.box_type,
           box_size=args.box_size,
           concentration=args.concentration,
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
