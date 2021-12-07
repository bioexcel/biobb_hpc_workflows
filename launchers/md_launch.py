#!/usr/bin/env python3

from pathlib import Path
import argparse
import re, os, sys
import shutil
import subprocess
import oyaml as yaml

def get_template_config_dict(config_yaml_path):
    with open(config_yaml_path) as config_yaml_file:
        return yaml.safe_load(config_yaml_file)

def read_sc_params(supercomputer_conf, supercomputer):
    with open(supercomputer_conf) as sc_yaml_file:
        params = yaml.safe_load(sc_yaml_file)

        # Check mandatory parameters
        #   workflows_path: [mandatory] Path to the BioBB_HPC workflows
        #   biobb_path: [mandatory] Path to the BioBB/PyCOMPSs Conda Pack
        #   num_cores_node: [mandatory] Number of cores in each computing node

        if supercomputer not in params:
            sys.exit("Supercomputer name not found in supercomputer config file (sc_conf.yml). Please check it and try again.")

        if 'workflows_path' not in params[supercomputer]:
            sys.exit("Mandatory parameter workflows_path not found in supercomputer config file (sc_conf.yml). Please check it and try again.")

        if 'biobb_path' not in params[supercomputer]:
            sys.exit("Mandatory parameter biobb_path not found in supercomputer config file (sc_conf.yml). Please check it and try again.")

        if 'num_cores_node' not in params[supercomputer]:
            sys.exit("Mandatory parameter num_cores_node not found in supercomputer config file (sc_conf.yml). Please check it and try again.")

        return params[supercomputer]

def launch(input_structures, queue, num_nodes, compss_version, md_length,
           base_dir, compss_debug, time, task_timeout, output_dir, job_name, mpi_nodes,
           gmxlib, system, ff, concentration, box_type, box_size,
           supercomputer, supercomputer_conf, mpibin):

    params = read_sc_params(supercomputer_conf, supercomputer)

    # Check optional parameters
    #   modules: [optional] Computer-specific modules to be loaded
    #   mpi_env: [optional] MPI environment. Available options: SLURM / MPI / OMPI
    #   mpi_flags: [optional] MPI extra flags
    #   project_name: [optional] Account project name

    modules = ''
    if 'modules_load' in params:
        modules = params['modules_load']

    modules_unload = ''
    if 'modules_unload' in params:
        modules_unload = params['modules_unload']

    project_name = ''
    if 'project_name' in params:
        project_name = params['project_name']

    #queue = "default"
    if 'queue' in params:
        queue = params['queue']

    partition = "default"
    if 'partition' in params:
        partition = params['partition']

    mpi_env = "SLURM"
    if 'mpi_env' in params:
        mpi_env = params['mpi_env']

    mpi_flags = ''
    if 'mpi_flags' in params:
        mpi_flags = params['mpi_flags']

    extra_env = ()
    if 'extra_env' in params:
        extra_env = params['extra_env'].split(',')

    sc_cfg = supercomputer
    if 'cfg_path' in params:
        sc_cfg = params['cfg_path']

    base_dir = Path(params['workflows_path'])

    template_py_path = base_dir.joinpath('workflows', 'MD', 'md_list.py')
    template_yaml_path = base_dir.joinpath('workflows', 'MD', 'md_list.yaml')

    # Create working dir path
    work_dir = Path('.')

    currentDir = os.getcwd()
    working_dir_path = work_dir.joinpath(currentDir, 'MDs', 'md_set')

    if output_dir:
        if output_dir.startswith('/'):
            working_dir_path = Path(output_dir).resolve()
        else:
            working_dir_path = work_dir.joinpath(currentDir,output_dir)

    working_dir_path.mkdir(parents=True, exist_ok=True)

    # Check if it's the first launch
    run_number = 0
    run_dir = working_dir_path.joinpath("wf_md")
    config_yaml_path = working_dir_path.joinpath(f"md.yaml")
    wf_py_path = working_dir_path.joinpath(f"md.py")
    prolog_path = working_dir_path.joinpath(f"prolog.sh")
    launch_path = working_dir_path.joinpath(f"launch.sh")
    if not job_name:
        job_name = "mdlaunch_job"

    job_name_orig = job_name
    while run_dir.exists():
        run_number += 1
        run_dir = working_dir_path.joinpath(f"wf_md_{str(run_number)}")
        config_yaml_path = working_dir_path.joinpath(f"md_{str(run_number)}.yaml")
        wf_py_path = working_dir_path.joinpath(f"md_{str(run_number)}.py")
        prolog_path = working_dir_path.joinpath(f"prolog_{str(run_number)}.sh")
        launch_path = working_dir_path.joinpath(f"launch_{str(run_number)}.sh")
        job_name = f"{job_name_orig}_{str(run_number)}"

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

    # MPI binary
    config_dict['step7_mdrun_min']['properties']['mpi_bin'] = mpibin
    config_dict['step10_mdrun_nvt']['properties']['mpi_bin'] = mpibin
    config_dict['step12_mdrun_npt']['properties']['mpi_bin'] = mpibin
    config_dict['step14_mdrun_md']['properties']['mpi_bin'] = mpibin

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

    # Creating prolog (env_script) file
    with open(prolog_path, 'w') as prolog_file:
        prolog_file.write(f"#!/bin/bash\n\n")
        prolog_file.write(f"# BioBB + PyCOMPSs environment\n")
        prolog_file.write(f"source {params['biobb_path']}/bin/activate\n\n")
        if modules_unload:
            prolog_file.write(f"# Machine-specific modules environment (unload)\n")
            prolog_file.write(f"module unload {modules_unload}\n\n")
        if modules:
            prolog_file.write(f"# Machine-specific modules environment (load)\n")
            prolog_file.write(f"module load {modules}\n\n")
        prolog_file.write(f"# Multinode MPI environment\n")
        prolog_file.write(f"export TASK_COMPUTING_NODES={mpi_nodes}\n")
        prolog_file.write(f"export TASK_COMPUTING_UNITS={params['num_cores_node']}\n")
        prolog_file.write(f"export MULTINODE_MPI_ENV={mpi_env}\n")
        if mpi_flags:
            prolog_file.write(f"export MULTINODE_MPI_EXTRA_FLAGS=\"{mpi_flags}\"\n")
        prolog_file.write(f"\n")
        if extra_env:
            prolog_file.write(f"# Extra environment variables\n")
            for new_env in extra_env:
                prolog_file.write(f"export {new_env}\n")
        prolog_file.write(f"export TASK_TIME_OUT={task_timeout}\n")

    # Create launch
    with open(launch_path, 'w') as launch_file:
        launch_file.write(f"#!/bin/bash\n\n")
        launch_file.write(f"# BioBB + PyCOMPSs environment\n")
        launch_file.write(f"source {params['biobb_path']}/bin/activate\n\n")
        launch_file.write(f"enqueue_compss ")
        if compss_debug:
            launch_file.write(f"-d --keep_workingdir ")

        if num_nodes == 1 or num_nodes == mpi_nodes :
            launch_file.write(f"--worker_in_master_cpus={params['num_cores_node']} ")
        else:
            launch_file.write(f"--worker_in_master_cpus=0 ")

        if project_name:
            launch_file.write(f"--project_name={project_name} ")
        launch_file.write(f"--job_name={job_name}  --num_nodes={num_nodes} \
--exec_time={str(time)} --base_log_dir=$PWD --worker_working_dir=$PWD \
--master_working_dir=$PWD --network=ethernet --qos={queue}  \
--sc_cfg={sc_cfg}.cfg --queue={partition} \
--env_script={prolog_path} {wf_py_path} --config {config_yaml_path} ")
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
    parser.add_argument('-q', '--queue', required=False, default='default', type=str, help="(bsc_ls) [bsc_ls|debug]")
    parser.add_argument('-t', '--time', required=False, default=120, type=int, help="(120) [integer] Time in minutes")
    parser.add_argument('-tt', '--task_timeout', required=False, default=864000, type=int, help="(864000) [integer] Task timeout in seconds (default 10 days)")
    parser.add_argument('-nn', '--num_nodes', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-cv', '--compss_version', required=False, default='2.6.1', type=str, help="(2.6.1) [version_name]")
    parser.add_argument('-d', '--compss_debug', required=False, help="Compss debug mode", action='store_true')
    parser.add_argument('-l', '--md_length', required=False, default=10, type=int, help="(10ns) [integer] MD length in nanoseconds")
    parser.add_argument('-mpi', '--mpi_nodes', required=False, default=1, type=int, help="(1) [integer] Number of MPI nodes to be used per MD simulation")
    parser.add_argument('--base_dir', required=False, default='.', type=str, help="('.') [path_to_base_dir]")
    parser.add_argument('-o', '--output_dir', required=False, default='', type=str, help="Output dir name: If output_dir is absolute it will be respected if it's a relative path: /base_dir/MDs/runs/output_dir', if output_dir not exists, the name is autogenerated.")
    parser.add_argument('-jn', '--job_name', required=False, default='', type=str, help="Job name if it not exists, the name is autogenerated.")
    parser.add_argument('-gl', '--gromacs_lib', required=False, default='.', type=str, help="Gromacs lib, path where to find the force field libraries.")
    parser.add_argument('-sc', '--supercomputer', required=True, default='mn', type=str, help="Supercomputer name or id, included in the supercomputer-specific configuration file (sc_conf parameter).")
    parser.add_argument('-sc_conf', '--supercomputer_conf', required=True, default='sc_conf.yml', type=str, help="Supercomputer-specific parameters, such as MPI library or modules.")
    parser.add_argument('--mpi_bin', required=False, default='srun', type=str, help="MPI binary (e.g. srun, mpirun)")
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
           task_timeout=args.task_timeout,
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
           supercomputer=args.supercomputer,
           supercomputer_conf=args.supercomputer_conf,
           mpibin=args.mpi_bin,
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
