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


def launch(mutation, wt_str, queue, num_nodes, compss_version, md_length,
           base_dir, compss_debug, time, output_dir, job_name, mpi_nodes, cumulative, gmxlib):

    base_dir = Path(base_dir)
    #pth_path = Path.home().joinpath('.local', 'lib', 'python3.6', 'site-packages', 'biobb.pth')

    if cumulative == True :
        template_py_path = base_dir.joinpath('workflows', 'MD', 'md_add_muts_wt.py')
    else:
        template_py_path = base_dir.joinpath('workflows', 'MD', 'md_muts_sets.py')

    template_yaml_path = base_dir.joinpath('workflows', 'MD', 'md_muts_sets.yaml')
    #template_py_path = base_dir.joinpath('workflows', 'MD', 'md_add_muts_wt.py')

    # Create working dir path
    work_dir = Path('.')
    working_dir_path = work_dir.joinpath('MDs', 'md_muts_set')

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
    config_dict['mutations'] = mutation 
    config_dict['input_pdb'] = wt_str
    config_dict['step13_grompp_md']['properties']['mdp']['nsteps'] = int((md_length*1000)/0.002)
    config_dict['step2_pdb2gmx']['properties']['gmxlib'] = gmxlib
    config_dict['step5_grompp_genion']['properties']['gmxlib'] = gmxlib
    config_dict['step7_grompp_min']['properties']['gmxlib'] = gmxlib
    config_dict['step9_grompp_nvt']['properties']['gmxlib'] = gmxlib
    config_dict['step11_grompp_npt']['properties']['gmxlib'] = gmxlib
    config_dict['step13_grompp_md']['properties']['gmxlib'] = gmxlib

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
    parser = argparse.ArgumentParser(description="Workflow to model,setup and run MD simulations for a set of mutations.")
    parser.add_argument('-m', '--mutation', required=True, help="Mutations set in WT;A:Arg6Gln,A:Asn13Lys;A:Glu31Asn,A:Lys43Asn format")
    parser.add_argument('-wt', '--wt_structure', required=True, default='wt.pdb', type=str, help="(wt.pdb) [Path to the WT structure]")
    parser.add_argument('-q', '--queue', required=False, default='bsc_ls', type=str, help="(bsc_ls) [bsc_ls|debug]")
    parser.add_argument('-t', '--time', required=False, default=120, type=int, help="(120) [integer] Time in minutes")
    parser.add_argument('-nn', '--num_nodes', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-cv', '--compss_version', required=False, default='2.6.1', type=str, help="(2.6.1) [version_name]")
    parser.add_argument('-d', '--compss_debug', required=False, help="Compss debug mode", action='store_true')
    parser.add_argument('-l', '--md_length', required=False, default=10, type=int, help="(10ns) [integer] MD length in nanoseconds")
    parser.add_argument('-c', '--cumulative', required=False, default='False', type=str, help="(False) [Boolean] Mutations will be accumulated")
    parser.add_argument('-mpi', '--mpi_nodes', required=False, default=1, type=int, help="(1) [integer] Number of MPI nodes to be used per MD simulation")
    parser.add_argument('--base_dir', required=False, default='/apps/BIOBB/workflows', type=str, help="('.') [path_to_base_dir]")
    parser.add_argument('-o', '--output_dir', required=False, default='', type=str, help="Output dir name: If output_dir is absolute it will be respected if it's a relative path: /base_dir/MDs/runs/output_dir', if output_dir not exists, the name is autogenerated.")
    parser.add_argument('-jn', '--job_name', required=False, default='', type=str, help="Job name if it not exists, the name is autogenerated.")
    parser.add_argument('-gl', '--gromacs_lib', required=False, default='.', type=str, help="Gromacs lib, path where to find the force field libraries.")
    args = parser.parse_args()

    # Specific call of each building block
    launch(mutation=args.mutation,
           wt_str=args.wt_structure,
           queue=args.queue,
           time=args.time,
           num_nodes=args.num_nodes,
           compss_version=args.compss_version,
           compss_debug=args.compss_debug,
           md_length=args.md_length,
           cumulative=args.cumulative,
           mpi_nodes=args.mpi_nodes,
           output_dir=args.output_dir,
           job_name=args.job_name,
           gmxlib=args.gromacs_lib,
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
