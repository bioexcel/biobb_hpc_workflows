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

def launch(mutation, wt_str, queue, num_nodes, compss_version, md_length, ff,
           base_dir, compss_debug, time, output_dir, job_name, mpi_nodes,
           cumulative, gmxlib,
           supercomputer, supercomputer_conf):

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

    queue = "default"
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

    base_dir = Path(params['workflows_path'])

    if cumulative == True :
        template_py_path = base_dir.joinpath('workflows', 'MD', 'md_add_muts_wt.py')
    else:
        template_py_path = base_dir.joinpath('workflows', 'MD', 'md_muts_sets.py')

    template_yaml_path = base_dir.joinpath('workflows', 'MD', 'md_muts_sets.yaml')
    #template_py_path = base_dir.joinpath('workflows', 'MD', 'md_add_muts_wt.py')

    # Create working dir path
    work_dir = Path('.')

    currentDir = os.getcwd()
    working_dir_path = work_dir.joinpath(currentDir, 'MDs', 'md_muts_set')

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
    config_dict['mutations'] = mutation
    config_dict['input_pdb'] = wt_str
    config_dict['step13_grompp_md']['properties']['mdp']['nsteps'] = int((md_length*1000)/0.002)
    config_dict['step2_pdb2gmx']['properties']['gmxlib'] = gmxlib
    config_dict['step5_grompp_genion']['properties']['gmxlib'] = gmxlib
    config_dict['step7_grompp_min']['properties']['gmxlib'] = gmxlib
    config_dict['step9_grompp_nvt']['properties']['gmxlib'] = gmxlib
    config_dict['step11_grompp_npt']['properties']['gmxlib'] = gmxlib
    config_dict['step13_grompp_md']['properties']['gmxlib'] = gmxlib

    # Force Field
    config_dict['step2_pdb2gmx']['properties']['force_field'] = ff

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
        prolog_file.write(f"export MULTINODE_MPI_ENV={params['mpi_env']}\n")
        if mpi_flags:
            prolog_file.write(f"export MULTINODE_MPI_EXTRA_FLAGS=\"{mpi_flags}\"\n")
        prolog_file.write(f"\n")
        if extra_env:
            prolog_file.write(f"# Extra environment variables\n")
            for new_env in extra_env:
                prolog_file.write(f"export {new_env}\n")

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
        if project_name:
            launch_file.write(f"--project_name={project_name} ")
        launch_file.write(f"--job_name={job_name}  --num_nodes={num_nodes} \
--exec_time={str(time)} --base_log_dir=$PWD --worker_working_dir=$PWD \
--master_working_dir=$PWD --network=ethernet --qos={queue}  \
--sc_cfg={supercomputer}.cfg --worker_in_master_cpus=0 --queue={partition} \
--env_script={prolog_path} {wf_py_path} --config {config_yaml_path} ")
        launch_file.write(f"\n")

    subprocess.call(f"bash {launch_path}", shell=True)


def main():
    parser = argparse.ArgumentParser(description="Workflow to model,setup and run MD simulations for a set of mutations.")
    parser.add_argument('-m', '--mutation', required=True, help="Mutations set in this format: WT+A:Arg6Gln,A:Asn13Lys+A:Glu31Asn,A:Lys43Asn where different MD simulations are separated by the '+' sign, whereas the set of mutations to be applied in each simulation are separated by the ',' sign. For example, this input: WT+A:Arg6Gln,A:Asn13Lys+A:Glu31Asn,A:Lys43Asn will generate 3 MD simulations, the first one without mutations (WT), the second one with two mutations (A:Arg6Gln,A:Asn13Lys) and the third one with two mutations (A:Glu31Asn,A:Lys43Asn). ")
    parser.add_argument('-wt', '--wt_structure', required=True, default='wt.pdb', type=str, help="(wt.pdb) [Path to the WT structure]")
    parser.add_argument('-q', '--queue', required=False, default='default', type=str, help="(bsc_ls) [bsc_ls|debug]")
    parser.add_argument('-f', '--force_field', required=False, default='amber99sb-ildn', type=str, help="(amber99sb-ildn) Force Field to be used in the simulation [amber99sb-ildn|charmm27|gromos54a7|oplsaa]")
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
    parser.add_argument('-sc', '--supercomputer', required=True, default='mn', type=str, help="Supercomputer name or id, included in the supercomputer-specific configuration file (sc_conf parameter).")
    parser.add_argument('-sc_conf', '--supercomputer_conf', required=True, default='sc_conf.yml', type=str, help="Supercomputer-specific parameters, such as MPI library or modules.")
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
           ff=args.force_field,
           supercomputer=args.supercomputer,
           supercomputer_conf=args.supercomputer_conf,
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
