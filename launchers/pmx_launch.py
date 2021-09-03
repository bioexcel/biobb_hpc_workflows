#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import one_to_three
import re,os, sys
import shutil
import subprocess
import oyaml as yaml

def get_mutation_dict(mutation):
    if mutation.strip()[1].isdigit():
        pattern = re.compile(r"[A-Z]*:*(?P<wt>[a-zA-Z]{1})(?P<resnum>\d+)(?P<mt>[a-zA-Z]{1})")
        mut_dict = pattern.match(mutation.strip()).groupdict()
        mut_dict['wt'] = one_to_three(mut_dict['wt'].upper())
        mut_dict['mt'] = one_to_three(mut_dict['mt'].upper())
    else:
        pattern = re.compile(r"[A-Z]*:*(?P<wt>[a-zA-Z]{3})(?P<resnum>\d+)(?P<mt>[a-zA-Z]{3})")
        mut_dict = pattern.match(mutation.strip()).groupdict()

    return mut_dict

def forward_mutations(mutation,pmx_resnum):
    for_mut = []
    mut_list = mutation.split(",")
    for mut in mut_list:
        mut_dict = get_mutation_dict(mut)
        offset_pmx = int(mut_dict['resnum'])-pmx_resnum
#        mut_rev = f"{mut_dict['wt']}{str(offset_pmx)}{mut_dict['mt']}"
        mut_rev = f"{str(offset_pmx)}{mut_dict['mt']}"
        for_mut.append(mut_rev)
    return ','.join(for_mut)

def reverse_mutations(mutation,pmx_resnum):
    rev_mut = []
    mut_list = mutation.split(",")
    for mut in mut_list:
        mut_dict = get_mutation_dict(mut)
        offset_pmx = int(mut_dict['resnum'])-pmx_resnum
#        mut_rev = f"{mut_dict['mt']}{str(offset_pmx)}{mut_dict['wt']}"
        mut_rev = f"{str(offset_pmx)}{mut_dict['wt']}"
        rev_mut.append(mut_rev)
    return ','.join(rev_mut)

def three_to_one_mutation(mutation):
    mut_dict = get_mutation_dict(mutation)
    return f"{mut_dict.get('resnum')}{three_to_one(mut_dict.get('mt').upper())}"

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

def launch(mutation, pmx_resnum, wt_top, wt_trj, mut_top, mut_trj, queue, num_nodes, compss_version, ions, fe_nsteps, fe_length,
           base_dir, compss_debug, time, output_dir, fe_dt, num_frames, wt_trjconv_skip, mut_trjconv_skip,
           wt_start, wt_end, mut_start, mut_end, job_name,
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

    #if pmx_resnum == 0:
    #    pmx_resnum = int(get_mutation_dict(mutation)['resnum'])

    # Mutation code
    mutation_code = mutation.replace(":", "_").replace(",", "-")

    # GROMACS make_ndx group ids
    ndx1 = 20
    ndx2 = 21

    mut_list = mutation.split(",")
    nmuts = len(mut_list)
    #ions = 0

    ndx1 = ndx1 + (nmuts + ions)*2
    ndx2 = ndx2 + (nmuts + ions)*2

    base_dir = Path(base_dir)
    pth_path = Path.home().joinpath('.local', 'lib', 'python3.6', 'site-packages', 'biobb.pth')

    template_yaml_path = base_dir.joinpath('workflows', 'PMX', 'pmx_cv.yaml')
    template_py_path = base_dir.joinpath('workflows', 'PMX', 'pmx_cv_cufix_term.py')

    # Get input trajs
    traj_wt_tpr_path = wt_top
    traj_wt_xtc_path = wt_trj
    traj_mut_tpr_path = mut_top
    traj_mut_xtc_path = mut_trj
    print('WT trajs: ')
    print(traj_wt_tpr_path, traj_wt_xtc_path)
    print('\n')
    print(f'{mutation_code} trajs: ')
    print(traj_mut_tpr_path, traj_mut_xtc_path)
    print('\n')

    # Create working dir path

    long_name = f"{mutation_code}_{str(num_frames)}f_{str(fe_length)}ps"
    wt_start_str = str(wt_start)
    wt_end_str = str(wt_end)+'ps'
    if wt_start_str == '0' and wt_end_str == '0ps':
        wt_start_str = 'all'
        wt_end_str = 'all'
    mut_start_str = str(mut_start)
    mut_end_str = str(mut_end) + 'ps'
    if mut_start_str == '0' and mut_end_str == '0ps':
        mut_start_str = 'all'
        mut_end_str = 'all'

    work_dir = Path('.')

    currentDir = os.getcwd()
    working_dir_path = work_dir.joinpath(currentDir, 'PMX', long_name, f"{wt_start_str}to{wt_end_str}_{mut_start_str}to{mut_end_str}")

    if output_dir:
        if output_dir.startswith('/'):
            working_dir_path = Path(output_dir).resolve()
        else:
            working_dir_path = work_dir.joinpath(currentDir,output_dir)

    working_dir_path.mkdir(parents=True, exist_ok=True)

    # Check if it's the first launch
    run_number = 0
    run_dir = working_dir_path.joinpath("wf_pmx")
    config_yaml_path = working_dir_path.joinpath(f"pmx_biobb.yaml")
    wf_py_path = working_dir_path.joinpath(f"pmx_biobb.py")
    prolog_path = working_dir_path.joinpath(f"prolog_biobb.sh")
    launch_path = working_dir_path.joinpath(f"launch_biobb.sh")
    if not job_name:
        job_name = long_name

    job_name_orig = job_name
    while run_dir.exists():
        run_number += 1
        run_dir = working_dir_path.joinpath(f"wf_pmx_{str(run_number)}")
        config_yaml_path = working_dir_path.joinpath(f"pmx_biobb_{str(run_number)}.yaml")
        wf_py_path = working_dir_path.joinpath(f"pmx_biobb_{str(run_number)}.py")
        prolog_path = working_dir_path.joinpath(f"prolog_biobb_{str(run_number)}.sh")
        launch_path = working_dir_path.joinpath(f"launch_biobb_{str(run_number)}.sh")
        job_name = f"{job_name_orig}_{str(run_number)}"

    # Copy py file
    shutil.copyfile(template_py_path, wf_py_path)

    # Read yaml template file
    config_dict = get_template_config_dict(template_yaml_path)
    # Update config_dict
    config_dict['working_dir_path'] = str(run_dir)
    forward_mut = forward_mutations(mutation,pmx_resnum)
    reverse_mut = reverse_mutations(mutation,pmx_resnum)
    #config_dict['mutations']['stateA'] = f"{mutation_dict['wt']}{str(pmx_resnum)}{mutation_dict['mt']}"
    config_dict['mutations']['stateA'] = forward_mut
    config_dict['mutations']['stateB'] = reverse_mut
    #reverse_mut = reverse_mutations(mutation)
    #config_dict['mutations']['stateA'] = mutation
    #config_dict['mutations']['stateB'] = reverse_mut
    config_dict['input_trajs']['stateA']['input_tpr_path'] = str(traj_wt_tpr_path)
    config_dict['input_trajs']['stateA']['input_traj_path'] = str(traj_wt_xtc_path)
    config_dict['input_trajs']['stateB']['input_tpr_path'] = str(traj_mut_tpr_path)
    config_dict['input_trajs']['stateB']['input_traj_path'] = str(traj_mut_xtc_path)
    config_dict['step1_trjconv_stateA']['properties']['skip'] = wt_trjconv_skip
    config_dict['step1_trjconv_stateA']['properties']['start'] = wt_start
    config_dict['step1_trjconv_stateA']['properties']['end'] = wt_end
    config_dict['step1_trjconv_stateB']['properties']['skip'] = mut_trjconv_skip
    config_dict['step1_trjconv_stateB']['properties']['start'] = mut_start
    config_dict['step1_trjconv_stateB']['properties']['end'] = mut_end
    config_dict['step4_gmx_makendx']['properties']['selection'] = " a D*\\n0 & ! {}\\nname {} FREEZE".format(str(ndx1),str(ndx2))
    config_dict['step9_gmx_grompp']['properties']['mdp']['nsteps'] = fe_nsteps
    config_dict['step9_gmx_grompp']['properties']['mdp']['delta-lambda'] =  float(f'{1 / fe_nsteps:.0g}')

    with open(config_yaml_path, 'w') as config_yaml_file:
        config_yaml_file.write(yaml.dump(config_dict))

    # Creating prolog (env_script) file
    with open(prolog_path, 'w') as prolog_file:
        prolog_file.write(f"#!/bin/bash\n\n")
        prolog_file.write(f"# BioBB + PyCOMPSs environment\n")
        prolog_file.write(f"source {params['biobb_path']}/bin/activate\n\n")
        if modules:
            prolog_file.write(f"# Machine-specific modules environment (load)\n")
            prolog_file.write(f"module load {modules}\n\n")
        if modules_unload:
            prolog_file.write(f"# Machine-specific modules environment (unload)\n")
            prolog_file.write(f"module unload {modules_unload}\n\n")
        prolog_file.write(f"# Multinode MPI environment\n")
        prolog_file.write(f"export TASK_COMPUTING_NODES=2\n")
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
#        if num_nodes == 1 or num_nodes == mpi_nodes :
        if num_nodes == 1 :
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
    parser = argparse.ArgumentParser(description="Free Energy estimation upon amino acid modification with fast growth Thermodynamic Integration using GROMACS and pmx tools.")
    parser.add_argument('-m', '--mutation', required=True, help="Mutation in 'Leu858Arg' format")
    parser.add_argument('-prn', '--pmx_resnum', required=False, default=0, type=int, help="(0) [integer]")
    #parser.add_argument('-wr', '--wt_replica', required=False, default=1, type=int, help="(1) [integer]")
    #parser.add_argument('-mr', '--mut_replica', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-wt_top', '--wt_topology', required=True, default='wt.tpr', type=str, help="(wt.tpr) [Path to the WT topology]")
    parser.add_argument('-wt_trj', '--wt_trajectory', required=True, default='wt.xtc', type=str, help="(wt.xtc) [Path to the WT trajectory]")
    parser.add_argument('-mut_top', '--mut_topology', required=True, default='mut.tpr', type=str, help="(mut.tpr) [Path to the MUT topology]")
    parser.add_argument('-mut_trj', '--mut_trajectory', required=True, default='mut.xtc', type=str, help="(mut.xtc) [Path to the MUT trajectory]")
    parser.add_argument('-q', '--queue', required=False, default='bsc_ls', type=str, help="(bsc_ls) [bsc_ls|debug]")
    parser.add_argument('-t', '--time', required=False, default=120, type=int, help="(120) [integer] Time in minutes")
    parser.add_argument('-nn', '--num_nodes', required=False, default=1, type=int, help="(1) [integer]")
    parser.add_argument('-cv', '--compss_version', required=False, default='2.6.1', type=str, help="(2.6.1) [version_name]")
    parser.add_argument('-d', '--compss_debug', required=False, help="Compss debug mode", action='store_true')
    parser.add_argument('-i', '--ions', required=False, default=0, type=int, help="(0) [integer] Number of structural ions")
    parser.add_argument('-fe', '--fe_length', required=False, default=50, type=int, help="(50) [integer] Number of picoseconds")
    parser.add_argument('-nf', '--num_frames', required=False, default=100, type=int, help="(100) [integer] Number of frames to be extracted of trajectory")
    parser.add_argument('--mut_start_end_num_frames', required=False, default=10000, type=int, help="(10000) [integer] Total number of frames between start and end of the mutated trajectory")
    parser.add_argument('--wt_start_end_num_frames', required=False, default=10000, type=int, help="(10000) [integer] Total number of frames between start and end of the wt trajectory")
    parser.add_argument('--wt_start', required=False, default=0, type=int, help="(0) [integer] Time of first frame to read from WT trajectory (default unit ps).")
    parser.add_argument('--wt_end', required=False, default=0, type=int, help="(0) [integer] Time of last frame to read from WT trajectory (default unit ps).")
    parser.add_argument('--mut_start', required=False, default=0, type=int, help="(0) [integer] Time of first frame to read from MUT trajectory (default unit ps).")
    parser.add_argument('--mut_end', required=False, default=0, type=int, help="(0) [integer] Time of last frame to read from MUT trajectory (default unit ps).")
    parser.add_argument('--base_dir', required=False, default='/apps/BIOBB/workflows', type=str, help="('.') [path_to_base_dir]")
    parser.add_argument('-o', '--output_dir', required=False, default='', type=str, help="Output dir name: If output_dir is absolute it will be respected if it's a relative path: /base_dir/PMX/pmxlaunch/runs/output_dir', if output_dir not exists, the name is autogenerated.")
    parser.add_argument('-jn', '--job_name', required=False, default='', type=str, help="Job name if it not exists, the name is autogenerated.")
    parser.add_argument('--free_energy_dt', required=False, default=0.002, type=float, help="(0.002) [float] Integration time in picoseconds")
    parser.add_argument('-sc', '--supercomputer', required=True, default='mn', type=str, help="Supercomputer name or id, included in the supercomputer-specific configuration file (sc_conf parameter).")
    parser.add_argument('-sc_conf', '--supercomputer_conf', required=True, default='sc_conf.yml', type=str, help="Supercomputer-specific parameters, such as MPI library or modules.")
    args = parser.parse_args()

    # Specific call of each building block
    launch(mutation=args.mutation,
           pmx_resnum=args.pmx_resnum,
           wt_top=args.wt_topology,
           wt_trj=args.wt_trajectory,
           mut_top=args.mut_topology,
           mut_trj=args.mut_trajectory,
           queue=args.queue,
           time=args.time,
           num_nodes=args.num_nodes,
           compss_version=args.compss_version,
           ions=args.ions,
           compss_debug=args.compss_debug,
           fe_nsteps=int(args.fe_length/args.free_energy_dt),
           fe_length=args.fe_length,
           output_dir=args.output_dir,
           job_name=args.job_name,
           fe_dt=args.free_energy_dt,
           num_frames=args.num_frames,
           wt_trjconv_skip=args.wt_start_end_num_frames // args.num_frames,
           mut_trjconv_skip=args.mut_start_end_num_frames//args.num_frames,
           wt_start=args.wt_start,
           wt_end=args.wt_end,
           mut_start=args.mut_start,
           mut_end=args.mut_end,
           supercomputer=args.supercomputer,
           supercomputer_conf=args.supercomputer_conf,
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
