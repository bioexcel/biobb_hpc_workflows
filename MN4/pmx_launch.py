#!/usr/bin/env python3

from pathlib import Path
import argparse
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import one_to_three
import re
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
        mut_rev = f"{mut_dict['wt']}{str(offset_pmx)}{mut_dict['mt']}"
        for_mut.append(mut_rev)
    return ','.join(for_mut)

def reverse_mutations(mutation,pmx_resnum):
    rev_mut = []
    mut_list = mutation.split(",")
    for mut in mut_list:
        mut_dict = get_mutation_dict(mut)
        offset_pmx = int(mut_dict['resnum'])-pmx_resnum
        mut_rev = f"{mut_dict['mt']}{str(offset_pmx)}{mut_dict['wt']}"
        rev_mut.append(mut_rev)
    return ','.join(rev_mut)

def three_to_one_mutation(mutation):
    mut_dict = get_mutation_dict(mutation)
    return f"{three_to_one(mut_dict.get('wt').upper())}{mut_dict.get('resnum')}{three_to_one(mut_dict.get('mt').upper())}"

def get_template_config_dict(config_yaml_path):
    with open(config_yaml_path) as config_yaml_file:
        return yaml.safe_load(config_yaml_file)

def launch(mutation, pmx_resnum, wt_top, wt_trj, mut_top, mut_trj, queue, num_nodes, compss_version, fe_nsteps, fe_length,
           base_dir, compss_debug, time, output_dir, fe_dt, num_frames, wt_trjconv_skip, mut_trjconv_skip,
           wt_start, wt_end, mut_start, mut_end, job_name):

    #if pmx_resnum == 0:
    #    pmx_resnum = int(get_mutation_dict(mutation)['resnum'])

    # Mutation code
    mutation_code = mutation.replace(":", "_").replace(",", "-")

    # GROMACS make_ndx group ids
    ndx1 = 20
    ndx2 = 21

    mut_list = mutation.split(",")
    nmuts = len(mut_list)
    ions = 0

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
    work_dir = Path('.')

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

    working_dir_path = work_dir.joinpath('PMX', long_name, f"{wt_start_str}to{wt_end_str}_{mut_start_str}to{mut_end_str}")

    if output_dir:
        if output_dir.startswith('/'):
            working_dir_path = Path(output_dir).resolve()
        else:
            working_dir_path = work_dir.joinpath('PMX', output_dir)
    working_dir_path.mkdir(parents=True, exist_ok=True)

    # Check if it's the first launch
    run_number = 0
    run_dir = working_dir_path.joinpath("wf_pmx")
    config_yaml_path = working_dir_path.joinpath(f"pmx.yaml")
    wf_py_path = working_dir_path.joinpath(f"pmx.py")
    launch_path = working_dir_path.joinpath(f"launch.sh")
    if not job_name:
        job_name = long_name
    while run_dir.exists():
        run_number += 1
        run_dir = working_dir_path.joinpath(f"wf_pmx_{str(run_number)}")
        config_yaml_path = working_dir_path.joinpath(f"pmx_{str(run_number)}.yaml")
        wf_py_path = working_dir_path.joinpath(f"pmx_{str(run_number)}.py")
        launch_path = working_dir_path.joinpath(f"launch_{str(run_number)}.sh")
        job_name = f"{job_name}_{str(run_number)}"

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

    # Create launch
    with open(launch_path, 'w') as launch_file:
        launch_file.write(f"#!/bin/bash\n")
        launch_file.write(f"\n")
        launch_file.write(f"module purge\n")
        launch_file.write(f"\n")
        launch_file.write(f"module load ANACONDA/2019.10\n")
        #launch_file.write(f"source activate biobb\n")
        launch_file.write(f"module load biobb/covid_dev\n")
        #launch_file.write(f"module load ANACONDA/2018.12_py3\n")
        #launch_file.write(f"source activate biobb\n")
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
        launch_file.write(f"\n")
        launch_file.write(f"#TASK_TIME_OUT env var\n")
        launch_file.write(f"#3600 seconds = 1h\n")
        launch_file.write(f"export TASK_TIME_OUT=\"3600\"\n")
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
                          --network=ethernet \
                          --qos={queue}  \
                          {wf_py_path} \
                          --config {config_yaml_path} ")
        #if imaged_traj_available:
        #    launch_file.write(f"--imaged_traj_available ")
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
           base_dir=Path(args.base_dir)
           )


if __name__ == '__main__':
    main()
