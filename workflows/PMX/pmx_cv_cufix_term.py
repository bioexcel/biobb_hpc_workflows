#!/usr/bin/env python3

import os
import io
import re
import sys
import zipfile
import time
import argparse
import subprocess
from shutil import copyfile
from pycompss.api.api import compss_wait_on_file
from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.multinode import multinode
from pycompss.api.parameter import FILE_IN, FILE_OUT

# biobb common modules
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

# pycompss: biobb pmx modules
from biobb_adapters.pycompss.biobb_pmx.pmx.pmxmutate import pmxmutate
from biobb_adapters.pycompss.biobb_pmx.pmx.pmxgentop import pmxgentop
from biobb_adapters.pycompss.biobb_pmx.pmx.pmxanalyse import pmxanalyse

# pycompss: biobb md modules
from biobb_adapters.pycompss.biobb_md.gromacs.pdb2gmx import pdb2gmx
from biobb_adapters.pycompss.biobb_md.gromacs.make_ndx import make_ndx
from biobb_adapters.pycompss.biobb_md.gromacs.grompp import grompp
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun import mdrun

# pycompss: biobb analysis modules
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_image import gmximage
from biobb_adapters.pycompss.biobb_analysis.gromacs.gmx_trjconv_str_ens import gmxtrjconvstrens

# pycompss: biobb structure utils modules
from biobb_adapters.pycompss.biobb_structure_utils.utils.extract_atoms import extractatoms

@constraint(computing_units="48")
@multinode(computing_nodes="1")
@task(output_structure_path = FILE_IN, input_structure_path = FILE_IN, output_ndx_path= FILE_OUT, input_gro_path=FILE_IN, output_tpr_path=FILE_OUT, input_top_zip_path=FILE_IN, output_trr_path=FILE_OUT, output_gro_path=FILE_OUT, output_edr_path=FILE_OUT, output_xtc_path=FILE_OUT, output_log_path=FILE_OUT, on_failure="IGNORE")
def check_structure_and_run_ndx(ensemble, output_structure_path, input_structure_path, output_ndx_path, properties_makendx, input_gro_path, input_top_zip_path, output_tpr_path, properties_grompp, output_trr_path, output_gro_path, output_edr_path, output_xtc_path, output_log_path, properties_mdrun, **kwargs):
    try:
        dummy = bool(os.path.getsize(output_structure_path))
        if not dummy:
        #if ensemble == "stateA":
            from shutil import copyfile
            copyfile(input_gro_path, output_gro_path)
            fu.write_failed_output(output_ndx_path)
            fu.write_failed_output(output_tpr_path)
            fu.write_failed_output(output_trr_path)
            fu.write_failed_output(output_edr_path)
            fu.write_failed_output(output_xtc_path)
            fu.write_failed_output(output_log_path)
        else:
            # step4_gmx_makendx
            make_ndx(input_structure_path=input_structure_path, output_ndx_path=output_ndx_path, properties=properties_makendx)
            # step5_gmx_grompp
            grompp(input_gro_path = input_gro_path, input_top_zip_path = input_top_zip_path, output_tpr_path=output_tpr_path, input_ndx_path=output_ndx_path, properties=properties_grompp)
            # step6_gmx_mdrun
            mdrun(input_tpr_path=output_tpr_path, output_trr_path= output_trr_path, output_gro_path = output_gro_path, output_edr_path=output_edr_path, output_xtc_path = output_xtc_path, output_log_path = output_log_path, properties=properties_mdrun)

    except Exception as e:
       raise Exception("Error ocurred while checking the file containing dummy atoms")

#@task(script=FILE_IN, input_gro=FILE_IN, mod_gro=FILE_OUT, on_failure="IGNORE")
#def fix_gro(script, input_gro, mod_gro):
@task(input_gro=FILE_IN, mod_gro=FILE_OUT, on_failure="IGNORE")
def fix_term(input_gro, mod_gro):

    copyfile(input_gro, mod_gro)

    cmd_grep_nter = "grep H3 " + input_gro + " | cut -c1-8"
    cmd_grep_cter = "grep OC2 " + input_gro + " | cut -c1-8"
    nter = subprocess.Popen(cmd_grep_nter, stdout=subprocess.PIPE, shell=True)
    cter = subprocess.Popen(cmd_grep_cter, stdout=subprocess.PIPE, shell=True)

    nter_residues = nter.communicate()[0].decode('utf-8')
    for nter_old in nter_residues.splitlines():
        r = re.compile(" +([0-9]+)([a-zA-Z]+)")
        m = r.match(nter_old)
        if m:
            nter_new = 'N'.join(m.groups())
            len_nter_new = str(len(nter_old)+1)
            nter_new_formatted = ("{:>"+len_nter_new+"}").format(nter_new)
            cmd = "sed -i 's/"+nter_old+" /"+nter_new_formatted+"/g' " + mod_gro
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            process.wait()
            if process.returncode != 0 :
                print("--  Process output  --")
                print(out.decode("utf-8"))
                print("--  Process error  --")
                print(err.decode("utf-8"))
                raise Exception(" Process return code different from 0")

    cter_residues = cter.communicate()[0].decode('utf-8')
    for cter_old in cter_residues.splitlines():
        r = re.compile(" +([0-9]+)([a-zA-Z]+)")
        m = r.match(cter_old)
        if m:
            cter_new = 'C'.join(m.groups())
            len_cter_new = str(len(nter_old)+1)
            cter_new_formatted = ("{:>"+len_cter_new+"}").format(cter_new)
            cmd = "sed -i 's/"+cter_old+" /"+cter_new_formatted+"/g' " + mod_gro
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            out, err = process.communicate()
            process.wait()
            if process.returncode != 0 :
                print("--  Process output  --")
                print(out.decode("utf-8"))
                print("--  Process error  --")
                print(err.decode("utf-8"))
                raise Exception(" Process return code different from 0")

def main(config, system=None):
    start_time = time.time()
    conf = settings.ConfReader(config, system)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic()
    global_paths = conf.get_paths_dic()

    dhdl_paths_listA = []
    dhdl_paths_listB = []
    for ensemble, mutation in conf.properties['mutations'].items():

        ensemble_prop = conf.get_prop_dic(prefix=ensemble)
        ensemble_paths = conf.get_paths_dic(prefix=ensemble)

        # step0_image
        global_log.info(ensemble+" Step 0: gmx image: Imaging trajectories to remove PBC issues")
        ensemble_paths['step0_image']['input_traj_path'] = conf.properties['input_trajs'][ensemble]['input_traj_path']
        ensemble_paths['step0_image']['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        gmximage(**ensemble_paths["step0_image"], properties=ensemble_prop["step0_image"])

        # step1_trjconv
        global_log.info(ensemble+" Step 1: gmx trjconv: Extract snapshots from equilibrium trajectories")
        ensemble_paths['step1_trjconv_'+ensemble]['input_top_path'] = conf.properties['input_trajs'][ensemble]['input_tpr_path']
        gmxtrjconvstrens(**ensemble_paths['step1_trjconv_'+ensemble], properties=ensemble_prop['step1_trjconv_'+ensemble])

    for ensemble, mutation in conf.properties['mutations'].items():
        ensemble_prop = conf.get_prop_dic(prefix=ensemble)
        ensemble_paths = conf.get_paths_dic(prefix=ensemble)
        compss_wait_on_file(ensemble_paths['step1_trjconv_'+ensemble]["output_str_ens_path"])

        with zipfile.ZipFile(ensemble_paths['step1_trjconv_'+ensemble]["output_str_ens_path"]) as zip_f:
            unique_dir = os.path.abspath(fu.create_unique_dir(prefix=ensemble_prop['step1_trjconv_'+ensemble]['working_dir_path']+'/'+ensemble+'/'))
            zip_f.extractall(unique_dir)
            state_pdb_list = [os.path.join(unique_dir, name)for name in zip_f.namelist()]

        for pdb_path in state_pdb_list:
            pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
            prop = conf.get_prop_dic(prefix=os.path.join(ensemble, pdb_name))
            paths = conf.get_paths_dic(prefix=os.path.join(ensemble, pdb_name))

            # Fixing terminal residues
            newpdb = pdb_path + ".term.gro"
            fix_term(pdb_path, newpdb)

            paths['step1_pmx_mutate']['input_structure_path'] = newpdb

            # step1_pmx_mutate
            global_log.info(ensemble+" "+pdb_name+" Step 1: pmx mutate: Generate Hybrid Structure")
            prop['step1_pmx_mutate']['mutation_list'] = mutation
            pmxmutate(**paths["step1_pmx_mutate"], properties=prop["step1_pmx_mutate"])

            # step1.1_check_dummies
            global_log.info(ensemble + " " + pdb_name + " Step 1.1 Check for dummy atoms")
            extractatoms(**paths['step1.1_check_dummies'], properties=prop['step1.1_check_dummies'])
            #compss_wait_on_file(paths['step1.1_check_dummies']['output_structure_path'])
            #try:
            #    dummy = bool(os.path.getsize(paths['step1.1_check_dummies']['output_structure_path']))
            #except:
            #    global_log.info("Error ocurred while checking the file containing dummy atoms")
            #    global_log.info(sys.exc_info()[0])

            # step1.2_remove_ligand
            #global_log.info(ensemble + " " + pdb_name + " Step 1.2 Remove ligand")
            #remove_ligand(**paths['step1.2_remove_ligand'], properties=prop['step1.2_remove_ligand'])

            # step2_gmx_pdb2gmx
            global_log.info(ensemble + " " + pdb_name + " Step 2: gmx pdb2gmx: Generate Topology")
            pdb2gmx(**paths["step2_gmx_pdb2gmx"], properties=prop["step2_gmx_pdb2gmx"])

            # step2.1_sort_gro
            #global_log.info(ensemble + " " + pdb_name + " Step 2.1 Sort gro residues")
            #sort_gro_residues(**paths['step2.1_sort_gro'], properties=prop['step2.1_sort_gro'])

            # step2.2_lig_gmx_appendLigand
            #global_log.info(ensemble + " " + pdb_name +" Step 2.2_lig: gmx appendLigand: Append a ligand to a GROMACS topology")
            #append_ligand(**paths["step2.2_lig_gmx_appendLigand"], properties=prop["step2.2_lig_gmx_appendLigand"])

            # step3_pmx_gentop
            global_log.info(ensemble + " " + pdb_name +" Step 3: pmx gentop: Generate Hybrid Topology")
            pmxgentop(**paths["step3_pmx_gentop"], properties=prop["step3_pmx_gentop"])

            # step4_gmx_makendx, step5_gmx_grompp, step6_gmx_mdrun
            #global_log.info(ensemble + " " + pdb_name +" Check structure and Step 4-6 (Dummies): gmx make_ndx: Generate Gromacs Index file to select atoms to freeze, gmx grompp: Creating portable binary run file for energy minimization, gmx mdrun: Running energy minimization")
            #check_structure_and_run_ndx(ensemble, **paths["step4_gmx_makendx"], properties_makendx=prop["step4_gmx_makendx"], **paths["step5_gmx_grompp"], properties_grompp=prop["step5_gmx_grompp"], **paths["step6_gmx_mdrun"], properties_mdrun=prop["step6_gmx_mdrun"])

            # step4_gmx_makendx, step5_gmx_grompp, step6_gmx_mdrun
            global_log.info(ensemble + " " + pdb_name +" Step 4-6 (Check structure and ndx): gmx make_ndx: Generate Gromacs Index file to select atoms to freeze, gmx grompp: Creating portable binary run file for energy minimization, gmx mdrun: Running energy minimization")
            fu.create_dir(prop["step4_gmx_makendx"]['working_dir_path'] + "/" + ensemble + "/" + pdb_name + "/step4_gmx_makendx")
            fu.create_dir(prop["step5_gmx_grompp"]['working_dir_path'] + "/" + ensemble + "/" + pdb_name + "/step5_gmx_grompp")
            fu.create_dir(prop["step6_gmx_mdrun"]['working_dir_path'] + "/" + ensemble + "/" + pdb_name + "/step6_gmx_mdrun")
            check_structure_and_run_ndx(ensemble, output_structure_path=paths['step1.1_check_dummies']['output_structure_path'], **paths["step4_gmx_makendx"], properties_makendx=prop["step4_gmx_makendx"], **paths["step5_gmx_grompp"], properties_grompp=prop["step5_gmx_grompp"], **paths["step6_gmx_mdrun"], properties_mdrun=prop["step6_gmx_mdrun"])

            # step7_gmx_grompp
            global_log.info(ensemble + " " + pdb_name +" Step 7: gmx grompp: Creating portable binary run file for system equilibration")
            grompp(**paths["step7_gmx_grompp"], properties=prop["step7_gmx_grompp"])

            # step8_gmx_mdrun
            global_log.info(ensemble + " " + pdb_name +" Step 8: gmx mdrun: Running system equilibration")
            mdrun(**paths["step8_gmx_mdrun"], properties=prop["step8_gmx_mdrun"])

            # step9_gmx_grompp
            global_log.info(ensemble + " " + pdb_name +" Step 9: Creating portable binary run file for thermodynamic integration (ti)")
            grompp(**paths["step9_gmx_grompp"], properties=prop["step9_gmx_grompp"])

            # step10_gmx_mdrun
            global_log.info(ensemble + " " + pdb_name +" Step 10: gmx mdrun: Running thermodynamic integration")
            mdrun(**paths["step10_gmx_mdrun"], properties=prop["step10_gmx_mdrun"])

            if ensemble == "stateA":
                dhdl_paths_listA.append(paths["step10_gmx_mdrun"]["output_dhdl_path"])
            elif ensemble == "stateB":
                dhdl_paths_listB.append(paths["step10_gmx_mdrun"]["output_dhdl_path"])

    # Creating zip file containing all the dhdl files
    dhdlA_path = os.path.join(global_prop["step11_pmx_analyse"]['working_dir_path'], 'dhdlA.zip')
    dhdlB_path = os.path.join(global_prop["step11_pmx_analyse"]['working_dir_path'], 'dhdlB.zip')
    for dhdl_file in dhdl_paths_listA:
        compss_wait_on_file(dhdl_file)
    for dhdl_file in dhdl_paths_listB:
        compss_wait_on_file(dhdl_file)
    fu.zip_list(dhdlA_path, dhdl_paths_listA)
    fu.zip_list(dhdlB_path, dhdl_paths_listB)

    # step11_pmx_analyse
    global_log.info(ensemble+" Step 11: pmx analyse: Calculate free energies from fast growth thermodynamic integration simulations")
    global_paths["step11_pmx_analyse"]["input_a_xvg_zip_path"]=dhdlA_path
    global_paths["step11_pmx_analyse"]["input_b_xvg_zip_path"]=dhdlB_path
    pmxanalyse(**global_paths["step11_pmx_analyse"], properties=global_prop["step11_pmx_analyse"])

    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % conf.get_working_dir_path())
    global_log.info('  Config File: %s' % config)
    if system:
        global_log.info('  System: %s' % system)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Based on the official PMX tutorial")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)
