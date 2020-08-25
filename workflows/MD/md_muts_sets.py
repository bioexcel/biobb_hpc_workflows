#!/usr/bin/env python3

import os
import sys
import zipfile
import time
import argparse
import subprocess
from pycompss.api.api import compss_wait_on_file

# biobb common modules
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

# pycompss: biobb model modules
from biobb_adapters.pycompss.biobb_model.model.mutate_pc import mutate_pc

# pycompss: biobb md modules
from biobb_adapters.pycompss.biobb_md.gromacs.pdb2gmx_pc import pdb2gmx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.editconf_pc import editconf_pc
from biobb_adapters.pycompss.biobb_md.gromacs.solvate_pc import solvate_pc
from biobb_adapters.pycompss.biobb_md.gromacs.genion_pc import genion_pc
from biobb_adapters.pycompss.biobb_md.gromacs.make_ndx_pc import make_ndx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.grompp_pc import grompp_pc
from biobb_adapters.pycompss.biobb_md.gromacs.grompp_ndx_pc import grompp_ndx_pc
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun_pc import mdrun_pc
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun_dhdl_pc import mdrun_dhdl_pc

def main(config, system=None):
    start_time = time.time()
    conf = settings.ConfReader(config, system)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic()
    global_paths = conf.get_paths_dic()

    mutations = []
    for mut_set in conf.properties['mutations'].split('+'):
        mutations.append([mut_set])

    #for mutation_list in conf.properties['mutations']:
    for mutation_list in mutations:


        if "WT" in mutation_list:

            mut_prop = conf.get_prop_dic(prefix="WT")
            mut_paths = conf.get_paths_dic(prefix="WT")

            mutpdb = conf.properties['input_pdb']

        else:

    	    # Converting list to String
            sep = ","
            str_mut_list = sep.join(mutation_list)

            # Updated cumulative list of mutations
            global_log.info("Mutation list: " + str_mut_list)

            mutation_code = str_mut_list.replace(":", "_").replace(",", "-")

            mut_prop = conf.get_prop_dic(prefix=mutation_code)
            mut_paths = conf.get_paths_dic(prefix=mutation_code)

            mut_paths['step1_mutate']['input_pdb_path'] = conf.properties['input_pdb']
            mut_prop['step1_mutate']['mutation_list'] = str_mut_list

            global_log.info("step1_mutate: Modeling a particular residue mutation")
            mutate_pc(**mut_paths["step1_mutate"], properties=mut_prop["step1_mutate"])
            compss_wait_on_file(mut_paths["step1_mutate"]["output_pdb_path"])

            mutpdb = mut_paths['step1_mutate']['output_pdb_path']

        # Fixing Histidine residues to HIE
        # Must be done after mutation because biopython does not understand HIE/D/P residues

        #pdbhie = mut_paths['step1_mutate']['output_pdb_path'] + ".hie.pdb"
        pdbhie = mutpdb + ".hie.pdb"

        cmd = "sed 's/HIS/HIE/g' " + mutpdb + " > " + pdbhie
        subprocess.call(cmd, shell=True)

        mut_paths['step2_pdb2gmx']['input_pdb_path'] = pdbhie

        global_log.info("step2_pdb2gmx: Generate the topology")
        pdb2gmx_pc(**mut_paths["step2_pdb2gmx"], properties=mut_prop["step2_pdb2gmx"])

        global_log.info("step3_editconf: Create the solvent box")
        editconf_pc(**mut_paths["step3_editconf"], properties=mut_prop["step3_editconf"])

        global_log.info("step4_solvate: Fill the solvent box with water molecules")
        solvate_pc(**mut_paths["step4_solvate"], properties=mut_prop["step4_solvate"])

        global_log.info("step5_grompp_genion: Preprocess ion generation")
        grompp_pc(**mut_paths["step5_grompp_genion"], properties=mut_prop["step5_grompp_genion"])

        global_log.info("step6_genion: Ion generation")
        genion_pc(**mut_paths["step6_genion"], properties=mut_prop["step6_genion"])

        global_log.info("step7_grompp_min: Preprocess energy minimization")
        grompp_pc(**mut_paths["step7_grompp_min"], properties=mut_prop["step7_grompp_min"])

        global_log.info("step8_mdrun_min: Execute energy minimization")
        mdrun_pc(**mut_paths["step8_mdrun_min"], properties=mut_prop["step8_mdrun_min"])

        global_log.info("step9_grompp_nvt: Preprocess system temperature equilibration")
        grompp_pc(**mut_paths["step9_grompp_nvt"], properties=mut_prop["step9_grompp_nvt"])

        global_log.info("step10_mdrun_nvt: Execute system temperature equilibration")
        mdrun_pc(**mut_paths["step10_mdrun_nvt"], properties=mut_prop["step10_mdrun_nvt"])

        global_log.info("step11_grompp_npt: Preprocess system pressure equilibration")
        grompp_pc(**mut_paths["step11_grompp_npt"], properties=mut_prop["step11_grompp_npt"])

        global_log.info("step12_mdrun_npt: Execute system pressure equilibration")
        mdrun_pc(**mut_paths["step12_mdrun_npt"], properties=mut_prop["step12_mdrun_npt"])

        global_log.info("step13_grompp_md: Preprocess free dynamics")
        grompp_pc(**mut_paths["step13_grompp_md"], properties=mut_prop["step13_grompp_md"])

        global_log.info("step14_mdrun_md: Execute free molecular dynamics simulation")
        mdrun_pc(**mut_paths["step14_mdrun_md"], properties=mut_prop["step14_mdrun_md"])

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
    parser = argparse.ArgumentParser(description="Cumulative mutation MDs")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)
