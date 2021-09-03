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

# pycompss: biobb md modules
from biobb_adapters.pycompss.biobb_md.gromacs.pdb2gmx import pdb2gmx
from biobb_adapters.pycompss.biobb_md.gromacs.editconf import editconf
from biobb_adapters.pycompss.biobb_md.gromacs.solvate import solvate
from biobb_adapters.pycompss.biobb_md.gromacs.genion import genion
from biobb_adapters.pycompss.biobb_md.gromacs.make_ndx import make_ndx
from biobb_adapters.pycompss.biobb_md.gromacs.grompp import grompp
from biobb_adapters.pycompss.biobb_md.gromacs.mdrun import mdrun

def main(config, system=None):
    start_time = time.time()
    conf = settings.ConfReader(config, system)
    global_log, _ = fu.get_logs(path=conf.get_working_dir_path(), light_format=True)
    global_prop = conf.get_prop_dic()
    global_paths = conf.get_paths_dic()

    for structure in conf.properties['input_structures'].split(','):

        prefix_str = os.path.basename(structure)
        prefix_str = prefix_str.replace('.','_')

        mut_prop = conf.get_prop_dic(prefix=prefix_str)
        mut_paths = conf.get_paths_dic(prefix=prefix_str)

        global_log.info("Starting setup process for PDB: " + prefix_str)
        mut_paths['step1_pdb2gmx']['input_pdb_path'] = structure

        global_log.info("step1_pdb2gmx: Generate the topology")
        pdb2gmx(**mut_paths["step1_pdb2gmx"], properties=mut_prop["step1_pdb2gmx"])

        global_log.info("step2_editconf: Create the solvent box")
        editconf(**mut_paths["step2_editconf"], properties=mut_prop["step2_editconf"])

        global_log.info("step3_solvate: Fill the solvent box with water molecules")
        solvate(**mut_paths["step3_solvate"], properties=mut_prop["step3_solvate"])

        global_log.info("step4_grompp_genion: Preprocess ion generation")
        grompp(**mut_paths["step4_grompp_genion"], properties=mut_prop["step4_grompp_genion"])

        global_log.info("step5_genion: Ion generation")
        genion(**mut_paths["step5_genion"], properties=mut_prop["step5_genion"])

        global_log.info("step6_grompp_min: Preprocess energy minimization")
        grompp(**mut_paths["step6_grompp_min"], properties=mut_prop["step6_grompp_min"])

        global_log.info("step7_mdrun_min: Execute energy minimization")
        mdrun(**mut_paths["step7_mdrun_min"], properties=mut_prop["step7_mdrun_min"])

        global_log.info("step8_make_ndx: Generate GROMACS index file")
        make_ndx(**mut_paths["step8_make_ndx"], properties=mut_prop["step8_make_ndx"])

        global_log.info("step9_grompp_nvt: Preprocess system temperature equilibration")
        grompp(**mut_paths["step9_grompp_nvt"], properties=mut_prop["step9_grompp_nvt"])

        global_log.info("step10_mdrun_nvt: Execute system temperature equilibration")
        mdrun(**mut_paths["step10_mdrun_nvt"], properties=mut_prop["step10_mdrun_nvt"])

        global_log.info("step11_grompp_npt: Preprocess system pressure equilibration")
        grompp(**mut_paths["step11_grompp_npt"], properties=mut_prop["step11_grompp_npt"])

        global_log.info("step12_mdrun_npt: Execute system pressure equilibration")
        mdrun(**mut_paths["step12_mdrun_npt"], properties=mut_prop["step12_mdrun_npt"])

        global_log.info("step13_grompp_md: Preprocess free dynamics")
        grompp(**mut_paths["step13_grompp_md"], properties=mut_prop["step13_grompp_md"])

        global_log.info("step14_mdrun_md: Execute free molecular dynamics simulation")
        mdrun(**mut_paths["step14_mdrun_md"], properties=mut_prop["step14_mdrun_md"])

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
    parser = argparse.ArgumentParser(description="MD setup and run for a list of structures.")
    parser.add_argument('--config', required=True)
    parser.add_argument('--system', required=False)
    args = parser.parse_args()
    main(args.config, args.system)
