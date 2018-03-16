"""
Module that calls the necessary functions and runs the re-scoring algorithm.
"""

import argparse
import sys
import subprocess
import os
import pandas as pd

from mapper import mapper
import rescore

# TODO config file for msgf+ -> config file for ms2pip
# Path to MSGFPlus & ms2pip
MSGF_DIR = "/home/compomics/software/MSGFPlus"
MS2PIP_DIR = "ms2pip_c/"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run MSGF+ to get PSMs, MS2PIP to extract spectral features, and Percolator to rescore PSMs')
    parser.add_argument('spec_file', metavar='spectrum-file',
                        help='file containing MS2 spectra (MGF)')
    parser.add_argument('fasta_file', metavar='FASTA-file',
                        help='file containing protein sequences')
    parser.add_argument('-m', '--mods', metavar='FILE', action="store", default='',
                        dest='modsfile', help='Mods.txt file for MSGF+')
    parser.add_argument('-f', '--frag', metavar='frag_method', action="store", default='HCD',
                        dest='frag', help='fragmentation method (CID or HCD), default HCD')

    args = parser.parse_args()
    """
    # Run MSGF+
    rescore.run_msgfplus(MSGF_DIR, args.spec_file, args.spec_file,
                         args.fasta_file, args.modsfile, args.frag)

    # Convert .mzid to pin, for percolator. XXX is the decoy pattern from MSGF+
    convert_command = "msgf2pin -P XXX {}.mzid > {}.pin".format(args.spec_file.rstrip(".mgf"), args.spec_file.rstrip(".mgf"))
    sys.stdout.write("Converting .mzid file to pin file: {} \n".format(convert_command))
    sys.stdout.flush()
    subprocess.run(convert_command, shell=True)

    # "Proteins" column has tab-separated values which makes the file cumbersome
    #  to read. mapper.fix_pin_tabs replaces those tabs with ";"
    sys.stdout.write('Fixing tabs on pin file... ')
    sys.stdout.flush()
    mapper.fix_pin_tabs(args.spec_file.rstrip(".mgf") + ".pin")
    os.remove(args.spec_file.rstrip(".mgf") + ".pin")
    os.rename(args.spec_file.rstrip(".mgf") + "_fixed.pin", args.spec_file.rstrip(".mgf") + ".pin")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    # Percolator generates its own spectrum ID, but we want it to match the mgf
    # file's TITLE.
    sys.stdout.write("Adding TITLE to pin file... ")
    sys.stdout.flush()
    mapper.map_mgf_title(args.spec_file.rstrip(".mgf") + ".pin", args.spec_file.rstrip(".mgf") + ".mzid")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()
    """
    # Create & write PEPREC file from the pin file
    sys.stdout.write("Generating PEPREC files... ")
    sys.stdout.flush()
    rescore.make_pepfile(args.spec_file.rstrip(".mgf") + ".pin")
    os.rename(args.spec_file.rstrip(".mgf") + ".pin.PEPREC", args.spec_file.rstrip(".mgf") + ".PEPREC")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    # Run ms2pip
    # TODO change config file to use CID or HCD models based on args.frag
    ms2pip_command = "python {}ms2pipC.py {} -c {} -s {}".format(MS2PIP_DIR, args.spec_file.rstrip(".mgf") + ".PEPREC", MS2PIP_DIR + 'config.file', args.spec_file)
    sys.stdout.write("Running ms2pip: {} \n".format(ms2pip_command))
    sys.stdout.flush()
    subprocess.run(ms2pip_command, shell=True)

    sys.stdout.write("Calculating features from predicted spectra... ")
    sys.stdout.flush()
    rescore.calculate_features(args.spec_file.rstrip(".mgf") + "_pred_and_emp.csv", args.spec_file.rstrip(".mgf") + "_all_features.csv")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    sys.stdout.write("Generating pin files with different features... ")
    sys.stdout.flush()
    rescore.join_features(args.spec_file.rstrip(".mgf") + "_all_features.csv", args.spec_file.rstrip(".mgf") + ".pin")
    rescore.write_pin_files(args.spec_file.rstrip(".mgf") + "_all_features.csv", args.spec_file.rstrip(".mgf"))
    # os.remove(args.spec_file + ".pin")
    # os.remove(args.spec_file + "_all_features.csv")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    # Run Percolator with different feature subsets
    for subset in ['_rescore', '_percolator', '_all_features']:
        fname = args.spec_file.rstrip(".mgf") + subset
        percolator_cmd = "percolator {} -r {} -M {} -v 0 -U\n".format(fname + ".pin", fname + ".pout", fname + ".pout_dec")
        sys.stdout.write("Running Percolator: {} \n".format(percolator_cmd))
        subprocess.run(percolator_cmd, shell=True)
        rescore.format_output(fname+".pout", fname+".pout_dec", fname+"_output_plots.png", fig=False)
        sys.stdout.flush()

    sys.stdout.write('All done!\n')
