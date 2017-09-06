"""
Module that calls the necessary functions and run the re-scoring algorithm.
"""

import argparse
import sys
import subprocess
import pandas as pd

from mapper import mapper
import rescore

# TODO better file names
# TODO ms2pip could be updated to python3
# Path to MSGFPlus & ms2pip
MSGF_DIR = "/home/compomics/software/MSGFPlus"
# MS2PIP_DIR = "/home/compomics/local/rescore-ms2pip/ms2pip_c"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MSGF+ and Percolator')
    parser.add_argument('spec_file', metavar='spectrum-file',
                        help='file containing MS2 spectra (MGF,PKL,DTA,mzXML,mzDATA or mzML)')
    parser.add_argument('fasta_file', metavar='FASTA-file',
                        help='file containing protein sequences')
    parser.add_argument('-m', '--mods', metavar='FILE', action="store", default='',
                        dest='modsfile', help='Mods.txt file for MSGF+')
    parser.add_argument('-f', '--frag', metavar='frag_method', action="store", default='HCD',
                        dest='frag', help='fragmentation method (CID or HCD), default HCD')

    args = parser.parse_args()

    # Run MSGF+
    rescore.run_msgfplus(MSGF_DIR, args.spec_file + ".msgf_out", args.spec_file,
                 args.fasta_file, args.modsfile, args.frag)

    # Convert .mzid to pin, for percolator. XXX is the decoy pattern from MSGF+
    convert_command = "msgf2pin -P XXX %s.mzid > %s.pin" % (
        args.spec_file + ".msgf_out", args.spec_file + ".msgf_out")
    sys.stdout.write("Converting .mzid file to pin file: {} \n".format(convert_command))
    sys.stdout.flush()
    subprocess.run(convert_command, shell=True)

    # Add mgf TITLE column to pin file
    sys.stdout.write('Fixing tabs on pin file... ')
    sys.stdout.flush()
    mapper.fix_pin_tabs(args.spec_file + ".msgf_out.pin")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    sys.stdout.write('Parsing pin file... ')
    sys.stdout.flush()
    pin = pd.read_csv(args.spec_file + ".msgf_out_fixed.pin", header=0, skiprows=[1], sep='\t')
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    sys.stdout.write("Adding TITLE to pin file... ")
    sys.stdout.flush()
    # now YOU are the bottleneck!
    pin = mapper.map_mgf_title(pin, args.spec_file + ".msgf_out.mzid")
    pin.to_csv(args.spec_file + ".titles.pin", sep='\t', index=False)
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    # Create & write PEPREC file from the pin file
    sys.stdout.write("Generating PEPREC files... ")
    sys.stdout.flush()
    peprec = rescore.make_pepfile(args.spec_file + ".titles.pin")
    rescore.write_PEPREC(peprec, args.spec_file + ".titles.pin")
    sys.stdout.write('Done! \n')
    sys.stdout.flush()

    # Run ms2pip_rescore: ms2pip is written in python 2!!!
    ms2pip_command = "python {}/ms2pipC.py {} -c {} -s {} -R".format(MS2PIP_DIR, args.spec_file + ".titles.pin.PEPREC", MS2PIP_DIR + '/config.file', args.spec_file)
    # sys.stdout.write("Running ms2pip with the rescore option: {} \n".format(ms2pip_command))
    sys.stdout.write('Please run ms2pip with the following command: {}'.format(ms2pip_command))
    sys.stdout.flush()
    # subprocess.run(ms2pip_command, shell=True)

    # features = rescore.join_features(args.spec_file + '.titles.pin.PEPREC_rescore_features.csv', args.spec_file + ".titles.pin")
    # rescore.write_pin_files(features, args.spec_file)

    # Run Percolator
