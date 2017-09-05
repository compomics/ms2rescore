"""
Module that calls the necessary functions and run the re-scoring algorithm.
"""

import argparse
import sys
import subprocess

from mapper import mapper
import rescore

# TODO better file names
# TODO call Percolator from this script
# TODO ms2pip config file - set a default?
# Path to MSGFPlus & ms2pip - these should come from a config file
MSGF_DIR = "/home/compomics/software/MSGFPlus"
MS2PIP_DIR = "ms2pip_c" # redundant bc when you clone this, you clone ms2pip_c

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
    sys.stdout.write("Lazily parsing pin file: {} \n".format(args.spec_file + ".msgf_out.pin"))
    sys.stdout.flush()
    pin = mapper.lazy_pin_parser(args.spec_file + ".msgf_out.pin")

    sys.stdout.write("Adding TITLE to pin file \n")
    sys.stdout.flush()
    pin = mapper.map_mgf_title(pin, args.spec_file + ".msgf_out.pin")

    # Write pin file
    sys.stdout.write("Writing pin file \n")
    sys.stdout.flush()
    pin.to_csv(args.spec_file + ".titles.pin", sep='\t', index=False)

    # Create & write PEPREC file from the pin file
    sys.stdout.write("Generating PEPREC files \n")
    sys.stdout.flush()
    peprec = rescore.make_pepfile(args.spec_file + ".titles.pin")
    rescore.write_PEPREC(peprec, args.spec_file + ".titles.pin")

    # Run ms2pip_rescore
    ms2pip_command = "python {}/ms2pipC.py {} -c {someconfigfile} -s {}".format(MS2PIP_DIR, args.spec_file + ".titles.pin.PEPREC", config, args.spec_file, -R)
    sys.stdout.write("Running ms2pip with the rescore option: {}".format(ms2pip_command))
    sys.stdout.flush()
    subprocess.run(ms2pip_command, shell=True)

    features = rescore.join_features(args.pep_file + '.titles.pin.PEPREC_rescore_features.csv', args.pin)
    rescore.write_pin_files(features)

    # Run Percolator
