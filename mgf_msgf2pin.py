"""
Create a pin (Percolator input) file
First MSGF+ is used to search <spec_file> against <fasta_file>. For now this
will be a concatenated search, i.e. MSGF+ generates a decoy database and searches
 the spectra against both at the same time (one PSM per spectrum).
Second, the resulting mzid fileis converted into a pin file. Then a column with
the correspondence to mgf TITLE is added to the pin file.
"""

import subprocess
import sys
import argparse
from mapper import mapper # shouldn't have to do this, check __init__.py

def run_msgfplus(outfile, mgffile, fastafile, modsfile, frag='HCD'):
    """
    Runs MSGFPlus with some fixed settings. Requires path to output file, mgf,
    fasta and modifications file.
    """
    if frag == 'HCD':
        m = 3
        inst = 1
    elif frag == 'CID':
        m = 1
        inst = 0
    if modsfile != '':
        mods = '-mod {} '.format(modsfile)
    else:
        mods = ''
    msgf_command = "java -Xmx8000M -jar {}/MSGFPlus.jar {}-s {} -d {} -o {} -t \
        10ppm -tda 1 -m {} -inst {} -minLength 8 -minCharge 2 -maxCharge 4 -n \
        1 -addFeatures 1 -protocol 0 -thread 23".format(MSGF_DIR, mods, mgffile,
        fastafile, outfile + ".mzid", m, inst)
    sys.stdout.write("Running search with MSGF+: {}".format(msgf_command))
    sys.stdout.flush
    subprocess.run(msgf_command, shell=True)

if __name__ == "__main__":
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

    # Path to MSGFPlus - this should come from a config file
    MSGF_DIR = "/home/compomics/software/MSGFPlus"

    # Run MSGF+
    run_msgfplus(args.spec_file + ".target", args.spec_file,
                 args.fasta_file, args.modsfile, args.frag)

    # Convert .mzid to pin, for percolator. XXX is the decoy pattern from MSGF+
    convert_command = "msgf2pin -P XXX %s.mzid > %s.pin" % (
        args.spec_file + ".target", args.spec_file + ".target")
    sys.stdout.write("Converting .mzid file to pin file: {}".format(convert_command))
    sys.stdout.flush()
    subprocess.run(convert_command, shell=True)

    # Add mgf TITLE column to pin file
    pin = mapper.lazy_pin_parser(args.spec_file + ".target.pin")
    pin = mapper.map_mgf_title(pin, args.spec_file + ".target")

    pin.to_csv(args.spec_file + ".target.titles.pin", sep='\t', index=False)
