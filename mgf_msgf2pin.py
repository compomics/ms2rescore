"""
Create a pin (Percolator input) file
First MSGF+ is used to search <spec_file> against <fasta_file>. For now this
will be a concatenated search, i.e. MSGF+ generates a decoy database and searches
 the spectra against both at the same time (one PSM per spectrum).
Second, the resulting mzid fileis converted into a pin file. Then a column with
the correspondence to mgf TITLE is added to the pin file.
"""
#TODO this should run msgf+ (concat or not concat) and msgf2pin and mapper.
import os
import argparse

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

# Path to MSGFPlus
msgfdir = "/home/compomics/software/MSGFPlus"


def run_msgfplus(outfile, mgffile, fastafile, modsfile, frag='HCD'):
    """
    Runs MSGFPlus
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
        1 -addFeatures 1 -protocol 0 -thread 23".format(msgfdir, mods, mgffile,
        fastafile, outfile + ".mzid", m, inst)
    print(msgf_command)
    os.system(msgf_command)


# Run MSGF+
run_msgfplus(args.spec_file + ".target", args.spec_file,
             args.fasta_file, args.modsfile, args.frag)

# Convert .mzid to pin, for percolator. XXX is the decoy pattern from MSGF+
convert_command = "msgf2pin -P XXX %s.mzid > %s.pin" % (
    args.spec_file + ".target", args.spec_file + ".target")
print(convert_command)
os.system(convert_command)

# Add mgf TITLE column to pin file
command = "python mapper -m {} -p {}".format(
    args.spec_file + '.target.mzid', args.spec_file + '.target.pin')
print(command)
os.system(command)
