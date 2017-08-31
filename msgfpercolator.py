import sys
import os
import argparse

"""
Run MSGFPlus and percolator
"""

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

def run_MSGFplus(outfile, mgffile, fastafile, modsfile, frag='HCD'):
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
    msgf_command = "java -Xmx8000M -jar {}/MSGFPlus.jar {}-s {} -d {} -o {} -t 10ppm -tda 1 -m {} -inst {} -minLength 8 -minCharge 2 -maxCharge 4 -n 1 -addFeatures 1 -protocol 0 -thread 23".format(
        msgfdir, mods, mgffile, fastafile, outfile + ".mzid", m, inst)
    print(msgf_command)
    os.system(msgf_command)

# Run MSGF+
run_MSGFplus(args.spec_file + ".target", args.spec_file, args.fasta_file, args.modsfile, args.frag)

# Convert .mzid to pin, for percolator
convert_command = "msgf2pin -P XXX %s.mzid > %s.pin" % (args.spec_file + ".target", args.spec_file + ".target")
os.system(convert_command)

# Maps mgf spectrum to pin
# This can perhaps be replaced with a call to mapper
id_map = {}
tid = ""
with open(args.spec_file + ".target.mzid") as f:
    for row in f:
        if "<SpectrumIdentificationItem" in row:
            l = row.rstrip().split('id=')
            tid = l[1][1:-2]
        if 'name="spectrum title"' in row:
            l = row.rstrip().split('value=')[1].split(' ')[0]
            #id_map[tid] = l[1:-1]
            id_map[tid] = l[1:]

pin_map = {}
header = ""
with open("%s" % (args.spec_file + ".target.pin")) as f:
    row = f.readline()
    header = row.split('\t')[0:-1]
    for row in f:
        l = row.rstrip().split('\t')
        pin_map[l[0]] = l[0:len(header)]

# Run Percolator
command = "percolator -U %s > %s.out" % (args.spec_file + ".target.pin", args.spec_file + ".target.pin")
os.system(command)

fout2 = open(args.spec_file + ".msgfout", "w")
with open("%s.out" % (args.spec_file + ".target.pin")) as f:
    row = f.readline()
    fout2.write('\t'.join(header) + '\t' + row)
    for row in f:
        l = row.rstrip().split('\t')
        l0 = l[0]
        tmp = '_'.join(l[0].split('_')[-6:-3])
        if tmp in id_map:
            l[0] = id_map[tmp]
            prot = l[5]
            if len(l) > 6:
                for i in range(6, len(l)):
                    prot += "|" + l[i]
                l[5] = prot
                fout2.write('\t'.join(pin_map[l0]) + '\t' + '\t'.join(l[:6]) + '\n')
