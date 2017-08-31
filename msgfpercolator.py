import sys
import os

import argparse
import tempfile

import globals

msgfdir = "/home/compomics/local/MSGFPlus"


def run_MSGFplusHCD(outfile, mgffile, fastafile, modsfile):


    msgf_command = "java -Xmx8000M -jar %s/MSGFPlus.jar -mod %s -s %s -d %s \
            -o %s -t 10ppm -tda 1 -m 3 -inst 1 -minLength 8 \
            -minCharge 2 -maxCharge 4 -n 1 -addFeatures 1 -protocol 0 -thread 23" \
            % (msgfdir, modsfile, mgffile, fastafile, outfile + ".mzid")

    os.system(msgf_command)


def run_MSGFplusCID(outfile, mgffile, fastafile, modsfile):


    msgf_command = "java -Xmx8000M -jar %s/MSGFPlus.jar -mod %s -s %s -d %s \
            -o %s -t 10ppm -tda 1 -m 1 -inst 0 -minLength 8 \
            -minCharge 2 -maxCharge 4 -n 1 -addFeatures 1 -protocol 0 -thread 23" \
            % (msgfdir, modsfile, mgffile, fastafile, outfile + ".mzid")
    os.system(msgf_command)


parser = argparse.ArgumentParser(description='MSGF+')
parser.add_argument('spec_file', metavar='spectrum-file',
                    help='file containing MS2 spectra (MGF,PKL,DTA,mzXML,mzDATA or mzML)')
parser.add_argument('fasta_file', metavar='FASTA-file',
                    help='file containing protein sequences')
parser.add_argument('-m', '--mods', metavar='FILE', action="store",
                    dest='modsfile',
                    help='Mods.txt file for MSGF+')

args = parser.parse_args()

run_MSGFplusHCD(args.spec_file + ".target", args.spec_file,
                args.fasta_file, args.modsfile)
convert_command = "msgf2pin -P XXX %s.mzid > %s.pin" %
(args.spec_file + ".target", args.spec_file + ".target")
os.system(convert_command)

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
