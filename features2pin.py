"""
prepare pin files with the ms2pip and percolator features
"""
import pandas as pd
import argparse

# Arguments: percolator features, .mzid
parser = argparse.ArgumentParser(description="Get a pin file with Percolator and ms2pip features")
parser.add_argument('-t', dest='targets', help='Path to target search ms2pip features file')
parser.add_argument('-d', dest='decoys', help='Path to decoy search ms2pip features file')
parser.add_argument('-p', dest='pin', help='Path to pin file')

args = parser.parse_args()

# Read rescore_features.csv file and fillna
rescore_targets = pd.read_csv(args.targets)
rescore_targets = rescore_targets.fillna(0)

# If not concat searches, do that for target and decoy files
rescore_decoys = pd.read_csv(args.decoys)
rescore_decoys = rescore_decoys.fillna(0)

# join target and decoy tables
all_features = pd.concat([rescore_decoys.merge(pin[pin.Label == -1], left_on='spec_id', right_on='TITLE'),
                          rescore_targets.merge(pin[pin.Label == 1], left_on='spec_id', right_on='TITLE')])

# read pin file
pin = pd.read_csv(args.pin, sep='\t')

# columns to save
rescore_features = list(all_features.columns[4:42])
percolator_features = list(all_features.columns[45:-5])
percolator_default = percolator_features[:27]

# Writing files with appropriate columns
all_features[['SpecId', 'Label', 'ScanNr'] + rescore_features + ['Peptide', 'Proteins']].to_csv('{}_only_rescore.pin'.format(args.pin.split('.')[:-1]), sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_features + ['Peptide', 'Proteins']].to_csv('{}_all_percolator.pin'.format(args.pin.split('.')[:-1]), sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_default + ['Peptide', 'Proteins']].to_csv('{}_percolator_default.pin'.format(args.pin.split('.')[:-1]), sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_features + rescore_features + ['Peptide', 'Proteins']].to_csv('{}_all_features.pin'.format(args.pin.split('.')[:-1]), sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_default + rescore_features + ['Peptide', 'Proteins']].to_csv('{}_default_and_rescore.pin'.format(args.pin.split('.')[:-1]), sep='\t', index=False)
