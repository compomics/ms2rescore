"""
Create pin files that include the ms2pip features
"""

import argparse
import pandas as pd

def join_features(path_to_target_features, path_to_pin, path_to_decoy_features=None):
    """
    Create one big DataFrame with all the Percolator and rescore features
    """
    # read pin file - should not need the lazy pin parser as this pin already has the TITLE which means it was processed by mapper
    pin = pd.read_csv(path_to_pin, sep='\t')

    # Read rescore_features.csv file and fillna
    rescore_targets = pd.read_csv(path_to_target_features)
    rescore_targets = rescore_targets.fillna(0)

    # If not concat searches, do that for target and decoy files
    if path_to_decoy_features != None:
        rescore_decoys = pd.read_csv(args.decoys)
        rescore_decoys = rescore_decoys.fillna(0)

        # join target and decoy tables
        all_features = pd.concat([rescore_decoys.merge(pin[pin.Label == -1], left_on='spec_id', right_on='TITLE'), rescore_targets.merge(pin[pin.Label == 1], left_on='spec_id', right_on='TITLE')])
    else:
        all_features = rescore_targets.merge(pin, left_on='spec_id', right_on='TITLE')

    return all_features

def write_pin_files(all_features):
    """
    write five copies of the pin files:
    _only_rescore.pin with only the rescore features
    _all_percolator.pin with all the percolator features
    _percolator_default.pin with only the default percolator features
    _all_features.pin with all the rescore and percolator features
    _default_and_rescore.pin with the default percolator features and all of the rescore
    """
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

    return None

if __name__ == '__main__':
    # Arguments: percolator features, .mzid
    parser = argparse.ArgumentParser(description="Get a pin file with Percolator and ms2pip features")
    parser.add_argument('-t', dest='targets', help='Path to target search ms2pip features file')
    parser.add_argument('-d', dest='decoys', help='Path to decoy search ms2pip features file')
    parser.add_argument('-p', dest='pin', help='Path to pin file')

    args = parser.parse_args()

    features = join_features(args.targets, args.pin)
    write_pin_files(features)
