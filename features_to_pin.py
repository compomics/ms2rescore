# Read rescore_features.csv file and fillna
rescore_targets = pd.read_csv('../pyrococcus/pepfile_targets_rescore_features.csv')
rescore_targets = rescore_targets.fillna(0)

# If not concat searches, do that for target and decoy files
rescore_decoys = pd.read_csv('../pyrococcus/pepfile_decoys_rescore_features.csv')
rescore_decoys = rescore_decoys.fillna(0)

# join target and decoy tables
all_features = pd.concat([rescore_decoys.merge(pin[pin.Label == -1], left_on='spec_id', right_on='TITLE'),
                          rescore_targets.merge(pin[pin.Label == 1], left_on='spec_id', right_on='TITLE')])

# read pin file
pin = pd.read_csv('../pyrococcus/msgf_features.pin', sep='\t')

# columns to save
rescore_features = list(all_features.columns[4:42])
percolator_features = list(all_features.columns[45:-5])
percolator_default = percolator_features[:27]

# Writing files with appropriate columns
all_features[['SpecId', 'Label', 'ScanNr'] + rescore_features + ['Peptide', 'Proteins']].to_csv('../pyrococcus/only_rescore.pin', sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_features + ['Peptide', 'Proteins']].to_csv('../pyrococcus/all_percolator.pin', sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_default + ['Peptide', 'Proteins']].to_csv('../pyrococcus/percolator_default.pin', sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_features + rescore_features + ['Peptide', 'Proteins']].to_csv('../pyrococcus/all_features.pin', sep='\t', index=False)
all_features[['SpecId', 'Label', 'ScanNr'] + percolator_default + rescore_features + ['Peptide', 'Proteins']].to_csv('../pyrococcus/default_plus_rescore.pin', sep='\t', index=False)
