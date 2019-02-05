import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="With a (fixed) percolator pin file and an ms2pip pin file, make new pin files")
    parser.add_argument("perc_pin")
    parser.add_argument("ms2pip_pin")

    args = parser.parse_args()

    sel_feats = ['ionb_mse', 'ionb_abs_diff_Q1_norm', 'min_abs_diff',
     'ionb_max_abs_diff_norm', 'iony_pearson_norm', 'spec_pearson',
     'spec_pearson_norm', 'max_abs_diff', 'abs_diff_Q1', 'ionb_min_abs_diff',
     'ionb_abs_diff_Q3_norm', 'abs_diff_Q3', 'ionb_pearson_norm',
     'iony_abs_diff_Q2', 'ionb_mse_norm', 'ionb_abs_diff_Q1',
     'iony_min_abs_diff_norm', 'ionb_abs_diff_Q3', 'iony_pearson', 'cos_ionb_norm',
     'spec_mse', 'max_abs_diff_norm', 'iony_abs_diff_Q1', 'ionb_min_abs_diff_norm',
     'ionb_abs_diff_Q2_norm', 'ionb_max_abs_diff', 'spec_mse_norm', 'abs_diff_Q2_norm',
     'cos_iony_norm', 'ionb_abs_diff_Q2', 'dotprod_ionb', 'min_abs_diff_norm',
     'abs_diff_Q2', 'dotprod_norm', 'cos_iony', 'min_abs_diff_iontype',
     'ionb_mean_abs_diff', 'mean_abs_diff', 'max_abs_diff_iontype', 'abs_diff_Q1_norm',
     'dotprod', 'iony_abs_diff_Q3', 'cos_norm', 'iony_min_abs_diff']


    perc_pin = pd.read_csv(args.perc_pin, sep='\t')
    perc_pin['SpecId'] = perc_pin.TITLE
    perc_pin.pop('TITLE')
    perc_pin.set_index('SpecId', inplace=True)

    ms2pip_pin = pd.read_csv(args.ms2pip_pin, sep='\t')
    ms2pip_pin.set_index('SpecId', inplace=True)

    perc_feats = list(perc_pin.columns[2:-2])

    pin_both = pd.concat([perc_pin, ms2pip_pin], axis=1, join='inner')
    pin_both = pin_both.loc[:,~pin_both.columns.duplicated()]

    if pin_both.shape[0] == 0:
        print(pin_both.head())
        exit(1)

    fname = args.perc_pin.rstrip('.pin')

    pin_both.loc[:, ['Label', 'ScanNr'] + perc_feats + sel_feats + ['Peptide', 'Proteins']].to_csv(fname + '_both.pin', sep='\t', index=True)
    pin_both.loc[:, ['Label', 'ScanNr'] + perc_feats + ['Peptide', 'Proteins']].to_csv(fname + '_perc.pin', sep='\t', index=True)
    pin_both.loc[:, ['Label', 'ScanNr'] + sel_feats + ['Peptide', 'Proteins']].to_csv(fname + '_resc.pin', sep='\t', index=True)
