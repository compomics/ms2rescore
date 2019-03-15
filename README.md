# ReScore

Use features calculated from comparing experimental spectra with computationally generated spectra (see [MS2PIP](https://github.com/sdgroeve/ms2pip_c)) to re-score peptide identifications using [Percolator](https://github.com/percolator/percolator/).

This branch starts with a FASTA and an MGF file, runs a search using MS-GF+, generates the MS2PIP input file and runs MS2PIP. From the MS2PIP output the spectral comparison features are generated. Three pin files are written: one with the default features, one with the MS2PIP features and a third with all of them in a single file. Percolator is called on each of these files, creating two pout files (one for targets and one for decoys) for each pin file.

## Prerequisites

- [MS-GF+](https://omics.pnl.gov/software/ms-gf)
- [Percolator](https://github.com/percolator/percolator/)
- Python 3

  - numpy
  - pandas
  - scikit-learn
  - scipy
  - tqdm

### Installing

Clone this repository. This includes the submodules [mapper](https://github.com/anasilviacs/mapper/tree/0ee46adcbb20a118a8274908255cc8b3f95a51db) and [MS2PIP](https://github.com/anasilviacs/ms2pip_c/tree/6f037dc2d0797cd25061aaed8091d625123971e1). As such, the repository should be called as:

```
git clone --recurse-submodules git://github.com/anasilviacs/rescore.git
```

Go into the folder `ms2pip_c` and run the bash script `compile.sh`. For more details, please refer to the `ms2pip_c` repository.

## Usage

To run the entire pipeline from the beginning, you need a spectrum file in the [MGF format](http://www.matrixscience.com/help/data_file_help.html) and a protein database in the [FASTA format](https://zhanglab.ccmb.med.umich.edu/FASTA/).

Run the pipeline as follows:

```
python driver.py <mgf file> <fasta file> <config file>
```

- `<mgf file>` is the spectrum file
- `<fasta file>` is the database file
- `<config file>` is json with configurations for MS-GF+, MS2PIP and Percolator. A reference file is included in this repository, that can be adapted.

Only MS-GF+ is directly supported, but if you want to use a different search engine and re-score the obtained identifications, you can use the helper functions in `rescore.py` and `mapper.py` to get the necessary files to run MS2PIP, which are a `PEPREC` file and the spectrum file, and to generate the pin files with the different subset of features. Afterwards you can run Percolator as normal.

### Output

Several files are created when the entire pipeline is ran. Their names are all built from the `<mgf file>` name and are stored in that file's folder.

- `<mgf file>.mzid` is the identification file from the MS-GF+ search
- `<mgf file>.PEPREC` is the list of peptides that MS2PIP will predict spectra and calculate features from
- `<mgf file>.PEPREC.rescore_features.csv` are the features calculated by MS2PIP
- three `.pin` files, which are the standard Percolator input format:

  - `<mgf file>_rescore.pin` with the MS2PIP features
  - `<mgf file>_percolator.pin` with the default Percolator features
  - `<mgf file>_all_features.pin` with all the features

- three corresponding `.pout` files, with Percolator's output for target PSMs
- three corresponding `.pout_dec` files, with Percolator's output for decoy PSMs
- three corresponding `.png` files, which show the final distribution of scores, q-values and PEPs.
