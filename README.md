# ReScore

Use features calculated from comparing experimental spectra with computationally generated spectra (see [MS2PIP](https://github.com/sdgroeve/ms2pip_c)) to re-score peptide identifications using [Percolator](https://github.com/percolator/percolator/).

## Prerequisites

- [MS-GF+](https://omics.pnl.gov/software/ms-gf): currently, MS-GF+ is the only supported search engine
- [Percolator](https://github.com/percolator/percolator/)
- Python 3

  - numpy==1.13.3
  - pandas==0.20.3
  - scikit-learn==0.19.0
  - scipy==0.19.1
  - xmltodict==0.11.0

### Installing

Clone this repository. This includes the submodules [mapper](https://github.com/anasilviacs/mapper/tree/0ee46adcbb20a118a8274908255cc8b3f95a51db) and [MS2PIP](https://github.com/anasilviacs/ms2pip_c/tree/6f037dc2d0797cd25061aaed8091d625123971e1).

Go into the folder `ms2pip_c` and run the bash script `compile.sh` to compile the prediction models. To make sure this step went smoothly, you can install the python package `pytest`, go into the folder `tests` and type `pytest`. Successfully passing the tests means that everything went well.

## Usage

To run the entire pipeline from the beginning, you need a spectrum file in the [MGF format](http://www.matrixscience.com/help/data_file_help.html) and a protein database in the [FASTA format](https://zhanglab.ccmb.med.umich.edu/FASTA/). Add the path to the folder where `MSGFPlus.jar` is to the `driver.py` file:

```
MSGF_DIR = "/path/to/MSGFPlus"
```

Run the pipeline as follows:

```
python driver.py <mgf file> <fasta file> -m <mods file> -f <frag method>
```

- `<mgf file>` is the spectrum file
- `<fasta file>` is the database file
- `<mods file>` is a file with modifications as per the MSGFPlus specifications (optional)
- `<frag method>` is the fragmentation method (optional, HCD by default)

### Output

Several files are created when the entire pipeline is ran. Their names are all built from the `<mgf file>` name.

- `<mgf file>.mzid` is the identification file from the MS-GF+ search
- `<mgf file>.PEPREC` is the list of peptides that MS2PIP will predict spectra and calculate features from
- `<mgf file>.PEPREC.rescore_features.csv` are the features calculated by MS2PIP
- several `.pin` files, which are the standard Percolator input format:

  - `<mgf file>_all_percolator.pin` with all the Percolator features
  - `<mgf file>_only_rescore.pin` with the MS2PIP features
  - `<mgf file>_percolator_default.pin` with the default Percolator features
  - `<mgf file>_default_and_rescore.pin` with the default Percolator features and the MS2PIP features
  - `<mgf file>_all_features.pin` with all the features

- several `.pout` files, with Percolator's output for target PSMs
- several `.pout_dec` files, with Percolator's output for decoy PSMs
