# ReScore

Use features calculated from comparing experimental spectra with computationally generated spectra (see [MS2PIP](https://github.com/compomics/ms2pip_c)) to re-score peptide identifications using [Percolator](https://github.com/percolator/percolator/).

On this branch, scripts that can be used to generate the necessary files are also provided. Currently, we provide scripts for:
- [MS-GF+](https://omics.pnl.gov/software/ms-gf) (Start with `.mgf` and `.fasta` file, or from `.mzid` identifications file)
- [MaxQuant](https://www.maxquant.org/) (Start from `msms.txt` identification file and `.mgf` files. Be sure to export without FDR filtering!)


## Prerequisites
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
git clone --recurse-submodules git://github.com/compomics/rescore.git
```

Make sure you are in the correct branch :)

If you missed the required usage of the `--recurse-submodules` flag, you can go into the submodules folders' and run the following command:

```
git submodule update --init
```

Go into the folder `ms2pip_c` and run the bash script `compile.sh`. For more details, please refer to the `ms2pip_c` repository.

## Usage

Refer to the [wiki pages](https://github.com/compomics/rescore/wiki) for instructions on how to use our helper scripts to generate the necessary input files for the ReScore pipeline.

To run the entire pipeline from the beginning, you need a spectrum file in the [MGF format](http://www.matrixscience.com/help/data_file_help.html) and a peptide list in the [PEPREC format](https://github.com/anasilviacs/ms2pip_c/tree/7b89618f236c84aed3c171132f690556c757b6b5). The file should have an additional column, `Label`, where target PSMs have a value of `1` and decoys of `-1`, and `Protein`, with the protein identifiers the PSM is associated with (this is optional for the execution of the pipeline, but likely important for posterior analysis).

Run the pipeline as follows:

```
python driver.py <mgf file> <peprec file> <config file>
```

- `<mgf file>` is the spectrum file
- `<peprec file>` is the peptide list file
- `<config file>` is json with configurations for MS2PIP. A reference file that can be adapted is included in this repository.

An example configuration file is included (`config.json`). The modifications in the `<PEPREC file>` should match the modification names in `config.json`. This file already includes several modifications, but more can be added by following the same structure.


### Output

Several intermediate files are created when the entire pipeline is ran. Their names are all built based on the `<PEPREC file>` name and are stored in that file's folder. The most relevant files are the Percolator INput files, `.pin`. With this pipeline, three files are written: one with only MS2PIP-based features, one with search engine-based features, and a third where both sets of features are combined.

`Percolator` is ran on each `.pin`, using the settings which have been specified in the configurations file. For details on these settings and on how to run `Percolator` please refer to its [wiki pages](https://github.com/percolator/percolator/wiki).

From each `Percolator` execution, three files are generated:

- `<file>.pout` with the output regarding target PSMs
- `<file>.pout_dec` with the output regarding decoy PSMs
- `<file>.weights` where the internal feature weights used by `Percolator`'s scoring function are stored.
