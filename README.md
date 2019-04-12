# MS²ReScore

Use features calculated from comparing experimental spectra with computationally
generated spectra (see [MS²PIP](https://github.com/compomics/ms2pip_c)) to
rescore peptide identifications using
[Percolator](https://github.com/percolator/percolator/).

On this branch, multiple pipelines can be run, depending on your input format:
- [MaxQuant](https://www.maxquant.org/): Start from `msms.txt` identification
file and directory with `.mgf` files. Be sure to export without FDR filtering!

Work in progress:
- [MS-GF+](https://omics.pnl.gov/software/ms-gf): Start with `.mgf` and `.fasta`
file, or from `.mzid` identifications file)
- PEPREC: Start with [MS²PIP](https://github.com/compomics/ms2pip_c) `.peprec`
and `.mgf` file.

## Prerequisites
- [Percolator](https://github.com/percolator/percolator/) needs to be callable
with the `percolator` command (tested with version 3.02.1)
- Python 3 on Linux
  - numpy
  - pandas
  - scikit-learn
  - scipy
  - tqdm

## Installation
Clone this repository. This includes the submodules
[mapper](https://github.com/anasilviacs/mapper/tree/0ee46adcbb20a118a8274908255cc8b3f95a51db)
and [MS2PIP](https://github.com/anasilviacs/ms2pip_c/tree/6f037dc2d0797cd25061aaed8091d625123971e1).
As such, the repository should be called as:

```
git clone --recurse-submodules git://github.com/compomics/rescore.git
```

Go to the newly made directory and switch to the correct branch:
```
cd rescore
git checkout dev
```

If you missed the required usage of the `--recurse-submodules` flag, you can go
into the submodules folders' and run the following command:
```
git submodule update --init
```

Go into the folder `ms2pip_c` and run the bash script `compile.sh` or use the 
precompiled files. For more details, please refer to the `ms2pip_c` repository.

Install MS²ReScore with
```
pip install .
```

## Usage
### Command line interface
Run MS²ReScore as follows:
```
usage: ms2rescore [-h] [-o FILE] [-l LEVEL] config-file

MS²ReScore: Rescoring of PSMs with predicted MS² peak intensities.

positional arguments:
  config-file  json MS2ReScore configuration file. See README.md

optional arguments:
  -h, --help   show this help message and exit
  -o FILE      Name for output files (default: `ms2rescore_out`)
  -l LEVEL     Logging level (default: `info`)
  ```

### Config file
The main argument is a link to the config file. This JSON file should contain
all required information for MS²ReScore to be run properly. An example file for
each pipeline is provided in this repository.

The config file contains three main top level keys (`general`, `ms2pip` and
`percolator`) and a key for each pipeline (e.g. `maxquant`). 

#### General
- `pipeline`: pipeline to use (eg `"MaxQuant"`)
- `feature_sets`: list with feature sets to use for rescoring. Options are:
    - `all` = both search engine features and MS²PIP features
    - `ms2pip` = only MS²PIP features
    - `searchengine` = only search engine features (classic Percolator)
- `run_percolator`: bool, wether or not to call Percolator from the MS²ReScore
pipeline
- `keep_tmp_files`: Keep temporary files or not (e.g. MS²PIP output). These
files can be used for a more in-depth data-analysis. If set to true, only the
PIN files and the Percolator output are kept.
- `num_cpu`: number of CPU cores to use when using parallel processing
For example:
```json
"general":{
    "pipeline":"MaxQuant",
    "feature_sets":["all", "ms2pip", "searchengine"],
    "run_percolator":true,
    "keep_tmp_files":false,
    "num_cpu":"24"
}
```

#### MS2PIP
```json
"ms2pip":{
    "dir":"ms2pip_c",
    "frag":"HCD",
    "frag_error":0.02,
    "modifications":[
        {"name":"Acetyl", "unimod_accession":1, "mass_shift":42.0367, "amino_acid":null, "n_term":true, "fixed":false},
        {"name":"Oxidation", "unimod_accession":35, "mass_shift":15.9994, "amino_acid":"M", "n_term":false, "fixed":false},
        {"name":"Carbamidomethyl", "unimod_accession":4, "mass_shift":57.0513, "amino_acid":"C", "n_term":false, "fixed":true}
    ]
}
```

#### Percolator
```json
"percolator":{
    "trainFDR":0.01,
    "subset-max-train":200000
}
```

#### MaxQuant
```json
"maxquant_to_rescore":{
    "msms_file":"ms2rescore/tests/data/msms_sample.txt",
    "mgf_dir":"data/mgf",
    "modifications_mapping":{
        "ox":"Oxidation",
        "ac":"Acetyl",
        "cm":"Carbamidomethyl"
    },
    "fixed_modifications":{
        "C":"Carbamidomethyl"
    }
}
```
#### PEPREC
- `<mgf file>` is the spectrum file
- `<peprec file>` is the peptide list file

For this, you need a spectrum file in the [MGF format](http://www.matrixscience.com/help/data_file_help.html) and a peptide list in the [PEPREC format](https://github.com/compomics/ms2pip_c/#peprec-file). The file should have an additional column, `Label`, where target PSMs have a value of `1` and decoys of `-1`, and `Protein`, with the protein identifiers the PSM is associated with (this is optional for the execution of the pipeline, but likely important for posterior analysis).

## Output

Several intermediate files are created when the entire pipeline is ran. Their names are all built based on the provided output filename. The most relevant files are the Percolator INput files, `.pin`.

`Percolator` is ran on each `.pin`, using the settings which have been specified in the config file. For details on these settings and on how to run `Percolator` please refer to its [wiki pages](https://github.com/percolator/percolator/wiki).

From each `Percolator` execution, three files are generated:
- `<file>.pout` with the output regarding target PSMs
- `<file>.pout_dec` with the output regarding decoy PSMs
- `<file>.weights` where the internal feature weights used by `Percolator`'s scoring function are stored.


## Unit testing
```
python pytest --pyargs ms2rescore
```
or
```
pytest --pyargs ms2rescore
```