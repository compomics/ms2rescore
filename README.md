# MS²ReScore
[![GitHub release](https://img.shields.io/github/release-pre/compomics/ms2rescore.svg)](https://github.com/compomics/ms2rescore/releases)
[![Build Status](https://travis-ci.org/compomics/ms2rescore.svg?branch=dev)](https://travis-ci.org/compomics/ms2rescore)
[![GitHub](https://img.shields.io/github/license/compomics/ms2pip_c.svg)](https://www.apache.org/licenses/LICENSE-2.0)

Use features calculated from comparing experimental spectra with computationally
generated spectra (see [MS²PIP](https://github.com/compomics/ms2pip_c)) to
rescore peptide identifications using
[Percolator](https://github.com/percolator/percolator/).

On this branch, multiple pipelines can be run, depending on your input format:
- [MaxQuant](https://www.maxquant.org/): Start from `msms.txt` identification
file and directory with `.mgf` files. Be sure to export without FDR filtering!
- [MSGFPlus](https://omics.pnl.gov/software/ms-gf): Start with `.mgf` and `.fasta`
file, or from `.mzid` identifications file

## Prerequisites
- Python 3.7 on Linux
- If the option `run_percolator` is set to True, [Percolator](https://github.com/percolator/percolator/) needs to be callable
with the `percolator` command (tested with version 3.02.1)

## Installation
Clone this repository. This includes the submodules
[mapper](https://github.com/anasilviacs/mapper/tree/0ee46adcbb20a118a8274908255cc8b3f95a51db)
and [MS2PIP](https://github.com/anasilviacs/ms2pip_c/tree/6f037dc2d0797cd25061aaed8091d625123971e1).
As such, the repository should be called as:

```
git clone --recurse-submodules git://github.com/compomics/ms2rescore.git
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

### Configuration
The main argument is a path to the config file. This JSON file should contain
all required information for MS²ReScore to be run properly. Example files for
each pipeline are provided in the GitHub repository.

The config file contains three main top level keys (`general`, `ms2pip` and
`percolator`) and a key for each pipeline (e.g. `maxquant`). 

#### General
- `pipeline`: pipeline to use (currently `MaxQuant` or `MSGFPlus`)
- `feature_sets`: list with feature sets to use for rescoring. Options are:
    - `all` = both search engine features and MS²PIP features
    - `ms2pip` = only MS²PIP features
    - `searchengine` = only search engine features (classic Percolator)
- `run_percolator`: Whether or not to call Percolator from the MS²ReScore
pipeline. If false, the end result is a Percolator PIN file.
- `keep_tmp_files`: Keep temporary files or not (e.g. MS²PIP output). These
files can be used for a more in-depth data-analysis. If set to true, only the
PIN files and the Percolator output are kept.
- `show_progress_bar`: Whether or not to display a tqdm progress bar.
- `num_cpu`: Number of CPU cores to use.

For example:
```json
"general":{
  "pipeline":"MSGFPlus",
  "feature_sets":["all", "ms2pip", "searchengine"],
  "run_percolator":true,
  "keep_tmp_files":false,
  "show_progress_bar":true,
  "num_cpu":24
}
```

#### MS2PIP
These settings are passed to MS²PIP (see [github.com/compomics/ms2pip_c](https://github.com/compomics/ms2pip_c) for more info).
- `dir`: Path to MS²PIP folder, containing ms2pipC.py
- `model`: MS²PIP model to use (e.g. `HCD`, see [MS²PIP models](https://github.com/compomics/ms2pip_c#mspip-models) for more info)
- `frag_error`: MS² mass error tolerance in Da
- `Modifications`: 
    - `name`: as used in e.g. MaxQuant `modifications_mapping` (see below)
    - `unimod_accession`: Required for parsing MSGFPlus output (see [unimod.org](http://www.unimod.org/) for correct accession numbers)
    - `mass_shift`: Mono-isotopic mass shift
    - `amino_acid`: Amino acid on which the modification occurs, or `null` if
    e.g. N-terminal modification
    - `n_term`: Whether or not the modification is N-terminal (C-terminal
    modifications are not yet supported)

For example
```json
"ms2pip":{
    "dir":"ms2pip_c",
    "frag":"HCD",
    "frag_error":0.02,
    "modifications":[
        {"name":"Acetyl", "unimod_accession":1, "mass_shift":42.0367, "amino_acid":null, "n_term":true},
        {"name":"Oxidation", "unimod_accession":35, "mass_shift":15.9994, "amino_acid":"M", "n_term":false},
        {"name":"Carbamidomethyl", "unimod_accession":4, "mass_shift":57.0513, "amino_acid":"C", "n_term":false}
    ]
}
```

#### Percolator
Command line options directly passed to Percolator (see the [Percolator wiki](https://github.com/percolator/percolator/wiki/Command-line-options) for more info). For
example:

```json
"percolator":{
    "trainFDR": 0.01
}
```
In this case, `--trainFDR 0.01` is passed to Percolator.

#### MaxQuant
The MaxQuant pipeline starts with an `msms.txt` file and a directory containing
MGF files. To convert Raw files to MGF, please use the
[CompOmics ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser/),
as this ensures correct parsing of the spectrum titles. Make sure to run
MaxQuant without FDR filtering (set to 1)!  
Tested with MaxQuant
v1.6.2.3.
- `msms_file`: Path to msms.txt file.
- `mgf_dir`: Path to directory containing MGF files.
- `modifications_mapping`: Maps MaxQuant output to MS²PIP modifications list.
Keys must contain MaxQuant's two-letter modification codes and values must match
one of the modifications listed in the MS²PIP configuration (see
[MS2PIP config](#MS2PIP)).
- `fixed_modifications`: Must list all modifications set as fixed during the
MaxQuant search (as this is not denoted in the msms.txt file). Keys refer to the
amino acid, values to the modification name used in the MS²PIP configuration.

For example:
```json
"maxquant_to_rescore":{
  "msms_file":"examples/id/msms.txt",
  "mgf_dir":"examples/mgf",
  "modifications_mapping":{
    "ox":"Oxidation",
    "ac":"Acetyl",
    "cm":"Carbamidomethyl",
    "de":"Deamidated",
    "gl":"Gln->pyro-Glu"
  },
  "fixed_modifications":{
    "C":"Carbamidomethyl"
}
```

#### MSGFPlus
The MSGFPlus pipeline can either include the search (start from a fasta and MGF
file) or start from an MSGFPlus mzid file. In the latter case, be sure to add
`-addFeatures 1` when running MSGFPlus, as this is required for Percolator.
In this pipeline, next to `percolator`, the `msgf2pin` command also needs to be
callable.

- `run_search`: Whether or not to run the search in this pipeline. If true, the
path to the MSGFPlus jar file is required. If false, the path to the mzid file
is required.
- `mgf_file`: Path to the MGF file.
- `mzid_file`: Path to the mzid file (only required if `run_search` is false)
- `search_params`:
    - `jar_file`: Path to MSGFPlus jar file.
    - `fasta_file`: Path to fasta search database. Does not need to include
    decoy sequences; these are added automatically.
    - `frag`: Fragmentation method (e.g. `HCD`).
    - `path_to_modsfile`: Path to MSGFPlus modifications config file.
    - `min_length`, `min_charge`, `max_charge` and `ms1_tolerance`: respective
    search settings for MSGFPlus.

For example:
```json
"msgfplus": {
  "run_search": false,
  "mgf_file": "examples/mgf/20161213_NGHF_DBJ_SA_Exp3A_HeLa_1ug_7min_15000_02.mgf",
  "mzid_file": "examples/id/msgfplus.mzid",
  "search_params":{
    "jar_file": "MSGFPlus.jar",
    "fasta_file": "examples/fasta/uniprot-proteome-human-contaminants.fasta",
    "frag": "HCD",
    "path_to_modsfile": "examples/parameters/msgfplus_modifications.txt",
    "min_length": 8,
    "min_charge": 2,
    "max_charge ": 4,
    "ms1_tolerance": "10ppm"
}
```

## Output

Several intermediate files are created when the entire pipeline is run. Their
names are all built based on the provided output filename. Depending on the
keep_tmp_files setting and whether or not Percolator is run, the following
output files can be expected:

For each feature set (`all`, `ms2pip` and/or `searchengine`):
- `<file>.pin` Percolator IN file
- `<file>.pout` Percolator OUT file with target PSMs
- `<file>.pout_dec` Percolator OUT file with decoy PSMs
- `<file>.weights` Internal feature weights used by Percolator's scoring
function.
