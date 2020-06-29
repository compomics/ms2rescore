<img src="https://github.com/compomics/ms2rescore/raw/dev/img/ms2rescore_logo.png" width="150" height="150" />
<br/><br/>

![Python](https://img.shields.io/badge/python-%3E3.6-blue?style=flat-square)
[![GitHub release](https://img.shields.io/github/release-pre/compomics/ms2rescore.svg?style=flat-square)](https://github.com/compomics/ms2rescore/releases)
[![Build Status](https://img.shields.io/travis/compomics/ms2rescore?style=flat-square)](https://travis-ci.org/compomics/ms2rescore)
[![GitHub issues](https://img.shields.io/github/issues/compomics/ms2rescore?style=flat-square)](https://github.com/compomics/ms2rescore/issues)
[![GitHub](https://img.shields.io/github/license/compomics/ms2rescore.svg?style=flat-square)](https://www.apache.org/licenses/LICENSE-2.0)

Sensitive PSM rescoring with predicted spectra using
[MS²PIP](https://github.com/compomics/ms2pip_c),
[DeepLC](https://github.com/compomics/deeplc), and
[Percolator](https://github.com/percolator/percolator/).

---

- [About MS²ReScore](#about-msrescore)
- [Installation](#installation)
- [Usage](#usage)
  - [Command line interface](#command-line-interface)
  - [Configuration file](#configuration-file)
  - [Output](#output)

---

## About MS²ReScore
Sensitive PSM rescoring with predicted spectra using
[MS²PIP](https://github.com/compomics/ms2pip_c),
[DeepLC](https://github.com/compomics/deeplc), and
[Percolator](https://github.com/percolator/percolator/). This allows more
peptide identifications at a lower false discovery rate.

MS²ReScore takes identifications from multiple search engines:
- [MaxQuant](https://www.maxquant.org/): Start from `msms.txt` identification
  file and directory with `.mgf` files. (Be sure to export without FDR
  filtering!)
- [MSGFPlus](https://omics.pnl.gov/software/ms-gf): Start with an `.mzid`
  identifications file and corresponding `.mgf`.
- [X!Tandem](https://www.thegpm.org/tandem/): Start with an X!Tandem `.xml`
  identifications file and corresponding `.mgf`.
- [PeptideShaker](http://compomics.github.io/projects/peptide-shaker): Start with a
  PeptideShaker Extended PSM Report and corresponding `.mgf` file.

If you use MS²ReScore for your research, please cite the following article:
> **Accurate peptide fragmentation predictions allow data driven approaches to replace and improve upon proteomics search engine scoring functions.** Ana S C Silva, Robbin Bouwmeester, Lennart Martens, and Sven Degroeve. _Bioinformatics_ (2019) [doi:10.1093/bioinformatics/btz383](https://doi.org/10.1093/bioinformatics/btz383)

To replicate the experiments described in this article, check out the [pub branch](https://github.com/compomics/ms2rescore/tree/pub) of this repository.

## Installation
MS²ReScore requires:
- Python >3.6 on Linux (or on [WSL](https://docs.microsoft.com/en-us/windows/wsl)
on Windows 10)
- If the option `run_percolator` is set to `True`, [Percolator](https://github.com/percolator/percolator/) needs to be callable
with the `percolator` command (tested with [version 3.02.1](https://github.com/percolator/percolator/releases/tag/rel-3-02-01))
- For some pipelines, the Percolator converters, such as tandem2pin need to be
  installed.

Clone or download this repository and install:
```
pip install .
```

## Usage
### Command line interface
Run MS²ReScore as follows:
```
usage: ms2rescore [-h] [-m FILE] [-c FILE] [-o FILE] [-l LEVEL]
                  identification_file

MS²ReScore: Sensitive PSM rescoring with predicted MS² peak intensities.

positional arguments:
  identification_file  Path to identification file (mzid, msms.txt, tandem
                       xml)

optional arguments:
  -h, --help           show this help message and exit
  -m FILE              Path to MGF file (default: derived from identifications
                       file). Not applicable to MaxQuant pipeline.
  -c FILE              Path to JSON MS²ReScore configuration file. See
                       README.md for more info. (default: config.json)
  -o FILE              Name for output files (default: derive from
                       identification file
  -l LEVEL             Logging level (default: `info`)
  ```

### Configuration file
It is important to configure MS²ReScore to your use case with the config
file. The config file is written in JSON. Example files for each pipeline are
provided in the GitHub repository.

The config file contains three top level categories (`general`, `ms2pip` and
`percolator`) and an additional category for each pipeline (e.g. `maxquant`).

#### General
- `pipeline`: pipeline to use (currently `MaxQuant`, `MSGFPlus`, or `XTandem`)
- `feature_sets`: list with feature sets to use for rescoring. Options are:
    - `all` = both search engine features and MS²PIP features (recommended)
    - `ms2pip` = only MS²PIP features
    - `searchengine` = only search engine features (classic Percolator)
- `run_percolator`: Whether or not to call Percolator from the MS²ReScore
pipeline. If false, the end result is a Percolator PIN file.
- `keep_tmp_files`: Keep temporary files or not (e.g. MS²PIP output). These
files can be used for a more in-depth data-analysis. If set to true, only the
PIN files and the Percolator output are kept.
- `show_progress_bar`: Whether or not to display a tqdm progress bar.
- `num_cpu`: Number of CPU cores to use. -1 to use all available.

For example:
```json
"general":{
  "pipeline":"MSGFPlus",
  "feature_sets":["all", "ms2pip", "searchengine"],
  "run_percolator":true,
  "keep_tmp_files":false,
  "show_progress_bar":true,
  "num_cpu":8
}
```

#### MS2PIP
These settings are passed to MS²PIP (see [github.com/compomics/ms2pip_c](https://github.com/compomics/ms2pip_c) for more info) and are required for some pipelines to parse
the modification labels correctly.
- `model`: MS²PIP model to use (e.g. `HCD`, see [MS²PIP models](http://compomics.github.io/projects/ms2pip_c#specialized-prediction-models) for more info)
- `frag_error`: MS² mass error tolerance in Da
- `modifications`:
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
    "model":"HCD",
    "frag_error":0.02,
    "modifications":[
        {"name":"Acetyl", "unimod_accession":1, "mass_shift":42.0367, "amino_acid":null, "n_term":true},
        {"name":"Oxidation", "unimod_accession":35, "mass_shift":15.9994, "amino_acid":"M", "n_term":false},
        {"name":"Carbamidomethyl", "unimod_accession":4, "mass_shift":57.0513, "amino_acid":"C", "n_term":false}
    ]
}
```

#### Percolator
This section is only required if `run_percolator` is set to `true` and can contain
command line options, which are directly passed to Percolator (see the
[Percolator wiki](https://github.com/percolator/percolator/wiki/Command-line-options)
for more info). For example:

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
as this ensures correct parsing of the spectrum titles.

:warning: **Make sure to run MaxQuant without FDR filtering (set to 1)!**

Tested with MaxQuant v1.6.2.3.
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
The MSGFPlus pipeline starts from an MSGFPlus `mzid` file. In this pipeline,
next to `percolator`, the `msgf2pin` command also needs to be callable. No extra
section in the config file is required.

:warning: **Be sure to run MSGFPlus as a concatenated target-decoy search, with the `-addFeatures 1` flag.**

#### X!Tandem
The X!Tandem pipeline starts with the identifications `.xml` file.  No extra
section in the config file is required.

### Output
Several intermediate files are created when the entire pipeline is run. Their
names are all built based on the provided output filename. Depending on the
pipeline, the `keep_tmp_files` setting and whether or not Percolator is run, the
following output files can be expected:

For each feature set (`all`, `ms2pip` and/or `searchengine`):
- `<file>.pin` Percolator IN file
- `<file>.pout` Percolator OUT file with target PSMs
- `<file>.pout_dec` Percolator OUT file with decoy PSMs
- `<file>.weights` Internal feature weights used by Percolator's scoring
function.


## Development
Conda environment:
```
conda create -y --name ms2rescore-dev -c bioconda -c defaults -c conda-forge python=3.7 percolator deeplc pyteomics numpy pandas scipy scikit-learn tqdm xmltodict pylint pydocstyle pytest seaborn matplotlib biopython numexpr tomlkit
```

Activatate the environment:
```
conda activate ms2rescore-dev
```

Install MS²PIP, which is not in any conda channel:
```
pip install ms2pip
```

Install MS²ReScore in "editable" mode:
```
pip install --editable .
```
