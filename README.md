<img src="https://github.com/compomics/ms2rescore/raw/dev/img/ms2rescore_logo.png" width="150" height="150" />
<br/><br/>

[![GitHub release](https://img.shields.io/github/release-pre/compomics/ms2rescore.svg?style=flat-square)](https://github.com/compomics/ms2rescore/releases)
[![PyPI](https://flat.badgen.net/pypi/v/ms2rescore)](https://pypi.org/project/ms2rescore/)
[![GitHub Workflow Status](https://flat.badgen.net/github/checks/compomics/ms2rescore/master)](https://github.com/compomics/ms2rescore/actions/)
[![GitHub issues](https://img.shields.io/github/issues/compomics/ms2rescore?style=flat-square)](https://github.com/compomics/ms2rescore/issues)
[![GitHub](https://img.shields.io/github/license/compomics/ms2rescore.svg?style=flat-square)](https://www.apache.org/licenses/LICENSE-2.0)
[![Last commit](https://flat.badgen.net/github/last-commit/compomics/ms2rescore)](https://github.com/compomics/ms2rescore/commits/)
[![Twitter](https://flat.badgen.net/twitter/follow/compomics?icon=twitter)](https://twitter.com/compomics)


Sensitive peptide identification rescoring with predicted spectra using
[MS²PIP](https://github.com/compomics/ms2pip_c),
[DeepLC](https://github.com/compomics/deeplc), and
[Percolator](https://github.com/percolator/percolator/).

---

- [About MS²ReScore](#about-msrescore)
- [Installation](#installation)
- [Usage](#usage)
  - [Command line interface](#command-line-interface)
  - [Configuration file](#configuration-file)
  - [Notes for specific search engines](#notes-for-specific-search-engines)
  - [Output](#output)

---

## About MS²ReScore

MS²ReScore performs sensitive peptide identification rescoring with predicted
spectra using [MS²PIP](https://github.com/compomics/ms2pip_c),
[DeepLC](https://github.com/compomics/deeplc), and
[Percolator](https://github.com/percolator/percolator/). This results in more confident
peptide identifications, which allows you to get **more peptide IDs** at the same false
discovery rate (FDR) threshold, or to set a **more stringent FDR threshold** while still
retaining a similar number of peptide IDs. MS²ReScore is **ideal for challenging
proteomics identification workflows**, such as proteogenomics, metaproteomics, or
immunopeptidomics.

MS²ReScore uses identifications from a
[Percolator IN (PIN) file](https://github.com/percolator/percolator/wiki/Interface#tab-delimited-file-format),
or from the output of one of these search engines:

- [MaxQuant](https://www.maxquant.org/): Start from `msms.txt` identification
  file and directory with `.mgf` files. (Be sure to export without FDR
  filtering!)
- [MSGFPlus](https://omics.pnl.gov/software/ms-gf): Start with an `.mzid`
  identification file and corresponding `.mgf`.
- [X!Tandem](https://www.thegpm.org/tandem/): Start with an X!Tandem `.xml`
  identification file and corresponding `.mgf`.
- [PeptideShaker](http://compomics.github.io/projects/peptide-shaker): Start with a
  PeptideShaker Extended PSM Report and corresponding `.mgf` file.

If you use MS²ReScore for your research, please cite the following article:

> **Accurate peptide fragmentation predictions allow data driven approaches to replace
and improve upon proteomics search engine scoring functions.** Ana S C Silva, Robbin
Bouwmeester, Lennart Martens, and Sven Degroeve. _Bioinformatics_ (2019)
[doi:10.1093/bioinformatics/btz383](https://doi.org/10.1093/bioinformatics/btz383)

To replicate the experiments described in this article, check out the
[pub branch](https://github.com/compomics/ms2rescore/tree/pub) of this repository.

---

## Installation

[![install pip](https://flat.badgen.net/badge/install%20with/pip/green)](https://pypi.org/project/ms2rescore/)

MS²ReScore requires:
- Python 3.7 or higher on Linux, macOS, or
[Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl)
- If the option `run_percolator` is set to `True`,
[Percolator](https://github.com/percolator/percolator/) needs to be callable with the
`percolator` command (tested with
[version 3.02.1](https://github.com/percolator/percolator/releases/tag/rel-3-02-01))
- Some pipelines require the Percolator converters, such as `tandem2pin`, as well. These
are usually installed alongside Percolator.

Minimal installation:
```sh
pip install ms2rescore
```

Recommended installation, including DeepLC for retention time prediction:
```sh
pip install ms2rescore[deeplc]
```

We recommend using a [venv](https://docs.python.org/3/library/venv.html) or
[conda](https://docs.conda.io/en/latest/) virtual environment.


---

## Usage

### Command line interface

Run MS²ReScore from the command line as follows:

```
ms2rescore -c <path-to-config-file> -m <path-to-mgf> <path-to-identification-file>
```

Run `ms2rescore --help` to see all command line options.

### Configuration file

MS²ReScore can be further configured through a **JSON configuration file**. A correct
configuration is required to, for example, correctly parse the peptide modifications
from the search engine output. If no configuration file is passed, or some options are
not configured, the
[default values](https://github.com/compomics/ms2rescore/blob/master/ms2rescore/package_data/config_default.json)
for these settings will be used. Options passed from the command line will override
the configuration file. The full configuration is validated against a
[JSON Schema](https://github.com/compomics/ms2rescore/blob/master/ms2rescore/package_data/config_schema.json).

A full example configuration file can be found in
[ms2rescore/package_data/config_default.json](https://github.com/compomics/ms2rescore/blob/master/ms2rescore/package_data/config_default.json).

The config file contains three top level categories (`general`, `ms2pip` and
`percolator`) and an additional categories for specific search engines
(e.g. `maxquant`). The most important options in `general` are:
- **`pipeline`** *(string)*: Pipeline to use, depending on input format. Must be one of:
`['infer', 'pin', 'tandem', 'maxquant', 'msgfplus', 'peptideshaker']`. Default: `infer`.
- **`feature_sets`** *(array)*: Feature sets for which to generate PIN files and
optionally run Percolator. Default: `['all']`.
  - **Items** *(string)*: Must be one of:
  `['all', 'ms2pip_rt', 'searchengine', 'rt', 'ms2pip']`.

An overview of all options can be found in [configuration.md](https://github.com/compomics/ms2rescore/blob/master/configuration.md)

### Notes for specific search engines

- **MSGFPlus:** Run MSGFPlus in a concatenated target-decoy search, with the
`-addFeatures 1` flag.
- **MaxQuant:**
  - Run MaxQuant without FDR filtering (set to 1)
  - MaxQuant requires additional options in the configuration file:
    - `modification_mapping`: Maps MaxQuant output to MS²PIP modifications list.
Keys must contain MaxQuant's two-letter modification codes and values must match
one of the modifications listed in the MS²PIP configuration (see
[MS2PIP config](#MS2PIP)).
    - `fixed_modifications`: Must list all modifications set as fixed during the
MaxQuant search (as this is not denoted in the msms.txt file). Keys refer to the
amino acid, values to the modification name used in the MS²PIP configuration.

### Output
Several intermediate files are created when the entire pipeline is run. These can be
accessed by specifying the `tmp_dir` option. Depending on whether or not Percolator is
run, the following output files can be expected:

For each feature set (e.g. `all`, `ms2pip`, `searchengine`...):
- `<file>.pin` Percolator IN file
- `<file>.pout` Percolator OUT file with target PSMs
- `<file>.pout_dec` Percolator OUT file with decoy PSMs
- `<file>.weights` Internal feature weights used by Percolator's scoring function.
