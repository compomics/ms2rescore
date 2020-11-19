# MS²ReScore configuration

*MS²ReScore JSON configuration file.*

## Properties

- **`general`** *(object)*: General MS²ReScore settings. Cannot contain additional properties.
  - **`pipeline`** *(string)*: Pipeline to use, depending on input format. Must be one of: `['infer', 'pin', 'tandem', 'maxquant', 'msgfplus', 'peptideshaker']`. Default: `infer`.
  - **`feature_sets`** *(array)*: Feature sets for which to generate PIN files and optionally run Percolator. Default: `['all']`.
    - **Items** *(string)*: Must be one of: `['all', 'ms2pip_rt', 'searchengine', 'rt', 'ms2pip']`.
  - **`id_decoy_pattern`**: Pattern used to identify the decoy PSMs in identification file. Passed to `--pattern` option of Percolator converters. Default: `None`.
  - **`run_percolator`** *(boolean)*: Run Percolator within MS²ReScore. Default: `False`.
  - **`num_cpu`** *(number)*: Number of parallel processes to use; -1 for all available. Minimum: `-1`. Default: `-1`.
  - **`config_file`**: Path to configuration file.
  - **`identification_file`** *(string)*: Path to identification file.
  - **`mgf_path`**: Path to MGF file or directory with MGF files.
  - **`tmp_path`**: Path to directory to place temporary files.
  - **`output_filename`**: Path and root name for output files.
  - **`log_level`** *(string)*: Logging level. Must be one of: `['debug', 'info', 'warning', 'error', 'critical']`.
- **`ms2pip`** *(object)*: MS²PIP settings. Cannot contain additional properties.
  - **`model`** *(string)*: MS²PIP model to use (see MS²PIP documentation). Default: `HCD`.
  - **`frag_error`** *(number)*: MS2 error tolerance in Da. Minimum: `0`. Default: `0.02`.
  - **`modifications`** *(array)*: Array of peptide mass modifications.
    - **Items**: Refer to *#/definitions/modifications*.
- **`percolator`** *(object)*: Command line options directly passed to Percolator (see the Percolator wiki).
- **`maxquant_to_rescore`** *(object)*: Settings specific to the MaxQuant pipeline. Cannot contain additional properties.
  - **`mgf_dir`** *(string)*: Path to directory with MGF files.
  - **`modification_mapping`** *(object)*: Mapping of MaxQuant modification labels to modifications names for MS²PIP. Default: `{}`.
  - **`fixed_modifications`** *(object)*: Mapping of amino acids with fixed modifications to the modification name. Default: `{}`.
## Definitions

- **`modifications`** *(object)*: Peptide mass modifications, per amino acid. Cannot contain additional properties.
  - **`name`** *(string)*: Unique name for modification.
  - **`unimod_accession`** *(number)*: Unimod accession of modification.
  - **`mass_shift`** *(number)*: Mono-isotopic mass shift.
  - **`amino_acid`**: Amino acid one-letter code, or null if amino acid-agnostic (e.g. N-term acetylation).
  - **`n_term`** *(boolean)*: Modification is N-terminal.
