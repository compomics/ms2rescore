#############
Configuration
#############


Introduction
============
MS²Rescore can be configured through the command line interface (CLI), the graphical user interface
(GUI), or a JSON/TOML configuration file. The configuration file can be used to set options that
are not available in the CLI or GUI, or to set default values for options that are available in the
CLI or GUI.

If no configuration file is passed, or some options are not configured, the
`default values <https://github.com/compomics/ms2rescore/blob/main/ms2rescore/package_data/config_default.json>`_
for these settings will be used. Options passed from the CLI and the GUI will override
the configuration file. The full configuration is validated against a
`JSON Schema <https://github.com/compomics/ms2rescore/blob/main/ms2rescore/package_data/config_schema.json>`_.
A full example configuration file can be found in
`ms2rescore/package_data/config_default.json <https://github.com/compomics/ms2rescore/blob/main/ms2rescore/package_data/config_default.json>`_.
An overview of all options can be found below.


Configuring input files
=======================

In the configuration file, input files can be specified as follows:

.. tab:: JSON

  .. code-block:: json

    "psm_file": "path/to/psms.tsv",
    "psm_file_type": "infer",
    "spectrum_path": "path/to/spectra.mgf"

.. tab:: TOML

  .. code-block:: toml

    psm_file = "path/to/psms.tsv"
    psm_file_type = "infer"
    spectrum_path = "path/to/spectra.mgf"

See :ref:`Input files` for more information.


Parsing modification labels
===========================

MS²Rescore uses the `HUPO-PSI standardized ProForma v2 notation
<https://github.com/HUPO-PSI/ProForma/>`_ to represent modified peptides in a string
format. Unfortunately, most PSM file types coming from different proteomics search
engines use a custom modification notation.

For example, a MaxQuant ``Modified sequence`` would be parsed as follows: ``_AM(ox)SIVMLSM_`` 🠚
``AM[ox]SIVMLSM``. However, the label ``ox`` is not a resolvable modification, as it is not
present in any of the supported controlled vocabularies. Therefore, ``ox`` needs to be mapped to
``U:Oxidation``, where ``U`` denotes that the `Unimod <http://www.unimod.org/>`_ database is used
and ``Oxidation`` denotes the official Unimod name.

To correctly parse the various notations to ProForma, :py:mod:`ms2rescore` requires a configuration
:py:obj:`modification_mapping` which maps each specific search engine modification label to a valid
ProForma label.

Accepted ProForma modification labels in :py:mod:`psm_utils` (and by extension in
:py:mod:`ms2rescore`) are, in order of preference:

=================== ========================= =======================
 Type                Long format example       Short format example
=================== ========================= =======================
 PSI-MOD accession   MOD:00046                 M:00046
 PSI-MOD name        MOD:O-phospho-L-serine    M:O-phospho-L-serine
 Unimod accession    UNIMOD:21                 U:21
 Unimod name         UNIMOD:Phospho            U:Phospho
 Formula             Formula:HO3P              /
 Mass shift          +79.96633052075           /
=================== ========================= =======================

If a modification is not defined in any of the supported controlled vocabularies,
preferably provide the formula instead of a mass shift, as the mass shift can always
be calculated from the formula, but not vice-versa, and some feature generators (such as DeepLC)
require the modification formula.

And example of the :py:obj:`modification_mapping` could be:

.. tab:: JSON

  .. code-block:: json

    "modification_mapping": {
      "gl": "U:Gln->pyro-Glu",
      "ox": "U:Oxidation",
      "ac": "U:Acetylation",
      "de": "U:Deamidation"
    }

.. tab:: TOML

  .. code-block:: toml

    [ms2rescore.modification_mapping]
    "gl" = "Gln->pyro-Glu"
    "ox" = "Oxidation"
    "ac" = "Acetylation"
    "de" = "Deamidation"

.. tab:: GUI

  .. figure:: ../_static/img/gui-modification-mapping.png
    :width: 500px
    :alt: modification mapping configuration in GUI


Adding fixed modifications
==========================

Some search engines, such as MaxQuant, do not report fixed modifications that were part of the
search. To correctly rescore PSMs, fixed modifications that are not reported in the PSM file must
be configured separately. For instance:

.. tab:: JSON

  .. code-block:: json

    "fixed_modifications": {
      "C": "U:Carbamidomethyl"
    }

.. tab:: TOML

    .. code-block:: toml

      [ms2rescore.fixed_modifications]
      "Carbamidomethyl" = ["C"]

.. tab:: GUI

  .. figure:: ../_static/img/gui-fixed-modifications.png
    :width: 500px
    :alt: fixed modifications configuration in GUI



Mapping PSMs to spectra
=======================

Essential for MS²Rescore to function correctly is linking the search engine PSMs to the original
spectra. As spectrum file converters and search engines often modify spectrum titles, two options
are available to map PSMs to spectra: ``spectrum_id_pattern`` and ``psm_id_pattern``. Through these
two options, regular expression patterns can be defined that extract the same spectrum identifier
from the spectrum file and from the PSM file, respectively.

For example, if the spectrum file contains the following identifier in the MGF title field:

.. code-block:: text

  mzspec=20161213_NGHF_DBJ_SA_Exp3A_HeLa_1ug_7min_15000_02.raw: controllerType=0 controllerNumber=1 scan=2

and the PSM file contains the following identifier in the ``spectrum_id`` field:

.. code-block:: text

  20161213_NGHF_DBJ_SA_Exp3A_HeLa_1ug_7min_15000_02.raw.2.2

then the following patterns can be used to extract ``2`` from both identifiers:

.. tab:: JSON

  .. code-block:: json

    "spectrum_id_pattern": ".*scan=(\\d+)$",
    "psm_id_pattern": ".*\\..*\\.(.*)"

.. tab:: TOML

  .. code-block:: toml

      spectrum_id_pattern = '.*scan=(\d+)$'
      psm_id_pattern = ".*\..*\.(.*)"


Both options should match the entire string and require a single capture group (denoted by the
parentheses) to mark the section of the match that should be extracted.

.. warning::
  Regular expression patterns often contain special characters that need to be escaped. For example,
  the ``\`` should be escaped with an additional ``\`` in JSON, as is shown above. In TOML files,
  the full regex can be wrapped in single quotes to avoid excaping.

.. note::
  Find out more about regular expression patterns and try them on
  `regex101.com <https://regex101.com/>`_. You can try out the above examples at
  https://regex101.com/r/VhBJRM/1 and https://regex101.com/r/JkT79a/1.


Selecting decoy PSMs
====================

Usually, PSMs are already marked as target or decoy in the PSM file. When this is not the case,
it can usually be derived from the protein name. For example, if the protein name contains the
prefix ``DECOY_``, the PSM is a decoy PSM. The following option can be used to define a regular
expression pattern that extracts the decoy status from the protein name:

.. tab:: JSON

  .. code-block:: json

    "decoy_pattern": "DECOY_"

.. tab:: TOML

    .. code-block:: toml

      decoy_pattern = "DECOY_"



All configuration options
=========================

.. include:: ../config_schema.md
   :parser: myst_parser.sphinx_
