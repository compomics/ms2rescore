#############
Configuration
#############


[TODO Rewrite for MS²Rescore docs]

Parsing of modification labels
------------------------------

MS²Rescore uses the `HUPO-PSI standardized ProForma v2 notation
<https://github.com/HUPO-PSI/ProForma/>`_ to represent modified peptides in a string
format. Unfortunately, most PSM file types coming from different proteomics search
engines use a custom modification notation.

| For example, a MaxQuant ``Modified sequence`` would be parsed as follows:
| ``_AM(Oxidation (M))SIVM(Oxidation (M))LSM_`` 🠚 ``AM[Oxidation (M)]]SIVM[Oxidation (M)]]LSM``

However, the label ``Oxidation (M)`` is not a resolvable modification, as it is not
present in any of the supported controlled vocabularies:


For this conversion, ``Oxidation (M)`` needs to be mapped to ``U:Oxidation``, where ``U``
denotes that the `Unimod <http://www.unimod.org/>`_ database is used and ``Oxidation``
denotes the official Unimod name.

To correctly parse the various notations to ProForma, :py:mod:`psm_utils.io` readers
require :py:obj:`modification_definitions` that map each specific search engine
modification label to a valid ProForma label. :py:obj:`modification_definitions` is defined as a :py:obj:`list` of :py:obj:`dict`'s.
Each :py:obj:`dict` should contain the following key-value pairs:

    - ``site``: Amino acids or peptide termini where the modification occurs. Should be
      the IUPAC one-letter code for amino acid residues and `N-term` or `C-term` for
      terminal modifications. Multiple values can be separated with a pipe character
      (``|``). For example, ``M``, ``K|N-term``, or ``S|T|Y``.
    - ``search_engine_label``: Label for the modification as used by the specific search
      engine that produced the PSM file. For example, ``Oxidation (M)``, or ``ox``.
    - ``proforma_label``: ProForma v2 modification label. See the
      `ProForma documentation <https://github.com/HUPO-PSI/ProForma/>`_ for more info.


Accepted ProForma modification labels in :py:mod:`psm_utils` are, in order of
preference:

=================== ========================= =======================
 Type                Long format example       Short format example
=================== ========================= =======================
 Unimod accession    UNIMOD:21                 U:21
 Unimod name         UNIMOD:Phospho            U:Phospho
 PSI-MOD accession   MOD:00046                 M:00046
 PSI-MOD name        MOD:O-phospho-L-serine    M:O-phospho-L-serine
 Formula             Formula:HO3P              /
 Mass shift          +79.96633052075           /
=================== ========================= =======================

If a modification is not defined in any of the supported controlled vocabularies,
preferably provide the formula instead of a mass shift, as the mass shift can always
be calculated from the formula, but not vice-versa.

For example:

.. code-block:: python

    modification_definitions = [
        {
            "site": "K|N-term",
            "search_engine_label": "TMT6",
            "proforma_label": "UNIMOD:TMT6plex"
        },
        {
            "site": "S|T|Y",
            "search_engine_label": "Phospho",
            "proforma_label": "UNIMOD:21"
        },
        {
            "site": "M",
            "search_engine_label": "custom_modification_1",
            "proforma_label": "Formula:H1O2"
        }
    ]
