###########
Input files
###########

PSM file(s)
===========

The peptide-spectrum match (PSM) file is generally the output from a proteomics search engine.
This file serves as the main input to MS²Rescore. One or multiple PSM files can be provided at
once. Note that merging PSMs from different MS runs could have an impact on the correctness of
the FDR control.

Various PSM file types are supported. The type can be specified with the ``psm_file_type`` option.
Check the list of :py:mod:`psm_utils` tags in the
:external+psm_utils:ref:`supported file formats <supported file formats>` section. Depending on the
file extension, the file type can also be inferred from the file name. In that case,
``psm_file_type`` option can be set to ``infer``.

.. attention::
   As a general rule, MS²Rescore always needs access to **all target and decoy PSMs, without any
   FDR-filtering**. For some search engines, this means that the FDR-filter should be disabled or
   set to 100%.


Spectrum file(s)
================

Spectrum files are required for some feature generators. Both ``mzML`` and ``mgf`` formats are
supported. The ``spectrum_path`` option can be either a single file or a folder. If the
``spectrum_path`` is a folder, MS²Rescore will search for spectrum files in the directory according
to the run names in the PSM file.
