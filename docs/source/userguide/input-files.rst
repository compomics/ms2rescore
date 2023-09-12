###########
Input files
###########

PSM file
========

[todo]

As a general rule, MS²Rescore always needs access to all target and decoy PSMs, not
only the FDR-filtered targets.

The ``psm_file_type`` can be one of the :py:mod:`psm_utils` tags as listed in the
:external+psm_utils:ref:`supported file formats <supported file formats>`. Depending on the file
extension, the file type can also be inferred from the file name. In that case, ``psm_file_type``
option can be set to ``infer``.


Spectrum file(s)
================

[todo]

If the ``spectrum_path`` is a directory, MS²Rescore will search for spectrum files in the
directory according to the run names in the PSM file.
