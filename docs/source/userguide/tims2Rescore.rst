.. _timsrescore:

TIMS²Rescore User Guide
=======================

Introduction
------------

The `TIMS²Rescore` tool is a DDA-PASEF adapted version of `ms2rescore` that allows users to perform rescoring of peptide-spectrum matches (PSMs) acquired on Bruker instruments. This guide provides an overview of how to use `timsrescore` in `ms2rescore` effectively.

Installation
------------

Before using `timsrescore`, ensure that you have `ms2rescore` installed on your system. You can install `ms2rescore` using the following command:

.. code-block:: bash

    pip install ms2rescore

Usage
-----

To use `timsrescore`, follow these steps:

1. Prepare your input files:
    - Ensure that you have the necessary input files, including the PSM file spectrum files
    - Make sure that the PSM file format from a supported search engine or a standard format like .mzid(:external+psm_utils:ref:`supported file formats <supported file formats>`.)
    - Spectrum files can directly be given as .d or minitdf files from Bruker instruments or first converted to .mzML format.

2. Run `timsrescore`:
    - Open a terminal or command prompt.
    - Navigate to the directory where your input files are located.
    - Execute the following command:

      .. code-block:: bash

          timsrescore -p <path_to_psm_file> -s <path_to_spectrum_file> -o <path_to_output_file>

      Replace `<path_to_psm_file>`, `<path_to_tims_file>`, and `<path_to_output_file>` with the actual paths to your input and output files.
    - By default timsTOF specific models will be used for predictions. Optionally you can further configure settings through a configuration file.

3. Review the results:
    - Once the `timsrescore` process completes, you will find the rescoring results in the specified output file or if not specified in the same directory as the input files
    - If you want a detailed overview of the performance, you can either give the set `write_report` to `True` in the configuration file, use the `--write_report` option in the command line or run the following command:
  
  .. code-block:: bash

          ms2rescore-report <output_prefix>

    Replace `<output_prefix>` with the actual output prefix of the result files to the output file.

Additional Options
------------------

`ms2rescore` provides additional options to customize the `timsrescore` process. You can explore these options by running the following command:

.. code-block:: bash

    timsrescore --help

Conclusion
----------

In this guide, you learned how to use `timsrescore` in `ms2rescore` to perform rescoring of PSMs using TIMS data. By following the steps outlined above, you can leverage the power of TIMS to enhance the accuracy of your PSM scores.

For more information and advanced usage, refer to the `ms2rescore` documentation.
