**********************
Command line interface
**********************

Run MS²Rescore
==============

.. argparse::
   :module: ms2rescore.__main__
   :func: _argument_parser
   :prog: ms2rescore


Other commands
==============

Generate HTML report
--------------------
Generate a report from MS²Rescore result file(s):

.. code-block:: console

    python -m ms2rescore.report [OPTIONS] OUTPUT_PREFIX



Start graphical user interface
------------------------------
Start the graphical user interface. For more info, see :ref:`Graphical user interface`.

.. code-block:: console

    python -m ms2rescore.gui
