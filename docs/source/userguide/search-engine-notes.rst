#################################
Notes for specific search engines
#################################

MSGFPlus
========

- Run MSGFPlus in a concatenated target-decoy search, with the ``-addFeatures 1`` flag.


MaxQuant
========

- Run MaxQuant without FDR filtering (set to 1)
- Make sure to correctly configure both ``modification_mapping`` and ``fixed_modifications``.
  See :ref:`Parsing modification labels` for more information.
