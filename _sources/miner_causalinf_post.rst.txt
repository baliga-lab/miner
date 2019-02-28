The miner-causalinf-post tool
=============================

This utility generates results from the output of the miner-neo utility

You can see the tool's available options when you enter ``miner-causalinf-post -h``
at the command prompt:

.. highlight:: none

::

    usage: miner-causalinf-post [-h]
                                expfile mapfile eigengenes neoresults datadir
                                outdir

    miner-causalinf - MINER compute causal inference

    positional arguments:
      expfile     input matrix
      mapfile     identifier mapping file
      eigengenes  eigengenes file used as NEO input
      neoresults  NEO results directory
      datadir     data directory
      outdir      output directory

    optional arguments:
      -h, --help  show this help message and exit

Parameters in detail
--------------------

``miner-causalinf-post`` expects at least these 6 arguments:

  * **expfile:** The gene expression file a matrix in csv format.
  * **mapfile:** The gene identifier map file.
  * **eigengenes:** The eigengenes file used as input to miner-neo
  * **datadir:** the reult directory used for miner-neo
  * **outdir:** The directory where the results of this tool will be stored in
