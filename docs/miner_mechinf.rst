The miner-mechinf tool
======================

This utility computes the mechanistic inference.

You can see the tool's available options when you enter ``miner-mechinf -h``
at the command prompt:

.. highlight:: none

::

    usage: miner-mechinf [-h] [-mc MINCORR]
                         expfile mapfile revclusters datadir outdir

    miner-mechinf - MINER compute mechanistic inference

    positional arguments:
      expfile               input matrix
      mapfile               identifier mapping file
      revclusters           revised clusters JSON file
      datadir               data directory
      outdir                output directory

    optional arguments:
      -h, --help            show this help message and exit
      -mc MINCORR, --mincorr MINCORR
                            minimum correlation


Parameters in detail
--------------------

``miner-mechinf`` expects at least these 5 arguments:

  * **expfile:** The gene expression file a matrix in csv format.
  * **mapfile:** The gene identifier map file.
  * **revclusters:** The path to the revised clusters JSON file
  * **datadir:** The path to the data directory
  * **outdir:** The path to the output directory

In addition, you can specify the following optional arguments:

  * ``--mincorr`` or ``--mc``: the minimum correlation value.
