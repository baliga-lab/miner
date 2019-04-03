import os
import re
from setuptools import setup

NAME = 'isb_miner'
PACKAGES = ['miner']
DESCRIPTION = 'Mechanistic Inference of Non-edge Relationships'
LICENSE = 'LGPL V3'
URI = 'https://github.com/baliga-lab/miner'
AUTHOR = 'Baliga Lab, Institute for Systems Biology'
VERSION = '1.0.7'

KEYWORDS = ["isb", "miner", "mechanistic", "inference", "network", "gene", "regulatory", "biological"]

# See trove classifiers
# https://testpypi.python.org/pypi?%3Aaction=list_classifiers

CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules"
    ]
INSTALL_REQUIRES = ['numpy', 'scipy', 'pandas', 'sklearn', 'matplotlib', 'seaborn', 'lifelines']

if __name__ == '__main__':
    setup(name=NAME, description=DESCRIPTION,
          license=LICENSE,
          url=URI,
          version=VERSION,
          author=AUTHOR,
          author_email='mwall@systemsbiology.net',
          maintainer=AUTHOR,
          maintainer_email='mwall@systemsbiology.net',
          keywords=KEYWORDS,
          packages=PACKAGES,
          zip_safe=False,
          classifiers=CLASSIFIERS,
          install_requires=INSTALL_REQUIRES,
          scripts=['bin/miner-coexpr', 'bin/miner-mechinf',
                   'bin/miner-enrichment',
                   'bin/miner-bcmembers', 'bin/miner-subtypes',
                   'bin/miner-survival', 'bin/miner-causalinf-pre',
                   'bin/miner-causalinf-post', 'bin/miner-neo'])
