#!/usr/bin/env python
from setuptools import setup, find_packages
import re
from glob import glob


__author__ = "Jason Sahl"
__credits__ = ["Jason Sahl"]
__license__ = "GPL v3"
__version__ = "1.0"
__maintainer__ = "Jason Sahl"
__email__ = "jsahl@tgen.org"
__status__ = "Development"
 
long_description = """WG-FAST, genotype
phylogenetically from a NASP formatted SNP
matrix
"""

data_files = ['ext/standard-raxml/raxmlHPC-PTHREADS-SSE3',
              'ext/standard-raxml/raxmlHPC-SSE3']
data_files += glob("ext/*.jar")
data_files += glob("ext/*.fasta")
data_files += glob("ext/bbmap/*.sh")


setup(name='wgfast',
      version=__version__,
      description='whole genome focused array SNP typing',
      author=__maintainer__,
      author_email=__email__,
      url='https://github.com/jasonsahl/wgfast',
      long_description=long_description,
      packages=find_packages(),
      entry_points={'console_scripts': ['wgfast=wgfast.main:main']},
      include_package_data=True,
      package_data={'wgfast':['command.txt']},
      data_files=[
          ('wgfast/bin', data_files)
      ],
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research']
      )
