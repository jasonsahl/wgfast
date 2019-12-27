from setuptools import setup, find_packages

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



setup(name='wgfast',
      version=__version__,
      description='whole genome focused array SNP typing',
      author=__maintainer__,
      author_email=__email__,
      url='https://github.com/jasonsahl/wgfast',
      long_description=long_description,
      packages=find_packages(),
      install_requires=[
          'click'],
      entry_points={'console_scripts': [
            'wgfast=wgfast.main:main',
            'wgfast_prep=wgfast.tools.wgfast_prep:main',
            'wgfastdb=wgfastdb.main:main',
            'wgfast_genbankqc=wgfast.genbankqc.__main__:cli'
            # 'wgfast_extract_tree_names=wgfast.tools.extract_tree_names:main',
            # 'wgfast_subsample_reads_and_place=wgfast.tools.subsample_reads_and_place:main',
            # 'wgfast_subsample_snps_pearson=wgfast.tools.subsample_snps_pearson:main'
            ]},
      include_package_data=True,
      package_data={'wgfast': ['command.txt',
            'bin/*.fasta',
            'bin/*.jar']},
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research']
      )
