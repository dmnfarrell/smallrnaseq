from setuptools import setup
import sys,os

setup(
    name = 'smallnaseq',
    version = '0.1.0',
    description = 'Python package for micro rna-seq routines ',
    url='https://github.com/dmnfarrell/smallrnaseq',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien[at]gmail.com',
    packages = ['smallrnaseq'],
    package_data={'smallrnaseq': ['data/*']},
    install_requires=['pandas>=0.17',
                      'biopython>=1.5',
                      'HTSeq>0.6'],
    entry_points = {
        'console_scripts': [
            'mirdeep2=smallrnaseq.mirdeep2:main',
            'srnabench=smallrnaseq.srnabench:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 2.7',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: Apache Software License',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research'],
    keywords = ['mirna','sequencing','mirdeep2','biology'],
)
