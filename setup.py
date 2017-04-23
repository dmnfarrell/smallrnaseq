from setuptools import setup
import sys,os

setup(
    name = 'smallrnaseq',
    version = '0.2.0',
    description = 'Python package for micro rna-seq routines ',
    url='https://github.com/dmnfarrell/smallrnaseq',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien[at]gmail.com',
    packages = ['smallrnaseq'],
    package_data={'smallrnaseq': ['data/*','*.R']},
    install_requires=['numpy>=1.10',
                      'pandas>=0.17',
                      'HTSeq>=0.6',
                      'pybedtools>=0.7.9',
                      'seaborn>=0.7',
                      'scikit-learn>=0.18',
                      'bx-python>=0.5',
                      'forgi>=0.4'],
    entry_points = {
        'console_scripts': [
            'smallrnaseq=smallrnaseq.app:main',
            'mirdeep2=smallrnaseq.mirdeep2:main']
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
    keywords = ['rna','sequencing','mirdeep2','biology','scientific'],
)
