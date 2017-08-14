from setuptools import setup
import sys,os

inst_requires = ['pandas>=0.17',
                    'seaborn>=0.7',
                    'scikit-learn>=0.18',
                    'pybedtools>=0.7.9',
                    'pysam>=0.10.0',
                    'HTSeq>=0.6',
                    'bx-python>=0.5',
                    'forgi>=0.4']

major, minor, micro = sys.version_info[:3]
if major == '2':
    inst_requires.append('future')

setup(
    name = 'smallrnaseq',
    version = '0.4.0',
    description = 'Package for short RNA-seq analysis',
    long_description = 'smallrnaseq is a Python package for processing of small RNA seq data.',
    url='https://github.com/dmnfarrell/smallrnaseq',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['smallrnaseq'],
    package_data={'smallrnaseq': ['data/*','*.R']},
    install_requires=inst_requires,
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
