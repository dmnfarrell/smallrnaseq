from setuptools import setup
import sys,os

inst_requires = ['numpy>=1.10',
                    'pandas>=0.20',
                    'seaborn>=0.7',
                    'scikit-learn>=0.23.2',
                    'pyfaidx>=0.5.4',
                    'pysam>=0.10.0',
                    'HTSeq>=0.6',
                    'bx-python>=0.5',
                    'forgi==1.1',
                    'logging_exceptions']

major, minor, micro = sys.version_info[:3]
if major == '2':
    inst_requires.append('future')

setup(
    name = 'smallrnaseq',
    version = '0.6.0',
    description = 'Package for short RNA-seq analysis',
    long_description = 'smallrnaseq is a Python package for processing of small RNA seq data.',
    url='https://github.com/dmnfarrell/smallrnaseq',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['smallrnaseq'],
    package_data={'smallrnaseq': ['data/*.*','data/de_example/*','*.R']},
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
