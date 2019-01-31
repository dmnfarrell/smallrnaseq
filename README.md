# smallrnaseq - tool for analysis of small RNA sequencing data

<img align="right" src=https://raw.githubusercontent.com/dmnfarrell/smallrnaseq/master/img/logo.png width=150px>

smallrnaseq is a Python package that integrates some of the standard approaches for convenient pre-processing, quantification and analysis of small RNAs. It includes a [command line interface](https://github.com/dmnfarrell/smallrnaseq/wiki/Command-line-interface) used in conjunction with a configuration file. Hence knowledge of Python is not a necessity for using this software. The aim is to provide convenience but with enough configuration to allow some flexibility. Supported on linux and OSX.

Functionality includes:

* command line interface
* ability to use several different short read aligners, (currently bowtie or subread)
* counting of mapped reads to annotated sequences or genomic features
* counting of miRNAs and isomiRs
* novel miRNA prediction
* running multiple files in a batch
* use of several normalization methods

### Usage

Up to date documentation is at http://smallrnaseq.readthedocs.io

There is a video tutorial at https://youtu.be/m24cuLyTqg0

### Installation

```pip install smallrnaseq```

or 

```snap install smallrnaseq```

See the [installation page](http://smallrnaseq.readthedocs.io/en/latest/description.html#installation) for more details 

### Python Dependencies

* Numpy
* Pandas
* Matplotlib
* HTSeq
* Seaborn
* bx-python
* scikit-learn
