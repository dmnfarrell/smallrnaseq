<img align="right" src=https://raw.githubusercontent.com/dmnfarrell/smallrnaseq/master/img/logo.png width=150px>

# smallrnaseq -  Python package for analysis of small RNA sequencing data such as microRNAs

#### Background

This package is a collection of Python modules designed for general analysis of small RNA-seq data. The package integrates some of the standard approaches for convenient pre-processing, quantification and analysis of small RNAs. It includes a command line interface used in conjunction with a configuraiton file. Hence knowledge of Python is not a necessity for using this software. The aim is to provide convenience but with enough configuration to allow some flexibility

Functionality includes:

* command line interface
* ability to use several different short read aligners, (currently bowtie or subread)
* counting of mapped reads to annotated sequences or genomic features
* counting of mirnas and isomirs with mirbase 
* running multiple files in a batch
* comparison across samples
* use of several normalization methods
* differential expression of counted genes

#### Installation

pip install smallrnaseq

#### Usage

See the [wiki page](https://github.com/dmnfarrell/smallrnaseq/wiki) for documentation.

#### Optional software

1. [miRDeep2](https://www.mdc-berlin.de/8551903/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep "mirdeep2")
2. [sRNAbench](http://bioinfo5.ugr.es/sRNAbench/sRNAbench.php "sRNAbench")

##### Python Dependencies

* Numpy
* Pandas
* Matplotlib
* BioPython
* HTSeq
* Seaborn
