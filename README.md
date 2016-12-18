<img align="right" src=https://raw.githubusercontent.com/dmnfarrell/smallrnaseq/master/img/logo.png width=150px>

# smallrnaseq -  Scripts for analysis of small RNA sequencing data e.g. microRNAs

#### Background

This package is a collection of Python modules originally designed to be used in conjunction with miRDeep2 and sRNAbench. It also includes general routines for analysis of small RNA sequencing data in general. The modules may prove useful to those wishing a convenient interface for small RNA or miRNA data analysis.
 
Functionality includes:

* ability to use several different short read aligners, (currently bowtie or subread)
* counting of mapped reads to annotated sequences or genomic features
* running multiple files in a batch
* comparison across samples
* use of several normalization methods
* command line functions for common tasks without needeing Python code

#### Installation

Simply clone the git repository and add the folder into your Python path. Available on pip soon.

#### Usage

See the [wiki page](https://github.com/dmnfarrell/smallrnaseq/wiki) for documentation.

#### Optional software

1. [miRDeep2](https://www.mdc-berlin.de/8551903/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep "mirdeep2")
2. [sRNAbench](http://bioinfo5.ugr.es/sRNAbench/sRNAbench.php "sRNAbench")
3. [bowtie, short read aligner](http://bowtie-bio.sourceforge.net/index.shtml)

##### Python Dependencies
* Numpy
* Pandas
* Matplotlib
* BioPython
* HTSeq
* Seaborn
