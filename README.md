<img align="right" src=https://raw.githubusercontent.com/dmnfarrell/smallrnaseq/master/img/logo.png width=150px>

**smallrnaseq -  Scripts for analysis of small RNA sequencing data e.g. microRNAs**

#### Background

This repository is a collection of Python modules originally designed to be used in conjunction with miRDeep2 and sRNAbench. It also includes routines for analysis of small RNA sequencing data in general. The modules may prove useful to thos wishing a convenient interface for miRNA data analysis.

#### Installation

Simply clone the git repository and add the folder into your Python path. available on pip soon.

#### Usage

**mirdeep2.py**: This script allows you to run mirdeep2 and/or analyse the output.

If you are using Python the modules can be used in another script. This assumes the miRDeep folder is in your path. Examples:
```
from mirnaseq import mirdeep2
#run mirdeep2 on multiple fastq files, stored in inpath
mirdeep2.runMultiple(inpath, fasta=False)
#get mirdeep results as a Pandas DataFrame
df = mirdeep2.getResults(path)
#filter the mirnas according to various parameters
df = mirdeep2.filterExprResults(score=0, freq=30, meanreads=500)
```

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
