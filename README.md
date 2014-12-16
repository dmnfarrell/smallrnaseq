mirnaseq
========

#### Scripts for analysis of microRNA sequencing data, for use with mirDeep2 and sRNAbench

#### Background

This repository is a collection of Python modules partly designed to be used in conjunction with miRDeep2. This is a program for analysis of small RNA sequencing data. These scripts were used for all analysis in our paper. [reference]. The code is deposited here primarily to allow our work to be replicated. However the modules may prove useful to other miRDeep users.

#### Installation

Simply clone the git repository and add the folder into your Python path.

#### Usage

**mirdeep2.py**: This script allows you to run mirdeep2 and/or analyse the output.

If you are using Python the modules can be used in another script. This assumes he miRDeep folder is in your path. Examples:
```
from mirnaseq import mirdeep2
#run mirdeep2 on multiple fastq files, stored in inpath
mirdeep2.runMultiple(inpath, fasta=False)
#get mirdeep results as a Pandas DataFrame
df = mirdeep2.getResults(path)
#filter the mirnas according to various parameters
df = mirdeep2.filterExprResults(score=0, freq=30, meanreads=500)
```

#### Required software

1. [miRDeep2](https://www.mdc-berlin.de/8551903/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep "mirdeep2")
2. [sRNAbench](http://bioinfo5.ugr.es/sRNAbench/sRNAbench.php "sRNAbench")

##### Dependencies
* Numpy
* Pandas
* Matplotlib
* BioPython
