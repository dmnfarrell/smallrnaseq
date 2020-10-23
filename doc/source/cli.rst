Using smallrnaseq
=================

This page refers to using the Command Line Interface of smallrnaseq. For programmers using the API see code examples.

Installing the package provides the command `smallrnaseq` in your path. This allows users is a command line interface to the library without the need for any Python coding at all. It provides a set of pre-defined functions with parameters specified in a text configuration file. The primary input is one or more fastq files containing short read data from a small rna-seq experiment. Note: It is assumed they have already been adapter trimmed, if required.

Usage
-----

Usage largely involves setting up the config file and having your input files prepared. Running the command `smallrnaseq -c default.conf` will create a new config file for you to work from if none exists. Just edit this with a text editor and then to execute::

  smallrnaseq -c default.conf -r

Several other functions are available from the command line without the config file, i.e. to collapse or trim reads. Type `smallrnaseq -h` to get a list of options.

Configuration file settings
---------------------------

The advantage of configuration files is in avoiding long commands that have to be remembered or are prone to mistakes. Also the config files can be kept to recall what setting we used or to copy them for another set of files. The current options available in the file are shown below. The meaning of each option is explained explained below.  If you are unsure or don't require an option value, leave it at the default. Options can be excluded from the file completely and defaults will be used but it's recommended to just leave unused options blank. You can also comment out lines with a '#' at the start if you want it to be ignored. The [base] heading should always be present an indicates a section of the file. The [aligner] section is for setting alignment parameters on a per library basis if you need to. The [novel] section is for a few novel mirna prediction settings. The [de] section is used for  differential expression which can be ignored if you don't need that::

  [base]
  filenames =
  path =
  overwrite = 0
  index_path = indexes
  libraries =
  ref_fasta =
  features =
  output = results
  add_labels = 1
  aligner = bowtie
  mirna = 0
  species = hsa
  pad5 = 3
  pad3 = 5
  verbose = 1
  cpus = 1

  [aligner]
  default_params = -v 1 --best
  mirna_params = -v 1 -a --best --strata --norc

  [novel]
  score_cutoff = 0.8
  read_cutoff = 50

  [de]
  sample_labels =
  sep = ,
  sample_col =
  factors_col =
  conditions =
  logfc_cutoff =


Settings explained:

+---------------+---------------------------------+-------------------+
| name          | example value                   | meaning           |
+===============+=================================+===================+
| filenames     | test.fastq                      | input fastq       |
|               |                                 | file(s) with      |
|               |                                 | reads, comma      |
|               |                                 | separated         |
+---------------+---------------------------------+-------------------+
| path          | testfiles                       | folder containing |
|               |                                 | fastq files       |
|               |                                 | instead of        |
|               |                                 | listing files     |
+---------------+---------------------------------+-------------------+
| index_path    | indexes                         | location of       |
|               |                                 | bowtie or subread |
|               |                                 | indexes           |
+---------------+---------------------------------+-------------------+
| aligner       | bowtie                          | which aligner to  |
|               |                                 | use, bowtie or    |
|               |                                 | subread           |
+---------------+---------------------------------+-------------------+
| ref_fasta     | hg19                            | reference genome  |
|               |                                 | fasta file,       |
|               |                                 | optional          |
+---------------+---------------------------------+-------------------+
| libraries     | RFAM_human,mirbase-hsa          | names of          |
|               |                                 | annotated library |
|               |                                 | indexes to map to |
+---------------+---------------------------------+-------------------+
| features      | Homo_sapiens.GRCh37.75.gtf      | genome annotation |
|               |                                 | file. ONLY needed |
|               |                                 | for counting      |
|               |                                 | genomic features  |
+---------------+---------------------------------+-------------------+
| output        | smrna_results                   | output folder for |
|               |                                 | temp files        |
+---------------+---------------------------------+-------------------+
| add_labels    | 1                               | whether to add    |
|               |                                 | labels to replace |
|               |                                 | the file names in |
|               |                                 | the results (0 or |
|               |                                 | 1)                |
+---------------+---------------------------------+-------------------+
| mirna         | 0                               | run mirna         |
|               |                                 | counting workflow |
|               |                                 | (0 or 1)          |
+---------------+---------------------------------+-------------------+
| species       | bta                             | mirbase species   |
|               |                                 | to use            |
+---------------+---------------------------------+-------------------+
| pad5          | 3                               | 3’ flanking bases |
|               |                                 | to add when       |
|               |                                 | generating mature |
|               |                                 | mirbase sequences |
+---------------+---------------------------------+-------------------+
| pad3          | 5                               | 5’ flanking bases |
|               |                                 | to add            |
+---------------+---------------------------------+-------------------+
| verbose       | 1                               | print extra       |
|               |                                 | information       |
+---------------+---------------------------------+-------------------+
| cpus          | 1                               | number of threads |
|               |                                 | to use            |
+---------------+---------------------------------+-------------------+
| sample_labels | samplefile.txt                  | csv file with     |
|               |                                 | sample labels     |
+---------------+---------------------------------+-------------------+
| default_param | -v 1 –best                      | default alignment |
| s             |                                 | parameters        |
+---------------+---------------------------------+-------------------+
| mirna_params  | -v 1 -a –best –strata –norc     | default miRNA     |
|               |                                 | alignment         |
|               |                                 | parameters        |
+---------------+---------------------------------+-------------------+

Examples
--------

Mapping to one or more libraries of known RNAs
++++++++++++++++++++++++++++++++++++++++++++++

This will simply consist of setting the `libraries` option to match the names of files you have created aligner indexes for. You can download RNA annotations in fasta format from http://rnacentral.org/. The index files should be placed in the `index_path` folder. To create an index using bowtie you supply a fasta file and the indexes are placed in a folder called 'indexes'. You can move them to another folder if needed::

  smallrnaseq -b myrnas.fa

By default file names are replaced with short ids of the form s01, s02, etc. This also writes out a file called `sample_labels.csv` in the output folder. Set `add_labels = 0` if you don't want this behaviour.

Mapping to known miRNAs and finding novel miRs
++++++++++++++++++++++++++++++++++++++++++++++

Say we have a set of fastq files in the folder 'testfiles' that we want to count miRNAs in. We would set the options `mirna = 1` and `path = testfiles`. Note if just mapping to mirbase mature/precursor sequence you don't have to create an index file since it is generated automatically. If you are using any other libraries you should create the aligner index first. For novel miRNA discovery we need the reference genome and index for that.

Mapping to a reference genome with small RNA features
+++++++++++++++++++++++++++++++++++++++++++++++++++++

There are some advantages to using a full reference genome in that it allows us to see of the reads align outside the target transcripts and we might therefore exclude them. Also we can normalise the read counts by the sum of all mapped reads. This depends on what features you have in the gtf file. To count our reads using features in a reference genome, provide the `ref_genome` name which should correspond to the index name of the reference sequence. We also need to set `features`, a set of features stored in a gtf, gff or bed file. Remember that the coordinates in this file must correspond to the reference sequence. See Counting features for more information.

Outputs
-------

The main outputs are csv files with the counts for each sample in a column, also produced is a 'long form' csv file with a separate row for every sample. These csv files can be opened in a spreadsheet. Temporary files such as collapsed read files are placed in a 'temp' folder inside the output folder.

Example
~~~~~~~

When you run the mirna workflow on a set of samples, a file called mirbase_mature_counts.csv is created in the folder. The column names are mostly self explanatory. Each sample has a column for the raw counts and normalized counts. `mean_norm` is the normalized mean and `total_reads` is the sum of all raw read counts::

  name,db,sample1,sample2,sample3,sample1 norm,sample2 norm,sample3 norm,total_reads,mean_norm
  bta-miR-486,mirbase-bta,1444,1070,5579,176722.56,47917.6,127569.57,239793,284398.062
  bta-miR-122,mirbase-bta,693,10676,4,84812.14,478101.21,91.46,11409,149778.2825
  bta-miR-423-5p,mirbase-bta,337,100,2008,41243.42,4478.28,45914.98,4270,62353.636
  bta-miR-22-3p,mirbase-bta,315,1053,5655,38550.97,47156.29,129307.39,7462,45154.618
  bta-miR-92a,mirbase-bta,332,891,2484,40631.5,39901.48,56799.21,5989,39094.122
  bta-miR-192,mirbase-bta,659,656,1482,80651.08,29377.52,33887.45,2945,33098.118
  bta-miR-21-5p,mirbase-bta,453,227,1311,55439.97,10165.7,29977.36,2054,27358.996
  ..

Differential expression
-----------------------

This workflow is done using a separate command and via some extra options not shown above for clarity. To execute this type of analysis you must have done gene counting (e.g. miRNAs) and have the results in a csv file. The analysis is then run using::

  smallrnaseq -c default.conf -d

In the default file the additional options needed are in a separate [de] section. Most are blank by default. If you want to do differential expression analysis of the genes from your results the other main thing you need to provide is a **sample label file** that matches the filenames you have analysed to the conditions. You then choose which conditions you want to compare. The options are explained below.

+---------------+---------------------------------+-------------------+
| name          | example value                   | meaning           |
+===============+=================================+===================+
| count_file    | mirna_counts.csv                | csv file with     |
|               |                                 | gene counts       |
+---------------+---------------------------------+-------------------+
| sample_labels | labels.txt                      | csv/text file     |
|               |                                 | with sample       |
|               |                                 | labels            |
+---------------+---------------------------------+-------------------+
| sep           | ,                               | separator for     |
|               |                                 | sample labels     |
|               |                                 | text file         |
+---------------+---------------------------------+-------------------+
| sample_col    | file_name                       | column for the    |
|               |                                 | sample file names |
+---------------+---------------------------------+-------------------+
| factors_col   | status                          | column for factor |
|               |                                 | labels            |
+---------------+---------------------------------+-------------------+
| conditions    | control,infected                | conditions to     |
|               |                                 | compare from      |
|               |                                 | factor column     |
|               |                                 | (each replicate   |
|               |                                 | will have the     |
|               |                                 | same condition)   |
+---------------+---------------------------------+-------------------+
| logfc_cutoff  | 1.5                             | cutoff for log    |
|               |                                 | fold changes      |
|               |                                 | (lower are        |
|               |                                 | ignored)          |
+---------------+---------------------------------+-------------------+

This analysis is run using the R edgeR package, so R must be installed. See Installing R.

Example
+++++++

The sample label file below shows how we use the above options. Our filenames are in the Run_s column. We want to compare the conditions 3 and 6 months samples in the age_s (the 'factor') column. mirna_mature_counts.csv is the file where counts from our mirna analysis are stored. This should have a column for each sample. Note that you could use counts output from any program as long as they are csv and have the right column names, see the output file formats section.
It is assumed you are using replicates. Note that column names are case sensitive::

  Run_s	age_s	isolate_s
  SRR3457948	3 months	animal1
  SRR3457949	6 months	animal1
  SRR3457950	6 months	animal4
  SRR3457951	15 months	animal4
  SRR3457952	3 months	animal5
  SRR3457953	6 months	animal5
  SRR3457954	15 months	animal5
  SRR3457955	3 months	animal6
  SRR3457956	6 months	animal6
  ...

So the config file will look like this::

  [de]
  sample_labels = SraRunTable.txt
  sep = tab #tab delimiter
  count_file = mirna_mature_counts.csv
  sample_col = Run_s
  factors_col = age_s
  conditions = 3 months,6 months
  logfc_cutoff = 1.5


Adapter trimming
----------------

Since there are numerous programs to perform this task it is left to the user to perform trimming of the reads prior to input.

Aligners
--------

Traditional sequence alignment tools like BLAST are not well suited for next generation sequencing where one needs to align millions of short sequences very quickly. This has given rise to the development of a new class of short read aligners of which there are now dozens available. For small RNA-seq data alignment present certain specific challenges but standard aligners used for normal RNA-seq are usually adequate with the right parameters. The aligner and parameters will have an effect on the results, so trying more than one might be a good idea.

Currently smallrnaseq integrates bowtie (version 1) and subread. Though others can be added on request. These are free and easy to install on linux and OSX systems. On Ubuntu the following command installs both::

  apt install bowtie subread

Links
+++++

 * http://subread.sourceforge.net/
 * http://bowtie-bio.sourceforge.net/index.shtml

Alignment settings
------------------

The parameters used for the alignment/mapping procedure can be important
in the final counts produced, irrespective of the aligner used.

This can be a complex topic in itself and general users will be confused
by the many options. The command line tool for smallrnaseq takes the
simple approach of providing a default alignment parameter for general
mapping to libraries, another for mapping miRNAs and one for reference
genomes. All can be changed in the config file if needed. You can also
set custom parameters per library in the aligner section.

Bowtie
++++++

For general mapping ``-v 1 --best`` is used. ``-v 1`` reports read
mappings with up to one mismatch, options ``--best`` orders the mappings
from best to worse alignments.

In miRDeep2 when mapping to the mature miRNAs (miRBase sequences) for
mature quantification the following parameters are used::

 -v 1 -a --best --strata --norc

Here ``-a`` means report all valid alignments, options
``--best --strata`` orders the mappings from best to worse alignments
according to the strata definition of bowtie. ``--norc`` means do not
map reads to the reverse complement of the sequences.

For reference genome mapping miRDeep2 uses these parameters::

 -n 0 -e 80 -l 18 -a -m 5 --best --strata

Subread
+++++++

``-m 2 -M 1`` is the default for general alignment to libraries. If you
use subread you can check the parameters by typing
``subread-align --help`` at the command line, or refer to the website.

Configuration file
++++++++++++++++++

In the aligner section set your parameters. In the example below
bos_taurus is the name of the reference genome. We have also used custom
settings for mirna and another library of tRNAs.

::

  [aligner]
  mirna_params = -n 1 -l 20
  bos_taurus = -v 1 -k 50
  bosTau8-tRNAs = -v 0 --best

Code example
++++++++++++

If using the package in your python code, aligner parameters are set via
the aligners module. This is done before calling mapping routines such
as ``map_rnas``.

for example::

  from smallrnaseq import aligners
  aligners.BOWTIE_PARAMS = '-v 0 --best'
  aligners.SUBREAD_PARAMS = '-m 0 -M 1'

References
----------

-  Shi, J., Dong, M., Li, L., Liu, L., Luz-Madrigal, A., Tsonis, P. A.,
  … Liang, C. (2015). mirPRo–a novel standalone program for
  differential expression and variation analysis of miRNAs. Scientific
  Reports, 5, 14617. http://doi.org/10.1038/srep14617

-  Friedländer, M. R., Mackowiak, S. D., Li, N., Chen, W., & Rajewsky,
  N. (2012). miRDeep2 accurately identifies known and hundreds of novel
  microRNA genes in seven animal clades. Nucleic Acids Research, 40(1),
  37–52. http://doi.org/10.1093/nar/gkr688
