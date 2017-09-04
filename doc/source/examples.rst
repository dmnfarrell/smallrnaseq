Code Examples
=============

Counting microRNAs
++++++++++++++++++

Requires you provide at least one fastq file and have a short read aligner installed. You should also
specify the species being mapped to::

    import smallrnaseq as smrna
    res = smrna.map_mirbase(files=['test_1.fastq','test_2.fastq'], overwrite=True, aligner='bowtie',
                        species='hsa', pad5=3, pad3=5)

Counting isomiRs
++++++++++++++++

This method is used to count isomiRs using results from previously mapped reads. So a sam file is
required::

    smrna.count_isomirs(samfile, truecounts, species='bta')

Mapping to the genome
+++++++++++++++++++++

This requires a reference genome and a gtf file with miRNA features::

    featcounts = srseq.map_genome_features(['test_1.fastq'], 'bos_taurus', gtffile,
                                        outpath='ncrna_map', aligner='subread', merge=True)

Differential Expression
+++++++++++++++++++++++

Assuming we have all the raw files, they need to be adapter trimmed.
Optionally you can remove other ncrnas before counting your target rnas class,
though that may not be advisable.
The following code maps all the files to bovine mature miRNAs and counts the mapped genes,
then saves the results to a csv file which has the counts in one column per sample.
You can skip this if you already have the counts file::

    import pandas as pd
    import smallrnaseq as smrna
    from smallrnaseq import base, utils, de

    path = 'pathtodata'
    base.BOWTIE_INDEXES = 'bowtie_index'
    refs = ['mirbase-bta'] #name of bowtie index

    files = glob.glob(path+'/*.fastq')
    outpath = 'ncrna_map'
    #map to selected annotation files
    counts = smrna.map_rnas(files, refs, outpath, overwrite=True)
    R = smrna.pivot_count_data(counts, idxcols=['name','db'])
    R.to_csv('mirna_counts.csv',index=False)