#!/usr/bin/env python

"""
    mirnaseq unit tests
    Created September 2016
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os
import pandas as pd
import unittest
from . import config, base, utils, analysis, mirdeep2, aligners, novel

class BasicTests(unittest.TestCase):
    """Basic tests for mirnaseq"""

    def setUp(self):
        self.df = None
        self.testdir = 'testing'
        if not os.path.exists(self.testdir):
            os.mkdir(self.testdir)
        aligners.BOWTIE_INDEXES = 'testing/indexes'
        return

    def test_collapse_reads(self):
        """collapse reads test"""

        f = os.path.join(base.datadir, 'bovine_serum_sample.fastq')
        base.collapse_reads(f)
        return

    def test_build_index(self):
        filename = os.path.join(base.datadir, 'bosTau8-tRNAs.fa')
        aligners.build_bowtie_index(filename, 'testing/indexes')
        #aligners.build_subread_index(filename, 'testing/indexes')
        return

    def test_htseq(self):
        """htseq basic test for sam file reading"""

        import HTSeq
        samfile = os.path.join(base.datadir, 'test.sam')
        sam = HTSeq.SAM_Reader(samfile)
        f=[]
        for a in sam:
            if a.aligned == True:
                seq = a.read.seq.decode()
                f.append((seq,a.read.name,a.iv.chrom))
        df = pd.DataFrame(f, columns=['seq','read','name'])
        #print (df)
        return

    def test_pandas(self):
        """dataframe-fasta tests"""

        f = os.path.join(base.datadir, 'test_collapsed.fa')
        df = utils.fasta_to_dataframe(f)
        #utils.dataframe_to_fasta(df)
        return

    def test_read_aligned(self):
        """read in alignment/counts test"""

        samfile = os.path.join(base.datadir, 'test.sam')
        cfile = os.path.join(base.datadir, 'test_collapsed.fa')
        reads = utils.get_aligned_reads(samfile, cfile)
        return

    def test_count_aligned(self):
        """count aligned test"""

        samfile = os.path.join(base.datadir, 'test.sam')
        cfile = os.path.join(base.datadir, 'test_collapsed.fa')
        counts = base.count_aligned(samfile, cfile)
        #print (counts)
        return

    def test_count_features(self):
        """feature counting"""

        #countfile = os.path.join(outpath, '%s.csv' %label)
        #readcounts = pd.read_csv(countfile, index_col=0)
        #hits = count_features(samfile, features=exons, truecounts=readcounts)
        return

    def test_map_rnas(self):
        """mapping to libraries"""

        aligners.BOWTIE_PARAMS = '-v 0 --best'
        fastafile = os.path.join(base.datadir, 'bosTau8-tRNAs.fa')
        aligners.build_bowtie_index(fastafile, path='testing/indexes')
        path = os.path.join(self.testdir, 'ncrna_map')
        f = os.path.join(base.datadir, 'bovine_serum_sample.fastq')
        res = base.map_rnas([f], ['bosTau8-tRNAs'], path, overwrite=True, aligner='bowtie')
        return

    def test_map_mirnas(self):
        """mirna counting test"""

        f = os.path.join(base.datadir, 'bovine_plasma_sample.fastq')
        path = os.path.join(self.testdir, 'ncrna_map')
        #res = base.map_mirbase(files=[f], overwrite=True, outpath=path,
        #                       aligner='bowtie', species='bta')
        return

    def test_map_features(self):
        """Genomic feature mapping/counting"""

        aligners.BOWTIE_PARAMS = '-v 0 --best'
        #fastafile = os.path.join(base.datadir, 'bostau.fa')
        #base.build_bowtie_index(fastafile)
        path = os.path.join(self.testdir, 'ncrna_map')
        f = os.path.join(base.datadir, 'bovine_serum_sample.fastq')
        return

    def test_mirdeep(self):
        """mirdeep2 script test"""

        conffile = 'testing/test_mdp.conf'
        cp = config.parse_config(conffile)
        options = config.get_options(cp)
        #print (options)
        #mirdeep2.run_multiple(**options)
        return

def run():
    unittest.main()

if __name__ == '__main__':
    unittest.main()
