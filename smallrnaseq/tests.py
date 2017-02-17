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
from . import config, base, analysis, mirdeep2, srnabench

class BasicTests(unittest.TestCase):
    """Basic tests for mirnaseq"""

    def setUp(self):
        self.df = None
        self.testdir = 'testing'
        if not os.path.exists(self.testdir):
            os.mkdir(self.testdir)
        base.BOWTIE_INDEXES = 'indexes'
        return

    def test_count_features(self):
        """feature counting"""

        return

    def test_map_rnas(self):
        """Generic mapping to rna annotations"""

        base.BOWTIE_PARAMS = '-v 0 --best'
        fastafile = os.path.join(base.datadir, 'bosTau8-tRNAs.fa')
        base.build_bowtie_index(fastafile)
        path = os.path.join(self.testdir, 'ncrna_map')
        f = os.path.join(base.datadir, 'bovine_serum_sample.fastq')
        res = base.map_rnas([f], ['bosTau8-tRNAs'], path, overwrite=True, aligner='bowtie')
        return

    def test_map_mirnas(self):
        """mirna counting"""

        f = os.path.join(base.datadir, 'bovine_plasma_sample.fastq')
        path = os.path.join(self.testdir, 'ncrna_map')
        res = base.map_mirbase(files=[f], overwrite=True, outpath=path,
                               aligner='bowtie', species='bta')
        res = base.map_mirbase(files=[f], overwrite=True, outpath=path,
                               aligner='subread', species='bta')
        print (len(res))
        return

    def test_map_features(self):
        """Genomic feature mapping/counting"""

        base.BOWTIE_PARAMS = '-v 0 --best'
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
