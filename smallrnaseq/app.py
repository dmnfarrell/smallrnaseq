#!/usr/bin/env python

"""
    smallrnaseq entry point for command line functions
    Created November 2016
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

from __future__ import absolute_import, print_function
import sys, os, string, types
import glob
from smallrnaseq import config, base, analysis, utils, aligners, plotting

def run(opts):
    """Run batch mapping routines based on conf file"""

    #print (opts)
    fastafile = opts['filename']
    path = opts['path']
    out = opts['output']
    indexes = opts['indexes'].split(',')
    if path != '':
        files = glob.glob(os.path.join(path,'*.fastq'))
    elif fastafile != '':
        files = [fastafile]
    else:
        print ('ou should provide at least one file or folder')
        return

    aligners.BOWTIE_INDEXES = opts['index_path']
    aligners.BOWTIE_PARAMS = opts['bowtie_params']

    if opts['mirbase'] == 1:
        #count with mirbase
        print ('mapping to mirbase')
        base.map_mirbase(files, outpath=out, pad5=3,
                         aligner='bowtie', species=opts['species'],
                         add_labels=opts['add_labels'])
    elif opts['ref_genome'] != '':
        print ('mapping to reference genome')
        res = base.map_genome_features(files, opts['ref_genome'], opts['features'],
                                 outpath=out, aligner='bowtie')
        counts = base.pivot_count_data(res, idxcols=['name','gene_name','gene_biotype'])
        counts.to_csv( 'feature_counts.csv', index=False)
        print ('results saved to feature_counts.csv')
    else:
        #map to provided libraries
        print ('mapping to these libraries: %s' %indexes)
        res = base.map_rnas(files, indexes, out, aligner='bowtie',
                            add_labels=opts['add_labels'])
        if res is None:
            print ('empty data returned. did alignments run?')
            return
        print ('results saved to rna_counts.csv')
        x = base.get_fractions_mapped(res)
        print (x)
        plotting.plot_fractions(x)
    print ('intermediate files saved to %s' %out)
    return

def build_indexes(filename):
    path = 'indexes'
    aligners.build_bowtie_index(filename, path)
    aligners.build_subread_index(filename, path)
    return

def main():
    """Run the application from outside the module - used for
       deploying as frozen app"""

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                        help="config file to use", metavar="FILE")
    parser.add_option("-r", "--run", dest="run",  action="store_true",
                        default=False, help="run smallrna mapping")
    parser.add_option("-b", "--build", dest="build",
                        help="build an index for given file", metavar="FILE")
    parser.add_option("-f", "--infile", dest="infile",
                        help="collapse reads in input file", metavar="FILE")
    parser.add_option("-l", "--collapse", dest="collapse", action="store_true",
                        default=False, help="collapse reads in input file")
    parser.add_option("-t", "--trim", dest="trim",  type='string',
                        help="trim given adapter")

    opts, remainder = parser.parse_args()
    if opts.infile != None:
        if opts.trim != None:
            utils.trim_adapters(opts.infile, adapters=opts.trim, outfile='cut.fastq')
        if opts.collapse == True:
            base.collapse_reads(opts.infile)
    elif opts.build != None:
        build_indexes(opts.build)
    elif opts.run == True:
        if opts.config == None:
            print ('No config file provided.')
            config.write_default_config('default.conf', defaults=config.baseoptions)
            return
        cp = config.parse_config(opts.config)
        options = config.get_options(cp)
        print ('using the following options:')
        print ('----------------------------')
        config.print_options(options)
        run(options)
    else:
        print ('you need a config file')
        #run tests
        #from smallrnaseq.tests import BasicTests
        #import unittest
        #unittest.main()

if __name__ == '__main__':
    main()
