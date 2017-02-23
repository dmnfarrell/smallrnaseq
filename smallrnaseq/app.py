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
import pandas as pd
from smallrnaseq import config, base, analysis, utils, aligners, plotting, de

def run(opts):
    """Run predefined batch mapping routines based on options from
     a config file"""

    opts = config.check_options(opts)
    fastafile = opts['filename']
    path = opts['path']
    out = opts['output']
    temp_path = os.path.join(out,'temp') #path for temp files
    indexes = opts['indexes'].split(',')
    species=opts['species']

    if path != '':
        files = glob.glob(os.path.join(path,'*.fastq'))
    elif fastafile != '':
        files = [fastafile]
    else:
        print ('you should provide at least one file or folder')
        return

    aligners.BOWTIE_INDEXES = aligners.SUBREAD_INDEXES = opts['index_path']
    if opts['aligner_params'] != '':
        aligners.set_params(opts['aligner'], opts['aligner_params'])

    if not os.path.exists(opts['index_path']):
        print ('no such folder for indexes')
        return
    if opts['mirbase'] == 1:
        #count with mirbase
        print ('mapping to mirbase')
        res, counts = base.map_mirbase(files, outpath=temp_path, pad5=3,
                         aligner=opts['aligner'], species=species,
                         add_labels=opts['add_labels'])
        res.to_csv( os.path.join(out, 'mirbase_mature_found.csv'),index=False )
        counts.to_csv( os.path.join(out, 'mirbase_mature_counts.csv'), index=False )
        plot_results(res, out)
        #isomir counting
        iso, isocounts = base.map_isomirs(files, temp_path, species)
        iso.to_csv( os.path.join(out, 'isomirs_found.csv'),index=False )

    elif opts['ref_genome'] != '':
        print ('mapping to reference genome')
        res = base.map_genome_features(files, opts['ref_genome'], opts['features'],
                                       outpath=temp_path, aligner=opts['aligner'])
        counts = base.pivot_count_data(res, idxcols=['name','gene_name','gene_biotype'])
        res.to_csv( os.path.join(out, 'features_found.csv'), index=False )
        counts.to_csv( os.path.join(out, 'feature_counts.csv'), index=False)
        print ('results saved to feature_counts.csv')
    else:
        #map to provided libraries
        print ('mapping to these libraries: %s' %indexes)
        res, counts = base.map_rnas(files, indexes, temp_path, aligner=opts['aligner'],
                                    add_labels=opts['add_labels'])
        if res is None:
            print ('empty data returned. did alignments run?')
            return
        print ('results saved to rna_counts.csv')
        res.to_csv( os.path.join(out, 'rna_found.csv'),index=False)
        counts.to_csv( os.path.join(out, 'rna_counts.csv'), index=False )
        plot_results(res, out)
    print ('intermediate files saved to %s' %out)
    return

def plot_results(res, path):
    """Some results plots"""

    if res is None or len(res) == 0:
        return
    counts = base.pivot_count_data(res, idxcols=['name','ref'])
    x = base.get_fractions_mapped(res)
    print (x)
    fig = plotting.plot_fractions(x)
    fig.savefig(os.path.join(path,'fractions_mapped.png'))
    fig = plotting.plot_sample_counts(counts)
    fig.savefig(os.path.join(path,'total_per_sample.png'))
    fig = plotting.plot_read_count_dists(counts)
    fig.savefig(os.path.join(path,'distr_per_sample.png'))
    scols,ncols = base.get_column_names(counts)
    if len(scols)>1:
        fig=plotting.expression_clustermap(counts)
        fig.savefig(os.path.join(path,'expr_map.png'))
    return

def build_indexes(filename):
    path = 'indexes'
    aligners.build_bowtie_index(filename, path)
    aligners.build_subread_index(filename, path)
    return

def diff_expression(opts):
    """Diff expression workflow"""

    path = opts['output']
    labelsfile = opts['sample_labels']
    countsfile = opts['count_file']
    for f in [labelsfile,countsfile]:
        if not os.path.exists(f) or f == '':
            print ('no such file %s!' %f)
            print_help()
            return
    labels = pd.read_csv(labelsfile, sep=opts['sep'])
    counts = pd.read_csv(countsfile)

    #define sample/factor cols and conditions for de
    samplecol = opts['sample_col']
    factorcol = opts['factors_col']
    conds = opts['conditions'].split(',')
    print (conds)
    #get the samples needed for the required conditions we want to compare
    data = de.get_factor_samples(counts,
                                 labels, [(factorcol,conds[0]),(factorcol,conds[1])],
                                 samplecol=samplecol, index='name')
    #print (data[:4])
    res = de.run_edgeR(data=data, cutoff=1.5)
    res.to_csv(os.path.join(path,'de_genes.csv'))
    print (res)
    names = res.name

    #plot these genes with seaborn
    xorder=conds
    m = de.melt_samples(counts, labels, names, samplecol=samplecol)
    g = base.sns.factorplot('age_s','read count', data=m, col='name', kind="point",
                            s=10, lw=1, col_wrap=4, size=4, aspect=1.2,
                            legend_out=True,sharey=False, order=xorder)
    g.savefig(os.path.join(path,'de_genes.png'))
    return

def print_help():
    """generic help message"""

    print ('to run a workflow use smallrnaseq -c <config> -r')
    print ('see https://github.com/dmnfarrell/smallrnaseq/wiki')
    print ('for further information')
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
                        help="input file", metavar="FILE")
    parser.add_option("-l", "--collapse", dest="collapse", action="store_true",
                        default=False, help="collapse reads in input file")
    parser.add_option("-t", "--trim", dest="trim",  type='string',
                        help="trim given adapter in input file")
    parser.add_option("-d", "--de", dest="de",  action="store_true",
                        default=False, help="run DE analysis")

    opts, remainder = parser.parse_args()
    if opts.infile != None:
        if opts.trim != None:
            utils.trim_adapters(opts.infile, adapters=opts.trim, outfile='cut.fastq')
        if opts.collapse == True:
            base.collapse_reads(opts.infile)
    elif opts.build != None:
        build_indexes(opts.build)
    elif opts.config != None and os.path.exists(opts.config):
        cp = config.parse_config(opts.config)
        options = config.get_options(cp)
        print ('using the following options:')
        print ('----------------------------')
        config.print_options(options)
        if opts.run == True:
            run(options)
        elif opts.de == True:
            diff_expression(options)
        else:
            print_help()
    else:
        if opts.config != None:
            conffile = opts.config
        else:
            print ('No config file provided.')
            conffile = 'default.conf'
        config.write_default_config(conffile, defaults=config.baseoptions)
        print_help()
        return

if __name__ == '__main__':
    main()
