#!/usr/bin/env python

"""
    Module for core smallrnaseq functions
    Created July 2014
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
import sys, os, string, time
import types, re, subprocess, glob, shutil
import pylab as plt
import numpy as np
import pandas as pd
import configparser
try:
    import HTSeq
except:
    'HTSeq not present'

from . import utils
import matplotlib as mpl
import seaborn as sns
sns.set_style("ticks", {'axes.facecolor': '#F7F7F7',
                        'axes.grid': False,'legend.frameon':True})
sns.set_context("notebook", font_scale=1.3)
mpl.rcParams['savefig.dpi'] = 90

path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(path, 'data')
MIRBASE = os.path.join(datadir, 'miRBase_all.csv')
BOWTIE_INDEXES = None
BOWTIE_PARAMS = '-v 1 --best'
SUBREAD_INDEXES = None
SUBREAD_PARAMS = '-m 2 -M 2'
BWA_INDEXES = None

baseoptions = {'base': [('input',''),('adapter',''),('filetype','fastq'),
                    ('bowtieindex',''),('refgenome',''),('species','hsa'),
                    ('overwrite',1)]}

def write_default_config(conffile='default.conf', defaults={}):
    """Write a default config file"""

    if not os.path.exists(conffile):
        cp = createConfigParserfromDict(defaults, ['base'])
        cp.write(open(conffile,'w'))
        print ('wrote config file %s' %conffile)
    return conffile

def createConfigParserfromDict(data, sections, **kwargs):
    """Helper method to create a ConfigParser from a dict and/or keywords"""

    cp = configparser.ConfigParser()
    for s in sections:
        cp.add_section(s)
        if not data.has_key(s):
            continue
        for i in data[s]:
            name,val = i
            cp.set(s, name, val)
    #use kwargs to create specific settings in the appropriate section
    for s in cp.sections():
        opts = cp.options(s)
        for k in kwargs:
            if k in opts:
                cp.set(s, k, kwargs[k])
    return cp

def parse_config(conffile=None):
    """Parse a configparser file"""

    f = open(conffile,'r')
    cp = configparser.ConfigParser()
    try:
        cp.read(conffile)
    except Exception as e:
        print ('failed to read config file! check format')
        print ('Error returned:', e)
        return
    f.close()
    return cp

def get_options(cp):
    """Makes sure boolean opts are parsed"""

    options = cp._sections['base']
    for o in options:
        try:
            options[o] = cp.getboolean('base', o)
        except:
            pass
    return options

def first(x):
    return x.iloc[0]

def get_exons(gtf):
    """Get exon featues into GenomicArrayOfSets for HTSeq
       counting"""

    exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    for feature in gtf:
        if feature.type == "exon":
            if 'transcript_id' in feature.attr:
                exons[ feature.iv ] += feature.attr["transcript_id"]
            else:
                pass
    return exons

def gtf_to_dataframe(gtf, index='transcript_id'):
    """Convert gtf/gff features to a pandas dataframe"""

    recs=[]
    for f in gtf:
        r = {'id':f.name, 'chrom':f.iv.chrom, 'start':f.iv.start,
             'end': f.iv.end, 'strand':f.iv.strand}
        r.update(f.attr)
        recs.append(r)
    df = pd.DataFrame(recs)
    try:
        df = df.drop_duplicates(['exon_id'])
    except:
        print ('no exon_id field')
    return df

def count_aligned_features(samfile, features, truecounts=None):
    """Count reads in features from an alignment, if no truecounts we
       assume a non-collapsed file was used to map
       Args:
           samfile: mapped sam file
           features: annotation from e.g. gtf file
           truecounts: read counts from original (un-collapsed) file
       Returns: dataframe of genes with total counts
    """

    sam = HTSeq.SAM_Reader(samfile)
    if type(truecounts) is pd.DataFrame:
        truecounts = {r.seq: r['reads'] for i,r in truecounts.iterrows()}
    import collections
    counts = collections.Counter()
    for almnt in sam:
        seq = str(almnt.read)
        if truecounts is not None and seq in truecounts:
            c = truecounts[seq]
        else:
            c = 1
        if not almnt.aligned:
            counts[ "_unmapped" ] += c
            continue
        gene_ids = set()
        for iv, val in features[ almnt.iv ].steps():
            gene_ids |= val
        #print almnt.iv, almnt.read, readcounts[seq], gene_ids
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[ gene_id ] += c
        elif len(gene_ids) == 0:
            counts[ "_no_feature" ] += c
        else:
            counts[ "_ambiguous" ] += c
    result = pd.DataFrame(counts.items(),columns=['name','reads'])
    result['norm'] = result.reads/result.reads.sum()*1e6
    result = result.sort_values(by='reads',ascending=False)

    um = ['_no_feature','_unmapped']
    mapped = float(result[-result.name.isin(um)].reads.sum())
    total = result.reads.sum()
    print ('%s/%s reads counted, %.2f percent' %(mapped, total, mapped/total*100))
    return result

def merge_features(counts, gtf_file):
    """Merge counts with dataframe containing original gtf info"""

    gtf = HTSeq.GFF_Reader(gtf_file)
    df = gtf_to_dataframe(gtf)
    hits = counts.merge(df, left_on='name', right_on='transcript_id', how='left')
    #hits = hits.sort_values(by='reads',ascending=False)
    return hits

def pivot_count_data(df, idxcols=['name']):
    """Pivot read counts over samples containing multiple 'read' columns
       and get mean normalised read counts"""

    if not 'norm' in df.columns:
        df['norm'] = df.fraction*1e6
    x = pd.pivot_table(df, values=['reads','norm'], index=idxcols, columns=['label'])
    x = x.reset_index()
    #print (x[:3])
    x['total_reads'] = x.ix[:,'reads'].sum(1)
    x['mean_norm'] = x.ix[:,'norm'].apply(lambda r: r[r.nonzero()[0]].mean(),1)
    #flatten column index and rename normalised columns
    x.columns = [' '.join(col).strip() for col in x.columns.values]
    x = x.sort_values('mean_norm', ascending=False)
    return x

def get_column_names(df):
    """Get sample column names"""

    cols = [i for i in df.columns if (i.startswith('reads'))]
    ncols = [i for i in df.columns if (i.startswith('norm'))]
    return cols, ncols

def count_aligned(samfile, readcounts=None, by='name', norm_method='total library'):
    """Count short read alignments to any generic fasta index (no features)
       Args:
           samfile: mapped sam file
           readcounts: original read counts if used collapsed reads
           by: whether to group the counts by name (default) or sequence - 'seq'
    """

    sam = HTSeq.SAM_Reader(samfile)
    f=[]
    for a in sam:
        if a.aligned == True:
            f.append((a.read.seq,a.read.name,a.iv.chrom))
        #else:
        #    f.append((a.read.seq,a.read.name,'_unmapped'))

    counts = pd.DataFrame(f, columns=['seq','read','name'])
    if readcounts is not None:
        counts = counts.merge(readcounts, on='seq')
    if by == 'name':
        counts = ( counts.groupby('name')
                  .agg({'reads':np.sum,'seq':first})
                  .reset_index().sort_values('reads',ascending=False) )
    return counts

def map_rnas(files, indexes, outpath, collapse=True, adapters=None, aligner='bowtie',
             use_remaining=False, overwrite=True, verbose=False):
    """Map reads to one or more gene annotations, assumes adapters are removed
    Args:
        files: input fastq read files
        indexes: bowtie indexes of annotations/genomes
        adapters: if adapters need to be trimmed
        bowtieparams: parameters for bowtie
        overwrite: whether to overwrite temp files
    """

    if overwrite == True:
        print ('removing old temp files')
        remove_files(outpath,'*_mapped.sam')
        remove_files(outpath, '*_r.fa')
    cfiles = collapse_files(files, outpath)
    print (cfiles)
    result = []
    for cfile in cfiles:
        rem = None
        label = os.path.splitext(os.path.basename(cfile))[0]
        countfile = os.path.join(outpath, '%s.csv' %label)
        readcounts = pd.read_csv(countfile, index_col=0)
        total = readcounts.reads.sum()
        #print total
        for idx in indexes:
            if use_remaining == True and rem != None:
                query = rem
            else:
                query = cfile
            #print query
            samfile = os.path.join(outpath, '%s_%s.sam' %(label,idx))
            rem = os.path.join(outpath, label+'_r.fa')
            #print samfile
            if aligner == 'bowtie':
                bowtie_align(query, idx, outfile=samfile, remaining=rem, verbose=verbose)
            elif aligner == 'subread':
                subread_align(query, idx, samfile)
            counts = count_aligned(samfile, readcounts)
            counts['label'] = label
            counts['db'] = idx
            counts['fraction'] = counts.reads/total
            #print counts[:2]
            result.append(counts)
    result = pd.concat(result)
    print ('done')
    return result

def map_genome_features(files, ref, gtf_file, outpath='', aligner='bowtie',
                        overwrite=True):
    """Map multiple files to a genome with features and return/process hits.
       Can be used for miRNA discovery
       Args:
           ref: genome bowtie index name
           gtf_file: gtf file with features
           bowtie _index: path with bowtie indexes
    """

    if overwrite == True:
        print ('removing old temp files')
        remove_files(outpath,'*_mapped.sam')
        remove_files(outpath, '*_r.fa')

    gtf = HTSeq.GFF_Reader(gtf_file)
    #use exons for rna-seq
    exons = get_exons(gtf)

    cfiles = collapse_files(files, outpath)
    print (cfiles)
    result = []
    for cfile in cfiles:
        label = os.path.splitext(os.path.basename(cfile))[0]
        samfile = os.path.join(outpath, '%s_%s.sam' %(label,ref))
        if aligner == 'bowtie':
            bowtie_align(cfile, ref, outfile=samfile)
        elif aligner == 'subread':
            subread_align(cfile, ref, samfile)
        #get true read counts for collapsed file
        countfile = os.path.join(outpath, '%s.csv' %label)
        readcounts = pd.read_csv(countfile, index_col=0)
        #count
        hits = count_aligned_features(samfile, exons, readcounts)
        #print hits[:10]
        hits['label'] = label
        hits['genome'] = ref
        result.append(hits)
    result = pd.concat(result)
    return result

def assign_sample_ids(files):
    """Assign ids for each sample and save"""

    l=[]
    i=1
    for f in files:
        label = os.path.splitext(os.path.basename(f))[0]
        sid = 's%02d' %i
        l.append([sid,label,f])
        i+=1
    l = pd.DataFrame(l,columns=['id','label','filename'])
    return l

def remove_files(path, wildcard=''):
    files = glob.glob(os.path.join(path, wildcard))
    for f in files:
        os.remove(f)
    return

def cut_files(files, outpath, adapters):
    """Cut adapters from fastq files"""

    for f in files:
        cut = os.path.join(outpath,label+'.fa')
        if not os.path.exists(cut):
            trim_adapters(f, adapters, cut)
    return

def collapse_reads(infile, outfile=None, min_length=15):
    """Collapse identical reads and retain copy number
      puts all seqs in memory so needs to be optimized"""

    if outfile == None:
        outfile = os.path.splitext(infile)[0]+'_collapsed.fa'
    print ('collapsing reads %s' %infile)
    ext = os.path.splitext(infile)[1]
    if ext == '.fastq':
        fastfile = HTSeq.FastqReader(infile, "solexa")
    elif ext == '.fa' or ext == '.fasta':
        fastfile = HTSeq.FastaReader(infile)
    sequences = [(s.name, s.seq, s.descr) for s in fastfile]
    df = pd.DataFrame(sequences, columns=['id','seq','descr'])
    df['length'] = df.seq.str.len()
    df = df[df.length>=min_length]
    g = df.groupby('seq').agg({'seq':np.size})
    g = g.rename(columns={'seq': 'reads'})
    g = g.sort_values(by='reads',ascending=False).reset_index()
    g['id'] = g.apply(lambda x: 'seq_'+str(x.name),axis=1)
    utils.dataframe_to_fasta(g, outfile=outfile)
    g.to_csv(os.path.splitext(outfile)[0]+'.csv')
    print ('collapsed %s reads to %s' %(len(df),len(g)))
    x = df.length.value_counts()
    return x

def collapse_files(files, outpath, **kwargs):
    """Collapse reads and save counts as csv
        min_length: min length of reads to include
    """

    outfiles = []
    for f in files:
        label = os.path.splitext(os.path.basename(f))[0]
        cfile = os.path.join(outpath, label+'.fa')
        if not os.path.exists(cfile):
            collapse_reads(f, outfile=cfile, **kwargs)
        else:
            print ('found collapsed file')
        outfiles.append(cfile)
    return outfiles

def get_mirbase_sequences(species='hsa'):
    """Extract species specific sequences from mirbase file"""

    df = pd.read_csv(MIRBASE)
    a = df[['mature1','mature1_seq','species']].set_index('mature1').rename(columns={'mature1_seq':'seq'})
    b = df[['mature2','mature2_seq','species']].set_index('mature2').rename(columns={'mature2_seq':'seq'})
    df = pd.concat([a,b]).dropna().reset_index()
    df = df.drop_duplicates('index')
    return df[df.species==species]

def build_mirbase_index(species, kind='mature'):
    """Build species-specific mirbase bowtie index"""

    mirs = get_mirbase_sequences(species)
    print ('got %s sequences' %len(mirs))
    idxname = 'mirbase-'+species
    outfile = '%s.fa' %idxname
    utils.dataframe_to_fasta(mirs, seqkey='seq', idkey='index',
                            outfile=outfile)
    build_bowtie_index(outfile, 'indexes')
    build_subread_index(outfile, 'indexes')
    return idxname

def map_mirbase(files, species='bta', outpath='mirna_results', overwrite=False,
                norm_method=None, aligner='bowtie', **kwargs):
    """Map multiple fastq files to mirbase mature sequences and get
       count results into one file. Used for counting of known miRNAs.
       Species: three letter name of species using mirbase convention
    """

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    #make sample ids
    labels = assign_sample_ids(files)
    print (labels)
    #cut adapters

    #generate new mirbase bowtie index
    if aligner == 'bowtie':
        global BOWTIE_INDEXES
        BOWTIE_INDEXES = 'indexes'
    elif aligner == 'subread':
        global SUBREAD_INDEXES
        SUBREAD_INDEXES = 'indexes'
    db = build_mirbase_index(species)

    #now map to the mirbase index for all files
    res = map_rnas(files, [db], outpath, overwrite=overwrite, aligner=aligner, **kwargs)
    #merge labels with results
    res = res.merge(labels, on='label')
    res.to_csv('mirna_counts.csv')
    return res

def map_mirnas(files, ref, outpath='mirna_results', overwrite=False):
    """Map reads to detect mirnas in genome. This means optional removal
       of other rnas first and then alignment to the genome using aligner
       parameters tuned for reads ~17-25 bp."""

    return

def filter_expr_results(df, freq=0.5, meanreads=0, totalreads=50):
    c,normcols = getColumnNames(df)
    df = df[df.freq>freq]
    df = df[df['total']>=totalreads]
    df = df[df['mean_norm']>=meanreads]
    return df

def compare_expression_profiles(df, by='db', key='reads', threshold=1, path=None):
    """Scatter matrix of count values across samples and/or by label"""

    df = df[df.name!='_unmapped']
    df = df[df[key]>=threshold]
    X = pd.pivot_table(df, values=key, index=['name',by], columns=['label'])
    X = np.log(X).dropna()
    X = X.reset_index(1)
    g = sns.pairplot(X,hue=by)
    if path != None:
        plt.savefig(os.path.join(path, 'expr_scatter.png'))
    return

def get_fractions_mapped(df):
    """Process results of multiple mappings to get fractions
    of each annotations mapped"""

    x = df.groupby(['db','label']).agg({'fraction':np.sum})
    x = x.unstack(level=0)
    x.columns = x.columns.droplevel(0)
    x['unmapped'] = 1-x.sum(1)
    #x = x.T
    #x = x.reindex_axis((x).mean(1).sort_values().index)
    x = x.reset_index()
    return x

def plot_fractions(df, label=None, path=None):
    """Process results of multiple mappings to get fractions
    of each annotations mapped
    label: plot this sample only"""

    if len(df.columns) == 1:
        label = df.columns[0]
    if label != None:
        explode = [0.05 for i in range(len(df))]
        axs = df.plot(y=label,kind='pie',colormap='Spectral',autopct='%.1f%%',startangle=0,figsize=(6,6),
                   labels=None,legend=True,pctdistance=1.1,explode=explode,fontsize=10)
    else:
        l = df.T.plot(kind='bar',stacked=True,cmap='Spectral',figsize=(12,6))
        plt.legend(ncol=4)

    plt.tight_layout()
    if path == None:
        path='.'
    plt.savefig(os.path.join(path,'ncrna_persample.png'))
    return

def bwa_align(infile, ref=None, bowtie_index=None, outfile=None):
    """Align reads with bwa"""

    if bwa_index == None:
        bwa_index = BWA_INDEXES
    ref = os.path.join(bwaindexes, ref)
    label = os.path.splitext(os.path.basename(infile))[0]
    outfile = label+'_'+ref+'_bwa.sam'
    cmd1 = 'bwa aln -n 0 -t 2 %s %s > out.sai' %(ref,infile)
    cmd2 = 'bwa samse %s out.sai %s > %s' %(ref,infile,outfile)
    result = subprocess.check_output(cmd1, shell=True, executable='/bin/bash')
    result = subprocess.check_output(cmd2, shell=True, executable='/bin/bash')
    return

def build_bowtie_index(fastafile, path):
    """Build a bowtie index"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'bowtie-build -f %s %s' %(fastafile, name)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    files = glob.glob(name+'*.ebwt')
    if not os.path.exists(path):
        os.mkdir(path)
    for f in files:
        shutil.move(f, os.path.join(path, f))
        #shutil.move(f, path)
    return

def build_subread_index(fastafile, path):
    """Build an index for subread"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'subread-buildindex -o %s %s' %(name,fastafile)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    if not os.path.exists(path):
        os.mkdir(path)
    exts = ['.00.b.array','.00.b.tab','.files','.reads']
    files = [name+i for i in exts]
    for f in files:
        shutil.move(f, os.path.join(path, f))
    return

def bowtie_align(infile, ref, outfile=None, remaining=None, verbose=True):
    """Map reads using bowtie"""

    label = os.path.splitext(os.path.basename(infile))[0]
    outpath = os.path.dirname(os.path.abspath(infile))
    if outfile == None:
        outfile = label+'_'+ref+'_bowtie.sam'

    if BOWTIE_INDEXES == None:
        print ('base.BOWTIE_INDEXES variable not set')
        return
    os.environ["BOWTIE_INDEXES"] = BOWTIE_INDEXES
    params = BOWTIE_PARAMS
    if remaining == None:
        remaining = os.path.join(outpath, label+'_r.fastq')
    cmd = 'bowtie -f -p 2 -S %s --un %s %s %s > %s' %(params,remaining,ref,infile,outfile)
    if verbose == True:
        print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    if verbose == True:
        print (result)
    return remaining

def subread_align(infile, ref, outfile):
    """Align reads with subread"""

    if SUBREAD_INDEXES == None:
        print ('base.SUBREAD_INDEXES variable not set')
        return
    ref = os.path.join(SUBREAD_INDEXES, ref)
    params = '-t 0 --SAMoutput -T 2 %s' %SUBREAD_PARAMS
    cmd = 'subread-align %s -i %s -r %s -o %s' %(params, ref, infile, outfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    return

def featurecounts(samfile, gtffile):
    """Count aligned features with the featureCounts program.
        Returns: a dataframe of counts"""

    params = '-T 5 -t exon -g transcript_id'
    cmd = 'featureCounts %s -a %s -o counts.txt %s' %(params, gtffile, samfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    counts =  pd.read_csv('counts.txt', sep='\t', comment='#')
    counts = counts.rename(columns={samfile:'reads'})
    counts = counts.sort('reads', ascending=False)
    return counts
