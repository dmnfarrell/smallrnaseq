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

try:
    import HTSeq
except:
    'HTSeq not present'

from . import utils, aligners
import matplotlib as mpl
import seaborn as sns
sns.set_style("ticks", {'axes.facecolor': '#F7F7F7',
                        'axes.grid': False,'legend.frameon':True})
sns.set_context("notebook", font_scale=1.3)
mpl.rcParams['savefig.dpi'] = 90

path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(path, 'data')
MIRBASE = os.path.join(datadir, 'miRBase_all.csv')


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

def gtf_to_dataframe(gtf=None, gtf_file=None, index='transcript_id'):
    """Convert gtf/gff features to a pandas dataframe"""

    if gtf_file != None:
        gtf = HTSeq.GFF_Reader(gtf_file)
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

def count_features(samfile, features=None, gtffile=None, readcounts=None, merge=False):
    """Count reads in features from an alignment, if no truecounts we
       assume a non-collapsed file was used to map
       Args:
           samfile: mapped sam file
           gtffile: feature file
           features: annotations read from bed or gtf file
           truecounts: read counts from original (un-collapsed) file
           merge: whether to merge the gtf fields with the results
       Returns: dataframe of genes with total counts
    """

    if gtffile != None:
        gtf = HTSeq.GFF_Reader(gtffile)
        features = get_exons(gtf)
    sam = HTSeq.SAM_Reader(samfile)
    if type(readcounts) is pd.DataFrame:
        readcounts = {r.seq: r['reads'] for i,r in readcounts.iterrows()}
    import collections
    counts = collections.Counter()
    for almnt in sam:
        seq = str(almnt.read)
        if readcounts is not None and seq in readcounts:
            c = readcounts[seq]
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
    #result['norm'] = result.reads/result.reads.sum()*1e6
    result = result.sort_values(by='reads',ascending=False)

    um = ['_no_feature','_unmapped']
    mapped = float(result[-result.name.isin(um)].reads.sum())
    total = result.reads.sum()

    print ('%s/%s reads counted, %.2f percent' %(mapped, total, mapped/total*100))
    if merge == True and gtffile != None:
        result = merge_features(result, gtffile)
    return result

def merge_features(counts, gtffile):
    """Merge counts with dataframe containing original gtf info"""

    gtf = HTSeq.GFF_Reader(gtffile)
    df = gtf_to_dataframe(gtf)
    hits = counts.merge(df, left_on='name', right_on='transcript_id', how='inner')
    #hits = hits.sort_values(by='reads',ascending=False)
    return hits

def feature_counts_summary(counts):
    """Summary of feature counts by gene biotype"""

    counts['perc'] = counts.norm/1e6*100
    s = counts.groupby('gene_biotype').agg({'reads':np.sum,'perc':np.sum})
    ax = s.plot(y='perc',kind='barh',figsize=(5,3))
    return s

def get_top_genes(counts):

    df = counts.groupby('gene_name')\
                .agg({'reads':sum,'transcript_id':np.size})\
                .sort_values('reads',ascending=False)
    return df

def count_aligned(samfile, collapsed=None, readcounts=None, by='name', duplicates=False):
    """Count short read alignments from a sam or bam file. Designed to be used with
       collapsed reads with original read counts in fasta id.
       Args:
           samfile: mapped sam file
           collapsed: collapsed fasta with original read counts in name
           readcounts: dataframe with original read counts, optional
           by: whether to group the counts by 'name' (default) or 'seq'
    """

    if collapsed != None:
        readcounts = utils.read_collapsed_file(collapsed)

    counts = utils.get_aligned_reads(samfile, readcounts=readcounts)
    if duplicates == False:
        counts = counts.drop_duplicates('seq')
    total = readcounts.reads.sum()
    counts = counts.groupby('name').agg({'reads':np.sum}).reset_index()
    mapped = float(counts.reads.sum())
    if len(counts) > 0:
        print ('%s/%s reads counted, %.2f percent' %(mapped, total, mapped/total*100))
    else:
        print ('no counts found')
    return counts

def pivot_count_data(counts, idxcols='name', norm_method='library'):
    """Pivot read counts created by count_aligned over multiple samples
       and get normalised read counts.
       Args:
           counts: dataframe of raw count data with samples per column
           idxcols: name of index column
           norm_method: how to normalize the counts (returned in extra columns),
                        default is total library counts
       Returns: dataframe of raw /normalised read counts with column per sample
    """

    x = pd.pivot_table(counts, values='reads', index=idxcols, columns='label')
    #print assign_sample_ids(x.columns)
    if norm_method == 'library':
        n = total_library_normalize(x)
    elif norm_method == 'quantile':
        n = quantile_normalize(x)

    #x.columns = [i+'_' for i in x.columns]
    scols = x.columns
    ncols = n.columns = [i+' norm' for i in n.columns]
    x = x.join(n)
    #scols,ncols = get_column_names(x)
    x['total_reads'] = x[scols].sum(1)
    x['mean_norm'] = x[ncols].apply(lambda r: r[r.nonzero()[0]].mean(),1)
    x = x.reset_index()
    return x

def normalize_samples(counts, norm_method='library', rename=True):
    """Normalize over a matrix of samples explicitly, this will overwrite any 'norm'
       columns created previously when pivoting the count data
       Args:
            counts: dataframe of raw count data with samples per column
            rename: rename columns with 'norm' label and add to existing ones
       Returns: dataframe of raw /normalised read counts
    """

    x = counts
    if norm_method == 'library':
        n = total_library_normalize(x)
    elif norm_method == 'quantile':
        n = quantile_normalize(x)
    if rename == True:
        scols = x.columns
        ncols = n.columns = [i+' norm' for i in n.columns]
        x = x.join(n)
    else:
        x = n
    return x

def get_column_names(df):
    """Get count data sample column names"""

    ignore = ['total_reads','mean_norm']
    ncols = [i for i in df.columns if (i.endswith('norm')) and i not in ignore]
    cols = [i.split(' ')[0] for i in ncols if i not in ignore]
    return cols, ncols

def total_library_normalize(df):
    """Normalise by size of total reads"""

    df = df.copy()
    for col in df:
        df[col] = df[col]/df[col].sum()*1e6
    df = df.round(2)
    return df

def quantile_normalize(df):
    """Quantile normlization of counts for multiple samples.
       see https://github.com/ShawnLYU/Quantile_Normalize
    """

    df = df.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    df=df.round(2)
    return df

def upper_quartile_normalize(df):
    """Upper quartile noralisation"""

    u = x.apply(lambda r: r[r>r.quantile(0.75)])
    df = df/u.mean()
    return df*1e6

def deseq_normalize(df):
    """Compute count/geometric mean per sample"""

    df = df.copy()
    df = df.apply(lambda r: r/r.mean(),1)
    return df

def map_rnas(files, indexes, outpath, collapse=True, adapters=None, aligner='bowtie',
             norm_method='library', use_remaining=True, overwrite=False,
             samplelabels=None, params={}):
    """Map reads to one or more gene annotations, assumes adapters are removed
    Args:
        files: input fastq read files
        indexes: bowtie indexes of annotations/genomes
        adapters: if adapters need to be trimmed
        use_remaining: only align to remaining reads after each index
        overwrite: whether to overwrite temp files
        samplelabels: mapping of file names to short ids
        params: index specific aligner parameters
    """

    if not os.path.exists(aligners.BOWTIE_INDEXES):
        print ('index folder does not exist!')
        return
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    if overwrite == True:
        print ('removing old temp files')
        utils.remove_files(outpath,'*_mapped.sam')
        utils.remove_files(outpath, '*_r.fa')
    if collapse == True:
        cfiles = collapse_files(files, outpath)
    else:
        cfiles = files
    if len(cfiles)==0:
        print ('WARNING no files to align')
        return
    print (cfiles)
    result = []
    #store current aligner parameters as default
    default_params = aligners.get_current_params(aligner)

    for cfile in cfiles:
        rem = None
        filename = os.path.splitext(os.path.basename(cfile))[0]
        readcounts = utils.read_collapsed_file(cfile)
        total = readcounts.reads.sum()
        print (filename)
        for idx in indexes:
            print (idx)
            if idx in params:
                align_params = params[idx]
            else:
                align_params = default_params
            aligners.set_params(aligner, align_params)
            #print (align_params)
            if use_remaining == True and rem != None:
                query = rem
            else:
                query = cfile
            samfile = os.path.join(outpath, '%s_%s.sam' %(filename,idx))
            rem = os.path.join(outpath, filename+'_r.fa')

            if aligner == 'bowtie':
                aligners.bowtie_align(query, idx, outfile=samfile,
                                      remaining=rem, verbose=True)
            elif aligner == 'subread':
                aligners.subread_align(query, idx, samfile)

            counts = count_aligned(samfile, readcounts=readcounts)
            if len(counts) == 0:
                print ('WARNING: no counts found for %s.' %idx)
                continue
            if samplelabels != None:
                counts['label'] = samplelabels[filename]
            else:
                counts['label'] = filename
            counts['ref'] = idx
            counts['fraction'] = counts.reads/total
            #print (counts.reads.sum(), total)
            result.append(counts)
            print()
        print()
    if len(result) == 0:
        return
    result = pd.concat(result)
    counts = pivot_count_data(result, idxcols=['name','ref'], norm_method=norm_method)
    print ('done')
    return result, counts

def map_genome_features(files, ref, gtf_file, outpath='', aligner='bowtie',
                        overwrite=True):
    """Map multiple files to a genome with features and return/process hits.
       Can be used for miRNA discovery
       Args:
           ref: genome bowtie index name
           gtf_file: gtf or bed file with features
           bowtie _index: path with bowtie indexes
    """

    if overwrite == True:
        print ('removing old temp files')
        utils.remove_files(outpath,'*_mapped.sam')
        utils.remove_files(outpath, '*_r.fa')

    ext = os.path.splitext(gtf_file)[1]
    if ext == '.gtf' or ext == '.gff':
        features = HTSeq.GFF_Reader(gtf_file)
    elif ext == '.bed':
        features = HTSeq.BED_Reader(gtf_file)

    exons = get_exons(features)

    cfiles = collapse_files(files, outpath)
    print (cfiles)
    result = []
    for cfile in cfiles:
        label = os.path.splitext(os.path.basename(cfile))[0]
        samfile = os.path.join(outpath, '%s_%s.sam' %(label,ref))
        if aligner == 'bowtie':
            aligners.bowtie_align(cfile, ref, outfile=samfile)
        elif aligner == 'subread':
            aligners.subread_align(cfile, ref, samfile)
        #get true read counts for collapsed file
        countfile = os.path.join(outpath, '%s.csv' %label)
        readcounts = pd.read_csv(countfile, index_col=0)
        #count features
        hits = count_features(samfile, features=exons, truecounts=readcounts)
        hits['label'] = label
        hits['genome'] = ref
        result.append(hits)
    result = pd.concat(result)
    result = merge_features(result, gtf_file)
    return result

def get_base_names(files):
    names = [os.path.splitext(os.path.basename(f))[0] for f in files]
    return names

def assign_sample_ids(names, outfile='sample_labels.csv'):
    """Assign new ids for sample filenames e.g. files. Useful for
       replacing long file names with short ids.
       Returns: dict of filename/id values
    """

    i=1
    labels = {}
    for fname in names:
        n = os.path.splitext(os.path.basename(fname))[0]
        print (n)
        sid = 's%02d' %i
        labels[n] = sid
        i+=1
    l = pd.DataFrame.from_dict(labels,orient='index')
    l.columns = ['id']; l.index.name='filename'
    l.to_csv(outfile)
    return labels

def trim_files(files, outpath, adapters):
    """Trim adapters from fastq files"""

    for f in files:
        cut = os.path.join(outpath, f)
        if not os.path.exists(cut):
            trim_adapters(f, adapters, cut)
    return

def collapse_reads(infile, outfile=None, min_length=15, progress=False):
    """Collapse identical reads, retaining copy number in a csv file
       and writing collapsed reads to a new fasta file.
       Creates a collapsed fasta file of unique reads and a csv
       file with the """

    from itertools import islice
    if outfile == None:
        outfile = os.path.splitext(infile)[0]+'_collapsed.fa'
    print ('collapsing reads %s' %infile)
    ext = os.path.splitext(infile)[1]
    if ext == '.fastq':
        fastfile = HTSeq.FastqReader(infile, "solexa")
    elif ext == '.fa' or ext == '.fasta':
        fastfile = HTSeq.FastaReader(infile)
    else:
        print ('not fasta or fastq')
        return False

    i=0
    total = 0
    f = {}
    for s in fastfile:
        seq = s.seq.decode()
        if seq in f:
            f[seq]['reads'] += 1
        else:
            f[seq] = {'name':s.name, 'reads':1}
        total += 1

    df = pd.DataFrame.from_dict(f, orient='index')
    df.index.name = 'seq'
    df = df.reset_index()
    l = df.seq.str.len()
    df = df[l>=min_length]
    df = df.drop(['name'], 1)
    df = df.sort_values(by='reads', ascending=False).reset_index()
    df['read_id'] = df.index.copy()
    df['read_id'] = df.apply( lambda x: str(x.read_id)+'_'+str(x.reads), 1 )
    #print df[:10]
    utils.dataframe_to_fasta(df, idkey='read_id', outfile=outfile)
    #df.to_csv(os.path.splitext(outfile)[0]+'.csv', index=False)
    print ('collapsed %s reads to %s' %(total,len(df)))
    return True

def collapse_files(files, outpath, **kwargs):
    """Collapse reads and save counts
        min_length: min length of reads to include
    """

    outfiles = []
    for f in files:
        label = os.path.splitext(os.path.basename(f))[0]
        collapsedfile = os.path.join(outpath, label+'.fa')
        #countsfile = os.path.join(outpath, label+'.csv')
        if not os.path.exists(collapsedfile):# or not os.path.exists(countsfile):
            res = collapse_reads(f, outfile=collapsedfile, **kwargs)
            if res == False:
                continue
        else:
            print ('found collapsed file')
        outfiles.append(collapsedfile)
    return outfiles

def _get_mature(r, idkey='mature1', seqkey='mature1_seq', pad5=0, pad3=0):
    """Get mature sequences from mirbase file,
      this is a row-based dataframe function"""

    p = r.precursor
    name = r[idkey]
    m = r[seqkey]
    if pd.isnull(m):
        s = np.nan
    else:
        i = p.find(m)
        start = i-pad5
        if start < 0: start = 0
        s = p[start:i+len(m)+pad3]
    return pd.Series([name,s],index=['name','sequence'])

def get_mature_padded(df, **kwargs):

    m = df.apply(lambda x: _get_mature(x,**kwargs), 1)
    return m

def get_mirbase(species=None):
    """Get mirbase sequences"""

    df = pd.read_csv(MIRBASE)
    if species is not None:
        df = df[df.species==species]
    df.precursor = df.precursor.str.replace('U','T')
    df.mature1_seq = df.mature1_seq.str.replace('U','T')
    df.mature2_seq = df.mature2_seq.str.replace('U','T')
    return df

def get_mirbase_sequences(species='hsa', pad5=0, pad3=0, dna=False):
    """Extract species specific sequences from mirbase file.
       Args:
           species: 3-letter code for species
           n: bases to extend around ends of mature sequence
       Returns:
            dataframe with mature or hairpin sequences
    """

    df = pd.read_csv(MIRBASE)
    if species != None:
        df = df[df.species==species]

    #get both 5p and 3p seqs for each mirna
    m1 = df.apply(lambda x: _get_mature(x, 'mature1', 'mature1_seq', pad5, pad3), 1)
    m2 = df.apply(lambda x: _get_mature(x, 'mature2', 'mature2_seq', pad5, pad3), 1)
    df = pd.concat([m1,m2]).dropna().reset_index(drop=True)
    df = df[df.sequence.str.len()>2]
    df = df.drop_duplicates('name')
    if dna == True:
        df['sequence'] = df.sequence.str.replace('U','T')
    return df

def build_mirbase_index(species, aligner='bowtie', pad5=3, pad3=5,
                        kind='mature', index_path='indexes'):
    """Build species-specific mirbase bowtie index
       Args:
           species: 3-letter code for species
           aligner: which aligner to build for
           n: bases to extend around ends of mature sequence
    """

    if kind == 'mature':
        mirs = get_mirbase_sequences(species, pad5, pad3)
        idxname = 'mirbase-'+species
        outfile = '%s.fa' %idxname
        utils.dataframe_to_fasta(mirs, seqkey='sequence', idkey='name',
                                outfile=outfile)
    else:
        mirs = get_mirbase(species)
        idxname = 'precursors-'+species
        outfile = '%s.fa' %idxname
        utils.dataframe_to_fasta(mirs, seqkey='precursor', idkey='mirbase_id',
                                outfile=outfile)
        print ('got %s sequences' %len(mirs))

    if aligner == 'bowtie':
        aligners.build_bowtie_index(outfile, index_path)
    elif aligner == 'subread':
        aligners.build_subread_index(outfile, index_path)
    return idxname

def map_mirbase(files, species='bta', outpath='mirna_results', indexes=[],
                index_path='indexes', ref_genome='',
                overwrite=False, aligner='bowtie',
                pad5=3, pad=5, **kwargs):
    """Map multiple fastq files to mirbase mature sequences and get
       count results into one file. Used for counting of known miRNAs.
       Args:
           Species: three letter name of species using mirbase convention
           indexes: other libraries to align to before we count mirnas
    """

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    #generate new mirbase bowtie index
    if aligner == 'bowtie':
        aligners.BOWTIE_INDEXES = index_path
        if aligners.BOWTIE_PARAMS == None:
            aligners.BOWTIE_PARAMS = '-n 1 -l 20'
    elif aligner == 'subread':
        aligners.SUBREAD_INDEXES = index_path
        #SUBREAD_PARAMS = '-m 2 -M 2'
    midx = build_mirbase_index(species, aligner, pad)
    pidx = build_mirbase_index(species, aligner, kind='precursor')

    indexes.extend([midx, pidx])
    #add the reference genome to be aligned last if provided
    if ref_genome != '':
        indexes.append(ref_genome)
    #now map to the mirbase index for all files
    res, counts = map_rnas(files, indexes, outpath, overwrite=overwrite,
                           aligner=aligner, **kwargs)
    #output_read_stacks(files, midx)
    return res, counts

def map_isomirs(files, outpath, species, samplelabels=None):
    """Count mirna isomirs using previously aligned files"""

    idx = 'mirbase-'+species
    result = []
    for f in files:
        filename = os.path.splitext(os.path.basename(f))[0]
        samfile = os.path.join(outpath, '%s_%s.sam' %(filename,idx))
        #countsfile = os.path.join(outpath, '%s.csv' %filename)
        collapsed = os.path.join(outpath, '%s.fa' %filename)
        c = count_isomirs(samfile, collapsed, species)
        if samplelabels != None:
            c['label'] = samplelabels[filename]
        else:
            c['label'] = filename
        result.append(c)
    result = pd.concat(result)
    counts = pivot_count_data(result, idxcols=['name'])
    return result, counts

def count_isomirs(samfile, collapsed, species):
    """Count miRNA isomirs using aligned reads from a samfile and actual
       read counts from a csv file"""

    print (samfile, collapsed)
    canonical = get_mirbase_sequences(species, dna=True).set_index('name')
    #padded sequences so we can see where each read landed relative to canonical
    mirs = get_mirbase_sequences(species, pad5=6, pad3=6, dna=True).set_index('name')
    reads = utils.get_aligned_reads(samfile, collapsed)
    reads = reads.drop(['read_id','start','end'],1)
    return reads

def output_read_stacks(files, ref_fasta, outfile='read_stacks.txt'):
    """Output read stack files from sam alignment files and a ref sequence"""

    out = open(outfile, 'w')
    s=''
    for f in files:
        filename = os.path.splitext(os.path.basename(f))[0]
        samfile = os.path.join(outpath, '%s_%s.sam' %(filename,idx))
        countfile = os.path.join(outpath, '%s.csv' %filename)
        readcounts = pd.read_csv(countfile)
        reads = utils.get_aligned_reads(samfile, readcounts)
        s += filename+'\n'
        s += utils.print_read_stacks(reads, fastafile='mirbase-bta.fa')
    #print (s)
    out.write(s)
    out.close()
    return

def filter_expr_results(df, freq=0.5, meanreads=0, totalreads=50):

    c,normcols = getColumnNames(df)
    df = df[df.freq>freq]
    df = df[df['total']>=totalreads]
    df = df[df['mean_norm']>=meanreads]
    return df

def compare_expression_profiles(df, by='ref', key='reads', threshold=1, path=None):
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

    x = df.groupby(['ref','label']).agg({'fraction':np.sum})
    x = x.unstack(level=0)
    x.columns = x.columns.droplevel(0)
    x['unmapped'] = 1-x.sum(1)
    #x = x.T
    #x = x.reindex_axis((x).mean(1).sort_values().index)
    x = x.reset_index()
    return x
