#!/usr/bin/env python

"""Misc miRNA analysis routines
   Created June 2014
   Copyright (C) Damien Farrell
"""

import sys, os, string, types, re
import shutil, glob, collections
import itertools
import subprocess
import matplotlib
import pylab as plt
import numpy as np
import pandas as pd
import HTSeq
import base

def first(x):
    return x.iloc[0]

def readLengthDist(df,ax=None):

    df['length'] = df.seq.str.len()
    bins=np.linspace(1,df.length.max(),df.length.max())
    if ax!=None:
        df.hist('length',bins=bins,ax=ax)
        plt.title('read length distribution')
        plt.xlabel('length')
    return np.histogram(df.length,bins=bins)

def fastq2fasta(infile, rename=True):

    fastqfile = HTSeq.FastqReader(infile, "solexa")
    outfile = open(os.path.splitext(infile)[0]+'.fa','w')
    i=1
    for s in fastqfile:
        if rename==True:
            s.name=str(i)
        s.write_to_fasta_file(outfile)
        i+=1
    outfile.close()
    return

def summariseFastq(f, filetype='fastq'):

    if filetype=='fastq':
        ffile = HTSeq.FastqReader(f, "solexa")
    else:
        ffile = HTSeq.FastaReader(f)
    sequences = [(s.name,s.seq) for s in ffile]
    #df = pd.DataFrame(sequences,columns=['id','seq'])
    return df

def summariseReads(path):
    """Count reads in all files in path"""

    resultfile = 'read_stats.csv'
    files = glob.glob(os.path.join(path,'*.fastq'))
    vals=[]
    for f in files:
        label = os.path.splitext(os.path.basename(f))[0]
        s = summariseFastq(f)
        l=len(s)
        vals.append([label,l])
        print label, l
    df = pd.DataFrame(vals,columns=['path','total reads'])
    df.to_csv(resultfile)
    df.barh('total reads')
    plt.xlabel('total read count')
    plt.ylabel('sample')
    plt.savefig('read_stats.png')
    return

def collapseReads(infile, outfile='collapsed.fa'):
    """Collapse identical reads and retain copy number
      - may use a lot of memory"""

    print 'collapsing reads %s' %infile
    fastqfile = HTSeq.FastqReader(infile, "solexa")
    #fastafile = HTSeq.FastaReader(infile)
    sequences = [(s.name, s.seq, s.descr) for s in fastqfile]
    df = pd.DataFrame(sequences, columns=['id','seq','descr'])
    df['length'] = df.seq.str.len()
    g = df.groupby('seq').agg({'seq':np.size})
    g = g.rename(columns={'seq': 'descr'}).reset_index()
    g = g.sort('descr',ascending=False)
    g['id'] = g.apply(lambda x: 'seq_'+str(x.name),axis=1)
    base.dataframe2Fasta(g,outfile=outfile)
    g.to_csv(os.path.splitext(outfile)[0]+'.csv')
    print 'collapsed %s reads to %s' %(len(df),len(g))
    #bins=np.linspace(1,df.length.max(),df.length.max())
    #x=np.histogram(df.length,bins=bins)
    #x=pd.Series(l[0],index=l[1][1:])
    x=df.length.value_counts()
    return x

def trimAdapters(infile, adapters, outfile='cut.fastq'):
    """Trim adapters using cutadapt"""

    if os.path.exists(outfile):
        return
    adptstr = ' -a '.join(adapters)
    cmd = 'cutadapt -m 15 -O 5 -q 20 --discard-untrimmed -a %s %s -o %s' %(adptstr,infile,outfile)
    print cmd
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    #print result
    return

def removeKnownRNAs(path, outpath='data_RNAremoved'):
    """Map to annotated RNAs and put remaining reads in output dir"""

    index = 'bosTau6-tRNAs'
    params = '-v 0 --best'
    files = glob.glob(os.path.join(path,'*.fastq'))
    #files = ['test.fastq']
    for f in files:
        print f
        label = os.path.splitext(os.path.basename(f))[0]
        #trimAdapters(f)
        fastq2fasta('cut.fastq')
        rem = os.path.join(outpath, label+'.fa')
        #samfile = os.path.join(outpath, '%s_%s_mapped.sam' %(label,index))
        samfile = 'out.sam'
        base.bowtieMap('cut.fa', index, outfile=samfile, params=params, remaining=rem)

    return

def mapRNAs(files=None, path=None, indexes=[], adapters=None):
    """Map to various ncRNA annotations and quantify perc of reads mapping.
        The order of indexes will affect results.
        path: input path with read files
        indexes: bowtie indexes of annotated rna classes"""

    params = '-v 1 --best'
    if files == None:
        files = glob.glob(os.path.join(path,'*.fastq'))
    outpath = 'ncrna_map'
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    outfiles = []
    for f in files:
        print f
        label = os.path.splitext(os.path.basename(f))[0]
        cut = os.path.join(outpath,label+'_cut.fastq')
        cfile = os.path.join(outpath,label+'.fa')
        print cut,cfile
        if not os.path.exists(cut) and adapters!=None:
            trimAdapters(f, adapters, cut)
        else:
            cut = f
        collapseReads(cut, outfile=cfile)
        outfiles.append(cfile)

    #cfiles = glob.glob(os.path.join(outpath,'*.fa'))
    res=[]
    colnames = ['name','total']+indexes
    readsleft = []
    for cfile in outfiles:
        label = os.path.splitext(os.path.basename(cfile))[0]
        #get total reads by using copy no. of each unique seq
        counts = pd.read_csv(os.path.join(outpath, '%s.csv' %label))
        total = counts['descr'].sum()
        x=[label,total]
        rem = None
        #fig,axs = plt.subplots(5,1)
        #grid=axs.flat
        #bins=np.arange(15,45,1.)
        i=0
        for index in indexes:
            if rem != None:
                query = rem
            else:
                query = cfile
            rem = os.path.join(outpath, label+'_r.fastq')
            samfile = os.path.join(outpath, '%s_%s_mapped.sam' %(label,index))
            if not os.path.exists(samfile):
                print rem
                rem = base.bowtieMap(query, index, outfile=samfile, params=params, remaining=rem)
            sam = HTSeq.SAM_Reader(samfile)
            f = [(a.read.seq,a.read.name) for a in sam if a.aligned == True]
            if len(f)>0:
                found = pd.DataFrame(f, columns=['seq','name'])
                #get original counts for mapped reads and sum them
                nr = counts[counts.id.isin(found.name)]
                fc = nr['descr'].sum()
                perc = fc/float(total)
                print index, len(f), fc, total, round(perc,3)
                #found['length'] = found.seq.str.len()
                #found.length.hist(ax=grid[i],bins=bins)
                #grid[i].set_title(index)
                #plt.xlim((15,45))
                #plt.tight_layout()
                #i+=1
            else: perc = 0.0
            print '%s: %.4f of total reads aligned' %(index,perc)
            x.append(perc)
        res.append(x)

    df = pd.DataFrame(res,columns=colnames)
    df.to_csv('ncrna_mapped.csv',float_format='%.3f')
    #df = pd.read_csv('ncrna_mapped.csv')
    df=df.sort('total')
    df = df.drop('total',1)
    df['unmapped'] = 1-df.sum(1)
    df['name'] = df['name'].apply(lambda r: ' '.join(r.split('_')[2:5]))
    df = df.set_index('name')
    print df

    plt.figure(figsize=(8,8))
    df.mean().plot(kind='pie',colormap='Paired',autopct='%.2f',startangle=90,
                    labels=None,legend=True)
    plt.title('mean percentage small RNAs mapped by category')
    plt.tight_layout()
    plt.savefig('ncrna_means.png',dpi=80)

    '''l = mp.plot(kind='bar',stacked=True,cmap='Paired',figsize=(12,6))
    plt.xlabel('percent mapped')
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig('ncrna_bysample.png')

    plt.figure(figsize=(12,8))
    ax=base.sns.boxplot(mp,color='Paired')
    plt.ylabel('percent mapped')
    plt.title('mean percentage small RNAs mapped by category')
    plt.tight_layout()
    plt.savefig('ncrna_dists.png')

    ax=df[df.pool==1].plot('total','bosTau6-tRNAs',kind='scatter',s=40,color='blue')
    df[df.pool==2].plot('total','bosTau6-tRNAs',kind='scatter',s=40,color='red',ax=ax)
    plt.title('percentage mapped vs total reads')
    plt.legend(['pool 1','pool 2'])
    plt.xlabel('total reads')
    plt.savefig('ncrna_mappedvreads.png')'''
    plt.show()
    return

def compareMethods():
    """Compare 2 methods for subsets of samples"""

    path1 = 'results_mirdeep_rnafiltered'
    path2 = 'results_srnabench_rnafiltered'

    #compare means of filtered knowns
    df = mdp.getResults(path1)
    df = df[df.novel==False]
    mk = mdp.filterExprResults(df,meanreads=200,freq=0.8)
    k,n = srb.getResults(path2)
    #sk = k[(k['mean read count']>=10) & (k['freq']>=0.8)]
    sk = k[(k['read count']>=500)]
    x = pd.merge(mk,sk,left_on='#miRNA',right_on='name',how='inner',suffixes=['1','2'])
    fig = plt.figure(figsize=(12,6))
    ax=fig.add_subplot(121)
    base.venndiagram([mk['#miRNA'], sk['name']],['mirdeep2','srnabench'],ax=ax)
    #print sk[-sk.name.isin(mk['#miRNA'])][['name','total']]
    ax=fig.add_subplot(122)
    x.plot('total','read count',kind='scatter',ax=ax, logx=True,logy=True,alpha=0.8,s=40)
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", color='red')
    ax.set_xlabel('mirdeep2')
    ax.set_ylabel('srnabench')
    ax.set_title('total read count comparison')
    plt.tight_layout()
    plt.savefig('mirdeep_vs_srnabench.png')
    plt.show()
    return

def test():
    path = '/opt/mirnaseq/data/vegh_13'
    #files = ['/opt/mirnaseq/analysis/test.fastq']
    files = ['/opt/mirnaseq/data/vegh_13/SRR576286.fastq']
    bidx =  ['mirbase-mature','Rfam_btau','bosTau6-tRNAs','noncodev4_btau','bos_taurus_alt']
    #adapters for our data
    adapters = ['TGGAATTCTCGGGTGCCAAGG','GCATTGTGGTTCAGTGGTAGAATTCTCGC']
    #adapters for Vegh data
    adapters = ['TAGCTTATCAGACTGATGTTGA','AGATCGGAAGAGCACACGTCTGAACTCC']
    #mapRNAs(files=files, indexes=bidx, adapters=adapters)
    summariseReads(path)
    return

if __name__ == '__main__':
    test()
