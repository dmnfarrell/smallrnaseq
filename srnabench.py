#!/usr/bin/env python

"""Module for sRNAbench wrappers and utilities
   Created July 2014
   Copyright (C) Damien Farrell
"""

import sys, os, string, types, re
import shutil, glob, collections
import itertools
import subprocess
import pylab as plt
import numpy as np
import pandas as pd
import base

srnabenchoptions = {'base': [('input',''),('adapter','TGGAATTCTCGGGTGCCAAGG'),('filetype','fastq'),
                    ('bowtieindex',''),('refgenome',''),('species','hsa'),
                    ('mature',''), ('hairpin',''), ('other',''),('mirbase',os.getcwd()),
                    ('overwrite',1)]}

def getShortlabel(label):
    x=label.split('_')
    return x[2]+'_'+x[4]

def run(infile, outpath, overwrite=True, adapter=None):
    """Run sRNAbench for a fastq file"""

    label = os.path.splitext(os.path.basename(infile))[0]
    resfile = os.path.join(outpath, 'mature_sense.grouped')
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outdir = os.path.join(outpath, label)
    print 'running %s' %infile
    if os.path.exists(outdir):
        if overwrite == False:
            return None
        else:
            shutil.rmtree(outdir)
    srbpath = '/local/sRNAbench'
    ref = 'bos_taurus_alt'
    #libs = 'bosTau6-tRNAs.fa'
    cmd = ('java -jar %s/sRNAbench.jar dbPath=%s input=%s microRNA=bta' #libs=%s'
           ' species=%s output=%s predict=true plotMiR=true' %(srbpath,srbpath,infile,ref,outdir))
    if adapter != None:
        cmd += ' adapter=%s' %adapter
    print cmd
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print result
    return outdir

def runAll(path, outpath='runs', overwrite=False):
    """Run all fastq files in folder"""

    if os.path.isdir(path):
        files = glob.glob(os.path.join(path,'*.fastq'))
    else:
        files = [path]
    print 'running sRNAbench for %s files' %len(files)
    for f in files:
        res = run(f, outpath, overwrite=overwrite)
        if res == None:
            print 'skipped %s' %f
    return

def readResultsFile(path, infile='mature_sense.grouped', filter=True):
    """Parse .grouped file.
       filter True means remove duplicates of same mirna"""

    f = os.path.join(path,infile)
    #print f
    if not os.path.exists(f):
        return
    df = pd.read_csv(f, sep='\t')
    df['path'] = os.path.basename(path)
    #take highest value if duplicates (same mature)
    g = df.groupby('name').agg(np.max)
    g = g.reset_index()
    return g

def getResults(path, outpath=None):
    """Get single results"""

    if outpath!=None:
        os.chdir(outpath)
    k = readResultsFile(path, 'mature_sense.grouped')
    k['perc'] = k['read count']/k['read count'].sum()
    k = k.sort('read count',ascending=False)
    n = getNovel(path)
    n = n[(n['5pRC']>30) | (n['3pRC']>30)]
    n = n.sort(['5pRC','3pRC'],ascending=False)

    kcols = ['name','read count','unique reads']
    ncols = ['name2','5pRC','3pRC','chrom','chromStart','chromEnd']#,'5pSeq','3pSeq']
    print k[k['read count']>300][kcols]
    print n[ncols]
    fig,ax=plt.subplots(figsize=(8,6))
    ax.set_title('sRNAbench top 10')
    #k.set_index('name')['perc'][:10].plot(kind='pie',ax=ax,
    #                colormap='Set2',autopct='%.2f',startangle=90)
    k.set_index('name')['read count'][:10].plot(kind='barh',colormap='Set2',ax=ax,log=True)
    plt.tight_layout()
    fig.savefig('srnabench_summary_known.png',dpi=80)
    return k,n

def getMultipleResults(path):
    """Fetch results for all dirs and aggregate read counts"""

    k = []
    n = []
    outdirs = [os.path.join(path,i) for i in os.listdir(path)]
    cols = []
    c=1
    for o in outdirs:
        r = readResultsFile(o, 'mature_sense.grouped')
        x = getNovel(o)
        if n is not None:
            n.append(x)
        if r is not None:
            k.append(r)
        cols.append(c)
        c+=1
    print cols
    k = pd.concat(k)
    n = pd.concat(n)
    #combine known into useful format
    p = k.pivot(index='name', columns='path', values='read count')
    g = k.groupby('name').agg({'read count':[np.size,np.mean,np.sum]})
    g.columns = ['freq','mean read count','total']
    g['perc'] = g['total']/g['total'].sum()
    k = p.merge(g,left_index=True,right_index=True)
    k = k.reset_index()
    k = k.sort('mean read count',ascending=False)
    return k,n

def getNovel(path):
    """Parse novel.txt file if available"""

    f = os.path.join(path,'novel.txt')
    if not os.path.exists(f):
        return
    df = pd.read_csv(f, delim_whitespace=True, usecols=range(16), converters={'chrom':str})
    df['path'] = os.path.basename(path)
    return df

def summariseAll(k,n,outpath=None):
    """Summarise multiple results"""

    if outpath != None:
        os.chdir(outpath)
    ky1 = 'unique reads'
    ky2 =  'read count' #'RC'
    cols = ['name','freq','mean read count','total']
    print
    print 'found:'
    final = k[(k['mean read count']>=10) & (k['freq']>=20)]
    print final[cols]
    print '-------------------------------'
    print '%s total' %len(k)
    print '%s with >=5 mean reads' %len(k[k['mean read count']>=5])
    print '%s found in >=5 samples' %len(k[k['freq']>=20])
    print '%s found in 1 sample only' %len(k[k['freq']==1])
    print 'top 10 account for %2.2f' %k['perc'][:10].sum()

    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    k.set_index('name')['total'][:10].plot(kind='barh',colormap='Set2',ax=ax,log=True)
    plt.tight_layout()
    fig.savefig('srnabench_top_known.png')
    fig = plotReadCountDists(final)
    fig.savefig('srnabench_known_counts.png')
    print
    #print n[n.columns[:8]].sort('chromStart')
    n['loc'] = n.apply( lambda x: x.chrom+':'+str(x.chromStart)+'-'+str(x.chromEnd),1)
    #print n.groupby(['loc']).agg({'name':np.size,'5pRC':np.mean,'3pRC':np.mean,
    #                        'chrom':base.first,'chromStart':base.first}) #'3pSeq':base.first,
    return k

def plotReadCountDists(df,h=8):
    """Boxplots of read count distributions per miRNA"""

    w=int(h*(len(df)/60.0))+4
    fig, ax = plt.subplots(figsize=(w,h))
    df = df.set_index('name')
    df = df.drop(['freq','mean read count','total','perc'],1)
    t=df.T
    t.index = df.columns
    #base.sns.boxplot(t,linewidth=1.0,color='coolwarm_r',saturation=0.2,)
    t.plot(kind='box',color='black',grid=False,whis=1.0,ax=ax)
    ax.set_yscale('log')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylabel('read count')
    plt.tight_layout()
    return fig

def runDE(inpath):
    """Call DE routine"""

    return de

def runsrnabenchDE(inpath, l1, l2, cutoff=1.5):
    """DE via sRNABench"""

    files1 = getFilesfromMapping(inpath, l1)
    files2 = getFilesfromMapping(inpath, l2)
    files1 = [os.path.basename(i) for i in files1]
    files2 = [os.path.basename(i) for i in files2]

    lbl = 'mature' #'hairpin'
    output = 'DEoutput_srnabench'
    if os.path.exists(output):
        shutil.rmtree(output)
    grpstr = ':'.join(files1)+'#'+':'.join(files2)
    cmd = ('java -jar /local/sRNAbench/sRNAbenchDE.jar input=%s '
           'grpString=%s output=%s diffExpr=true minRC=1 readLevel=true seqStat=true '
           'diffExprFiles=%s_sense.grouped' %(inpath,grpstr,output,lbl))
    print l1,l2
    #print cmd
    print
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    resfile = glob.glob(os.path.join(output,'*.edgeR'))[0]
    print lbl
    df = pd.read_csv(resfile, sep='\t')
    df = df.reset_index()
    df = df.rename(columns={'index':'name'})
    df = df.sort('logFC',ascending=False)
    df = df[(df.FDR<0.05) & ((df.logFC>cutoff) | (df.logFC<-cutoff))]
    return df

def test():
    runDE('results_srnabench_combined', 'MAP TP0', 'MAP TP END')
    runsrnabenchDE('results_srnabench_combined', 'MAP TP0', 'MAP TP END')
    return

def main():
    #base.seabornsetup()
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-r", "--run", dest="run", action='store_true',
                           help="run predictions")
    parser.add_option("-i", "--input", dest="input",
                           help="input path or file")
    parser.add_option("-s", "--summarise", dest="summarise",
                           help="analyse results of multiple runs together")
    parser.add_option("-a", "--analyse", dest="analyse",
                           help="analyse results of single run")
    opts, remainder = parser.parse_args()
    pd.set_option('display.width', 800)
    if opts.run == True:
        runAll(opts.input)
    elif opts.summarise != None:
        k,n = getMultipleResults(opts.summarise)
        summariseAll(k,n)
    elif opts.analyse != None:
        getResults(opts.analyse)
    else:
        test()

if __name__ == '__main__':
    main()
