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

srbpath = '/local/sRNAbench'
srnabenchoptions = {'base': [('input',''),('outpath','srnabench_runs'),
                    ('adapter','TGGAATTCTCGGGTGCCAAGG'),('filetype','fastq'),
                    ('bowtieindex',''),('refgenome',''),('species','hsa'),
                    ('mature',''), ('hairpin',''), ('other',''), ('isomir','false'),
                    ('overwrite',1), ('matureMM',0), ('p',3)]}

def getShortlabel(label):
    x=label.split('_')
    return x[2]+'_'+x[4]

def run(infile, outpath='srnabench_runs', overwrite=True, adapter=None,
            ref='bos_taurus_alt', predict='false', **kwargs):
    """Run sRNAbench for a fastq file"""

    label = os.path.splitext(os.path.basename(infile))[0]
    resfile = os.path.join(outpath, 'mature_sense.grouped')
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outdir = os.path.join(outpath, label)
    print 'running %s' %infile

    if os.path.exists(outdir):
        if overwrite == False:
            return outdir
        else:
            shutil.rmtree(outdir)

    cmd = ('java -jar %s/sRNAbench.jar dbPath=%s input=%s microRNA=bta'
           ' species=%s output=%s predict=%s plotMiR=true matureMM=0 isoMiR=true'
           ' p=3' %(srbpath,srbpath,infile,ref,outdir,predict))
    if adapter != None:
        cmd += ' adapter=%s' %adapter
    print cmd
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print result
    return outdir

def runAll(path, outpath='runs', overwrite=False, filetype='fastq'):
    """Run all fastq files in folder"""

    if os.path.isdir(path):
        files = glob.glob(os.path.join(path,'*.'+filetype))
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

def plotResults(k):
    fig,ax=plt.subplots(figsize=(8,6))
    ax.set_title('sRNAbench top 10')
    k.set_index('name')['read count'][:10].plot(kind='barh',colormap='Set2',ax=ax,log=True)
    plt.tight_layout()
    fig.savefig('srnabench_summary_known.png',dpi=80)
    return

def getMultipleResults(path):
    """Fetch results for all dirs and aggregate read counts"""

    k = []
    n = []
    m = []
    outdirs = [os.path.join(path,i) for i in os.listdir(path)]
    labels = []
    c=1
    for o in outdirs:
        r = readResultsFile(o, 'mature_sense.grouped')
        x = getNovel(o)
        if x is not None:
            n.append(x)
        if r is not None:
            k.append(r)
        labels.append(c)
        c+=1
        iso = getIsomiRs(o)
        if m is not None:
            m.append(iso)

    k = pd.concat(k)
    if len(n)>0:
        n = pd.concat(n)
    else:
        n=None
    #combine known into useful format
    p = k.pivot_table(index='name', columns='path', values='read count')
    cols = p.columns
    samples = float(len(cols))
    g = k.groupby('name').agg({'read count':[np.size,np.mean,np.sum]})
    g.columns = ['freq','mean read count','total']
    g['perc'] = g['total']/g['total'].sum()
    g['freq'] = g.freq/float(samples)
    k = p.merge(g,left_index=True,right_index=True)
    k = k.reset_index()
    k = k.sort('mean read count',ascending=False)
    #combine isomirs
    if len(m)>0:
        m = pd.concat(m)
        m = m.pivot_table(index=['read','name'], columns='path', values='read count')
        m = m.fillna(0)
        m = m.reset_index()
        m['total'] = m.sum(1)
        m['mean read count'] = m[cols].mean(1)
        m['freq'] = m[cols].apply(lambda r: len(r.nonzero()[0])/samples,1)
        m=m[m.freq>0.3]
        m=m.sort(['total'],ascending=False)
        #print m[['name','total','freq']]
    return k,n,m

def getNovel(path):
    """Parse novel.txt file if available"""

    f = os.path.join(path,'novel.txt')
    if not os.path.exists(f):
        return
    df = pd.read_csv(f, delim_whitespace=True, usecols=range(16), converters={'chrom':str})
    df['path'] = os.path.basename(path)
    return df

def getIsomiRs(path):
    """Get isomiR results"""
    f = os.path.join(path,'miRBase_isoAnnotation.txt')
    if not os.path.exists(f):
        return
    df = pd.read_csv(f, sep='\t')
    df['path'] = os.path.basename(path)
    return df

def summariseAll(k,n,iso,outpath=None):
    """Summarise multiple results"""

    if outpath != None:
        os.chdir(outpath)
    ky1 = 'unique reads'
    ky2 =  'read count' #'RC'
    cols = ['name','freq','mean read count','total']
    print
    print 'found:'
    final = k[(k['mean read count']>=10) & (k['freq']>=.8)]
    print final[cols]
    print '-------------------------------'
    print '%s total' %len(k)
    print '%s with >=10 mean reads' %len(k[k['mean read count']>=10])
    print '%s found in 1 sample only' %len(k[k['freq']==1])
    print 'top 10 account for %2.2f' %k['perc'][:10].sum()

    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    k.set_index('name')['total'][:10].plot(kind='barh',colormap='Spectral',ax=ax,log=True)
    plt.tight_layout()
    fig.savefig('srnabench_top_known.png')
    fig = plotReadCountDists(final)
    fig.savefig('srnabench_known_counts.png')
    print
    if iso is not None:
        analyseIsomiRs(iso)
    return k

def analyseIsomiRs(iso):
    """Analyse isomiR results"""

    plt.close('all')
    print 'top isomiRs:'
    subcols = ['name','read','total','freq']
    print iso[subcols][:10]
    #top isomir per group
    top = iso.sort('total', ascending=False).groupby('name', as_index=False).first()
    top[subcols].to_csv('srnabench_top_isomirs.csv')
    #stats
    g=iso.groupby('name').agg({'total':[np.sum,np.size]})#.sort('sum',ascending=False)
    g.columns = g.columns.get_level_values(1)
    #g[:40].plot(kind='barh')
    fig,ax = plt.subplots(1,1)
    g.plot('size','sum',kind='scatter',logy=True,logx=True,alpha=0.8,s=50,ax=ax)
    ax.set_title('no. isomiRs per miRNA vs total reads')
    fig.savefig('srnabench_isomirs.png')
    #print g[:30]
    return

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

def runDE(inpath, l1, l2, cutoff=1.5):
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
    base.seabornsetup()
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-r", "--run", dest="run", action='store_true',
                           help="run predictions")
    parser.add_option("-i", "--input", dest="input",
                           help="input path or file")
    parser.add_option("-a", "--analyse", dest="analyse",
                           help="analyse results of runs")
    #parser.add_option("-s", "--summarise", dest="summarise",
    #                       help="analyse results of multiple runs together")
    parser.add_option("-c", "--config", dest="config",
                            help="config file")
    opts, remainder = parser.parse_args()
    pd.set_option('display.width', 800)
    if opts.run == True:
        if opts.config == None:
            base.writeDefaultConfig('srnabench.conf',defaults=srnabenchoptions)
            cp = base.parseConfig('srnabench.conf')#opts.config)
            options = cp._sections['base']
            print options
        runAll(opts.input, filetype='fa')
    elif opts.analyse != None:
        k,n,iso = getMultipleResults(opts.analyse)
        summariseAll(k,n,iso)
    else:
        test()

if __name__ == '__main__':
    main()
