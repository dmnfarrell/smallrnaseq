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
                    ('bowtieindex',''),('ref',''),  ('predict','false'),
                    ('mature',''), ('hairpin',''), ('other',''), ('isomir','false'),
                    ('overwrite',1), ('matureMM',0), ('p',3)]}
isoclasses = {'lv5pT':'5p trimmed',
                'lv5pE':'5p extended',
                'lv5p':'5p length variant',
                'lv3pT':'3p trimmed',
                'lv3pE':'3p extended',
                'lv3p':'3p length variant',
                'mv': 'multiple length variants'}

def getShortlabel(label):
    x=label.split('_')
    return x[2]+'_'+x[4]

def run(infile, outpath='srnabench_runs', overwrite=True, adapter='',
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
           ' species=%s output=%s predict=%s plotMiR=true matureMM=0 isoMiR=true' #hierarchical=false'
           ' p=3' %(srbpath,srbpath,infile,ref,outdir,predict))
    if adapter != '':
        cmd += ' adapter=%s' %adapter
    print cmd
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print result
    return outdir

def runAll(path, outpath='runs', filetype='fastq', **kwargs):
    """Run all fastq files in folder"""

    if os.path.isdir(path):
        files = glob.glob(os.path.join(path,'*.'+filetype))
    else:
        files = [path]
    print 'running sRNAbench for %s files' %len(files)
    for f in files:
        res = run(f, outpath, **kwargs)
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

def normaliseCols(df, cols):
    def n(s):
        return s/s.sum()
    print df[cols].apply(n, axis=0)
    return

def getResults(path):
    """Fetch results for all dirs and aggregate read counts"""

    k = []
    n = []
    m = []
    outdirs = [os.path.join(path,i) for i in os.listdir(path)]
    cols = []
    fmap = {} #filename mapping
    c=1
    for o in outdirs:
        r = readResultsFile(o, 'mature_sense.grouped')
        x = getNovel(o)
        if x is not None:
            #print len(x), o
            n.append(x)
        if r is not None:
            k.append(r)
        id='s'+str(c)
        cols.append(id)
        fmap[id] = os.path.basename(o)
        c+=1
        iso = getIsomiRs(o)
        if iso is not None:
            m.append(iso)
    fmap = pd.DataFrame(fmap.items(),columns=['col','filename'])
    fmap.to_csv('srnabench_colnames.csv',index=False)
    #combine known into useful format
    k = pd.concat(k)
    p = k.pivot_table(index='name', columns='path', values='read count')
    #cols = p.columns
    p.columns = cols
    samples = float(len(cols))
    g = k.groupby('name').agg({'read count':[np.size,np.mean,np.sum]})
    g.columns = ['freq','mean read count','total']
    g['perc'] = g['total']/g['total'].sum()
    g['freq'] = g.freq/float(samples)
    k = p.merge(g,left_index=True,right_index=True)
    k = k.reset_index()
    k = k.fillna(0)
    k = k.sort('mean read count',ascending=False)
    #normaliseCols(k, cols)
    #combine isomir results
    if len(m)>0:
        m = pd.concat(m)
        m = m.pivot_table(index=['read','name','isoClass','NucVar'],
                            columns='path', values='read count')
        m.columns = cols
        m = m.fillna(0)
        m = m.reset_index()
        m['total'] = m.sum(1)
        m['mean read count'] = m[cols].mean(1)
        m['freq'] = m[cols].apply(lambda r: len(r.nonzero()[0])/samples,1)
        m=m.sort(['total'],ascending=False)
    else: m = None
    #novel
    if len(n)>0:
        n = pd.concat(n)
        print n[n.columns[:8]].sort(['chrom','5pSeq'])
        n = n.pivot_table(index=['5pSeq','chrom'], columns='path', values='5pRC')
        print n[n.columns[:3]]
        #n.columns=cols
    else:
        n=None
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

def analyseResults(k,n,outpath=None):
    """Summarise multiple results"""

    if outpath != None:
        os.chdir(outpath)
    ky1 = 'unique reads'
    ky2 =  'read count' #'RC'
    cols = ['name','freq','mean read count','total','perc']
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
    fig,ax = plt.subplots(figsize=(10,6))
    k[k.columns-cols].sum().plot(kind='bar',ax=ax)
    fig.savefig('srnabench_total_persample.png')
    print
    return k

def analyseIsomiRs(iso,outpath=None):
    """Analyse isomiR results in detail"""

    if iso is None:
        return
    if outpath != None:
        os.chdir(outpath)
    subcols = ['name','read','isoClass','NucVar','total','freq']
    iso = iso.sort('total', ascending=False)
    #filter low abundance reads same as profiling
    iso = iso[(iso.total>10) & (iso.freq>0.6)]
    iso['length'] = iso.read.str.len()
    #parse classes info
    def parseisoinfo(r):
        '''parse srnabench hierarchical scheme'''
        s = r['isoClass'].split('|')
        lv = s[0]
        if len(s)>1:
            pos = int(s[-1].split('#')[-1])
        else:
            pos = 0
        return pd.Series(dict(pos=pos,variant=lv))
    x = iso.apply(parseisoinfo,1)
    iso = iso.merge(x,left_index=True,right_index=True)

    #get top isomir per mirRNA
    g = iso.groupby('name', as_index=False)
    t=[]
    for i,x in g:
        r = base.first(x)
        s = x.total.sum()
        perc = r.total/s
        t.append((r['name'],r.read,r.total,s,perc,np.size(x.total),r.variant))
    top = pd.DataFrame(t,columns=['name','read','counts','total','domisoperc','isomirs','variant'])
    top = top.sort('counts',ascending=False)
    top.to_csv('srnabench_isomirs_dominant.csv',index=False)
    print 'top isomiRs:'
    print top[:20]
    print '%s/%s with only 1 isomir' %(len(top[top.domisoperc==1]),len(top))
    print 'different dominant isomir:', len(top[top.variant!='exact'])/float(len(top))
    print top.domisoperc.mean()
    print
    #stats
    fig,ax = plt.subplots(1,1)
    top.plot('isomirs','total',kind='scatter',logy=True,logx=True,alpha=0.8,s=50,ax=ax)
    ax.set_title('no. isomiRs per miRNA vs total adundance')
    ax.set_xlabel('no. isomiRs')
    ax.set_ylabel('total reads')
    fig.savefig('srnabench_isomirs.png',dpi=150)
    #length dists of isomirs
    fig,ax = plt.subplots(1,1)
    x = iso[iso.name.isin(iso.name[:80])]
    bins=range(15,30,1)
    x.hist('length',bins=bins,ax=ax,by='name',sharex=True,color="gray")
    #base.sns.distplot(iso.length, kde=False)
    fig.suptitle('isomiR length distributions')
    fig.savefig('srnabench_isomir_lengths.png',dpi=150)
    #l = g.agg({'length':[np.median,np.std]})

    #print iso[subcols]
    plt.close('all')

    #diff = ['bta-miR-27a-3p','bta-miR-127','bta-miR-101','bta-miR-23b-3p']
    #d = iso[iso.name.isin(diff)][subcols+['pos','variant']]
    #for i,df in d.groupby('name'):
    #    print df
    #print iso[iso.pos<-5][subcols]

    #lv = iso[-iso.variant.str.contains('exact')]
    c=iso.variant.value_counts()
    #c=c[c>10]
    fig,ax = plt.subplots(1,1)
    c.plot(kind='pie',colormap='Spectral',ax=ax)
    '''c.plot(kind='bar',ax=ax)
    ax.set_ylabel('total')
    ax.set_xticklabels(c.index, minor=False, rotation=45)'''
    ax.set_title('isomiR class distribution')
    plt.tight_layout()
    fig.savefig('srnabench_isomir_classes.png',dpi=150)
    fig,axs = plt.subplots(2,2)
    grid=axs.flat
    bins=np.arange(-6,6,1)
    i=0
    for v in ['lv3p','lv5p','nta#A','nta#T']:
        ax=grid[i]
        iso[iso.variant==v].hist('pos',ax=ax,bins=bins)
        ax.set_title(v)
        ax.set_xticks(bins+0.5)
        ax.set_xticklabels(bins)
        i+=1
    fig.suptitle('isomiR variant lengths')
    fig.savefig('srnabench_isomir_variantlengths.png')
    #plt.show()
    iso.to_csv('srnabench_isomirs_all.csv',index=False)
    return top

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

    parser.add_option("-c", "--config", dest="config",
                            help="config file")
    opts, remainder = parser.parse_args()
    pd.set_option('display.width', 600)
    if opts.run == True:
        if opts.config == None:
            base.writeDefaultConfig('srnabench.conf',defaults=srnabenchoptions)
        cp = base.parseConfig(opts.config)
        options = cp._sections['base']
        print options
        runAll(opts.input, **options)
    elif opts.analyse != None:
        k,n,iso = getResults(opts.analyse)
        analyseResults(k,n)
        analyseIsomiRs(iso)
    else:
        test()

if __name__ == '__main__':
    main()
