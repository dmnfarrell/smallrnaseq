#!/usr/bin/env python

"""Module for miRDeep2 wrappers and utilities
   Created July 2014
   Copyright (C) Damien Farrell
"""

import sys, os, string, types, re, StringIO
import shutil, glob, collections
import itertools
import subprocess
import pylab as plt
import numpy as np
import pandas as pd
import HTSeq
import base

mirdeep2options = {'base': [('input',''),('adapter','TGGAATTCTCGGGTGCCAAGG'),('filetype','fastq'),
                    ('bowtieindex',''),('refgenome',''),('species','hsa'),
                    ('mature',''), ('hairpin',''), ('other',''),('mirbase',os.getcwd()),
                    ('overwrite',1)]}
mirdeepcols = ['#miRNA','read_count','mean_norm','miRDeep2 score','chr','seed','precursor',
                'precursor coordinate','freq','mirbase seed match','star read count','rfam alert',
                'consensus mature sequence','consensus star sequence',
                'consensus precursor sequence']

def createMiRBaseFiles(species,path):
    """Generate species specific mature/hairpin files for input to mirdeep"""

    mature = os.path.join(path, 'mature.fa')
    hairpin = os.path.join(path, 'hairpin.fa')
    names=[]
    for f in [mature,hairpin]:
        fname = os.path.splitext(f)[0]+'_'+species+'.fa'
        base.getSubsetFasta(f, labels=[species], outfile=fname)
        names.append(fname)
    print 'wrote mirbase files for species %s' %species
    return names

def createSampleMap(path, ext='fastq'):
    """Create filename mapping to run all samples at once.
       This is required."""

    os.chdir(path)
    files = sorted(glob.glob('*.'+ext))
    print files
    fname = 'combined.txt'
    i=1
    rows=[]
    for f in files:
        rows.append((os.path.basename(f), 's%02d'%i))
        i+=1
    res=pd.DataFrame(rows)
    res.to_csv(fname, index=False,sep=' ',header=False)
    return fname

def runMultiple(**kwargs):
    """Prepare and run mirdeep2"""

    if kwargs['filetype'] == 'fasta':
        ext='fa'
    else:
        ext='fastq'
    path = kwargs['input']
    #create filename/id mapping
    samplemap = createSampleMap(path, ext)
    #get mirbase subset for species if provided
    if kwargs['species'] != '':
        mature, hairpin = createMiRBaseFiles(kwargs['species'], kwargs['mirbase'])
        kwargs['mature'] = mature
        kwargs['hairpin'] = hairpin
    run(samplemap, **kwargs)
    return

def run(infile, refgenome, bowtieindex, mature='', hairpin='', other='',
        randfold=True, overwrite=False, filetype='fastq', adapter=None,
        clean=True, outpath=None, **kwargs):
    """Run all mirdeep2 steps including adapter trimming.
       Uses a config file even if we only have one sample."""

    label = os.path.splitext(os.path.basename(infile))[0]
    print 'running %s' %label
    os.environ["BOWTIE_INDEXES"] = os.path.dirname(bowtieindex)
    collapsed = 'collapsedreads.fa'
    if filetype=='fasta': params='-c'
    else: params='-e -h'
    if randfold == False: params+=' -c'
    other = 'none'

    #if mapping has been done already we can skip it
    if not os.path.exists('mapped.arf') or overwrite == True:
        try:
            os.remove('mapped.arf')
            os.remove(collapsed)
        except:
            pass
        cmd1 = ('mapper.pl %s -d %s -j -l 18 -m -k %s -s %s'
                ' -p %s -t mapped.arf -v' %(infile,params,adapter,collapsed,bowtieindex))
        print cmd1
        result = subprocess.check_output(cmd1, shell=True, executable='/bin/bash')
    else:
        print 'arf file found, skipping mapper step'

    #mirdeep core
    cmd2 = ('miRDeep2.pl %s %s mapped.arf'
           ' %s %s %s -z _%s'
           ' -d > report.log' %(collapsed,refgenome,mature,other,hairpin,label))
    print cmd2
    result = subprocess.check_output(cmd2, shell=True, executable='/bin/bash')
    #remove junk
    if clean == True:
        tempfiles = glob.glob('dir_*')
        for f in tempfiles:
            shutil.rmtree(f)
    #move results to dest folder
    if outpath != None:
        pass #move files..
    return

def quantifier(path, mature, precursor, star=None, collapsed='collapsedreads.fa'):
    """Run quantifier module using custom known mature/precursors"""

    current = os.getcwd()
    os.chdir(path)
    cmd = 'quantifier.pl -p %s -m %s -r %s -y novel -k -d' %(precursor,mature,collapsed)
    print cmd
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    os.chdir(current)
    return

def moveResults(dest):
    """Move results from mirdeep to a clean folder"""

    return

def getChromosome(x):
    val = x.split('_')[0]
    try:
        return '%02d' %int(val)
    except:
        return val

def readResultsFile(infile):
    """Get mirdeep results from the summary file"""

    if os.path.splitext(infile)[1] != '.csv':
        return
    scol = 'miRDeep2 score'
    df = pd.read_csv(infile, sep='\t',header=23)
    df = df.dropna()
    #remove junk
    idx = df[df['provisional id']=='tag id'].index[0]
    df = df.drop(idx)
    df['novel'] = np.where(df['miRBase miRNA']=='-',True,False)
    colstorename = {'example miRBase miRNA with the same seed':'mirbase seed match',
         'significant randfold p-value': 'randfold'}
    df = df.rename(columns=colstorename)
    df = df.convert_objects(convert_numeric=True)
    df['chr'] = df['provisional id'].apply(getChromosome)
    df['seed'] = df['consensus mature sequence'].apply(lambda x: x[1:8])
    return df

def getResults(path):
    """Process known and novel results from mirdeep run.
       Combines with expression data and removes redundant entries"""

    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    df = readResultsFile(resfile)

    #use quantifier module to get novel expression results from predicted precursors if not done already
    novelmature = os.path.join(path, 'novel_mature.fa')
    novelprecursor = os.path.join(path, 'novel_precursor.fa')
    if not os.path.exists(novelmature):
        novel = df[df.novel==True]
        mkey = 'consensus mature sequence'
        pkey = 'consensus precursor sequence'
        base.dataframe2Fasta(novel, mkey, 'provisional id', outfile=novelmature)
        base.dataframe2Fasta(novel, pkey, 'provisional id', outfile=novelprecursor)
        quantifier(path, os.path.abspath(novelmature), os.path.abspath(novelprecursor),
                    'collapsedreads.fa')

    #get expression results and filter by score by merging with predictions
    files = glob.glob(os.path.join(path,'miRNAs_expressed_all_samples*.csv'))
    res=[]

    for f in files:
        if 'novel' in f:
            key='provisional id'
        else:
            key='miRBase miRNA'
        q = pd.read_csv(f,sep='\t')
        samples = float(len(q.filter(regex="norm").columns))
        print 'samples: %s' %samples
        q['freq'] = q.filter(regex="norm").apply(lambda r: len(r.nonzero()[0])/samples,1)
        #apply 5p id so we can merge with results file and keep star seqs
        q['id'] = q['#miRNA'].apply(lambda x: x[:-2]+'5p' if str(x).endswith('3p') else x)
        q = q.merge(df,left_on='id',right_on=key).drop_duplicates('#miRNA')
        res.append(q)
    res = pd.concat(res)

    #get mean normalised count
    res['mean_norm'] = res.filter(regex="norm").apply(lambda r: r[r.nonzero()[0]].mean(),1)
    res = res.sort(['read_count'],ascending=False)
    #res['std'] = res.filter(regex="norm").std(1)
    #res['cv'] = res['std']/res['mean_norm']
    return res

def filterExprResults(n, cols=None, score=0, freq=0.5, meanreads=0, totalreads=50):
    """Additional filters for abundances/no. samples"""

    if cols is None:
        cols = [i for i in n.columns if (i.startswith('s') and len(i)<=3)]
    normcols = [i+'(norm)' for i in cols]
    n = n[(n['miRDeep2 score']>score)]# | (n['read_count']>10000)]
    n = n[n.freq>freq]
    n = n[n['read_count']>totalreads]
    n = n[n['mean_norm']>meanreads]
    n = n[n['randfold'].isin(['yes','-'])]
    #n = n[n['rfam alert']=='-']
    n = n.reset_index(drop=True)
    return n

def analyseResults(path, outpath=None, **kwargs):
    """General analysis of mirdeep results"""

    if outpath != None:
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        os.chdir(outpath)
    df = getResults(path)
    idcols,normcols = getColumnNames(df)
    known = df[df.novel==False]
    novel = df[df.novel==True]
    idmap = getFilesMapping(path)
    k = filterExprResults(known,score=0,freq=.7,meanreads=200)
    n = filterExprResults(novel,score=4,freq=.7,meanreads=200)
    cols = mirdeepcols
    core = pd.concat([k,n])
    base.dataframe2Fasta(core, 'consensus mature sequence', '#miRNA', 'mirdeep_core.fa')

    k[cols].to_csv('known_mirdeep.csv')
    n[cols].to_csv('novel_mirdeep.csv')
    base.createHtml(n[cols],'novel_mirdeep')
    k['perc'] = k['read_count']/k['read_count'].sum()

    print k[cols[:8]]
    print
    print n[cols[:9]]

    print 'mirdeep summary'
    print '-------------------------------'
    print '%s (%s novel) identified' %(len(df),len(novel))
    print 'quantifier results after filtering:'
    print '%s/%s known' %(len(k),len(known))
    print '%s/%s novel' %(len(n),len(novel))
    print 'top 10 known account for %2.2f' %k['perc'][:10].sum()
    print

    k=k.set_index('#miRNA')
    n=n.set_index('#miRNA')

    fig, ax = plt.subplots(figsize=(8,6))
    #k['perc'][:10].plot(kind='pie',colormap='Set2',autopct='%.2f',
    #                 startangle=90,labels=None,legend=True)
    k['read_count'][:10].plot(kind='barh',colormap='Set2',ax=ax,log=True)
    #ax.set_xscale('log')
    plt.title('miRDeep2 top 10')
    plt.tight_layout()
    fig.savefig('mirdeep_top_known.png',dpi=150)
    fig, ax = plt.subplots(figsize=(8,8))
    df.plot('freq','mean_norm',kind='scatter',ax=ax,logy=True,alpha=0.8)
    fig.savefig('mirdeep_freqsvcounts.png',dpi=150)
    fig=plt.figure()

    fig = plotReadCountDists(n,h=5)
    fig.savefig('mirdeep_novel_counts.png',dpi=150)
    fig = plotReadCountDists(k)
    fig.savefig('mirdeep_known_counts.png')

    '''fig=plotSampleDist(n, normcols)
    fig.savefig('mirdeep_novel_persamplecounts.png')
    fig=plotSampleDist(k, normcols)
    fig.savefig('mirdeep_known_persamplecounts.png')'''
    #perSampleDists(k)

    fig,ax = plt.subplots(figsize=(8,6))
    core[idcols].sum().plot(kind='bar',ax=ax)
    #core.columns = idmap.filename
    plt.title('total miRNA counts per sample (unnormalised)')
    plt.tight_layout()
    fig.savefig('mirdeep_total_persample.png')

    #found per chromosome
    fig = plt.figure(figsize=(8,6))
    x = core.sort('chr').groupby('chr').size()
    x.plot(kind='bar')
    plt.title('miRNA per chromosome')
    fig.savefig('mirdeep_chromosome_dist.png')
    #plt.show()
    #plt.close()
    return df,k,n

def plotReadCountDists(df,h=8):
    """Boxplots of read count distributions per miRNA - use seaborn?"""

    w=int(h*(len(df)/60.0))+4
    fig, ax = plt.subplots(figsize=(w,h))
    cols,normcols = getColumnNames(df)
    df = df[normcols]
    t=df.T
    t.index = cols
    #base.sns.boxplot(t,linewidth=1.0,color='coolwarm_r',saturation=0.2,)
    t.plot(kind='box',color='black',grid=False,whis=1.0,ax=ax)
    #base.sns.despine(trim=True)
    ax.set_yscale('log')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylabel('read count')
    plt.tight_layout()
    return fig

def getColumnNames(df):
    cols = [i for i in df.columns if (i.startswith('s') and len(i)<=3)]
    normcols = [i+'(norm)' for i in cols]
    return cols, normcols

def perSampleDists(df,cols=None,names=None):

    base.seabornsetup()
    #df=df[df.index.isin(names)]
    df = df[:2]
    cols,normcols = getColumnNames(df)
    df = df[normcols]
    df= df.reset_index()
    m = pd.melt(df,id_vars=['#miRNA'],
                   var_name='sample',value_name='read count')
    print m
    g = base.sns.factorplot('sample','read count', None, m, col='#miRNA', kind="bar",
                            col_wrap=3,size=3,aspect=1.5,legend_out=True,sharey=False)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)
    plt.tight_layout()
    plt.savefig('mirdeep_persample_counts.png')
    return

def plotSampleDist(df, cols=None, plots=12):
    """Plot distr. of read counts for a mirna per sample"""

    if len(df)>plots:
        df = df[:plots+1]
    c=3
    nr = int(np.ceil(len(df) / c))+1
    fig,axs=plt.subplots(nr,c,sharex=True,figsize=(12,8))
    grid=axs.flat
    i=0
    for n,r in df.iterrows():
        ax=grid[i]
        #name=r['#miRNA']
        if cols is None:
            r = r.filter(regex="norm")
        else:
            r = r[cols]
        pd.Series.plot(r,kind='bar',ax=ax,grid=False)
        ax.set_title(n)
        ax.set_xticks([])
        ax.set_xlabel('')
        i+=1
    plt.tight_layout()
    return fig

def getFilesMapping(path):
    """Get file<->mirdeep2 id mapping"""

    idmap = pd.read_csv(os.path.join(path, 'combined.txt'),sep=' ',
                header=None,names=['filename','id'])
    return idmap

def testQuantifier(path):
    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    df = readResultsFile(resfile)
    novelmature = os.path.join(path, 'novel_mature.fa')
    novelstar = os.path.join(path, 'novel_star.fa')
    novelprecursor = os.path.join(path, 'novel_precursor.fa')
    reads = '../results_mirdeep_combined/collapsedreads.fa'
    #base.dataframe2Fasta(df[df.novel==True], 'consensus star sequence', 'provisional id',
    #            outfile=novelstar)
    quantifier(path, os.path.abspath(novelmature), os.path.abspath(novelprecursor),
                    os.path.abspath(novelstar), reads)
    return

def checkQuantifierResults(path):
    """Check quantifier vs results file"""

    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    df = readResultsFile(resfile)
    files = glob.glob(os.path.join(path,'miRNAs_expressed_all_samples*.csv'))
    q = pd.read_csv(files[0],sep='\t')
    key='provisional id'
    m=q.merge(df,left_on='#miRNA',right_on=key).drop_duplicates('#miRNA')
    m.sc = m['miRDeep2 score']
    m['diff'] = m['read_count']-m['total read count']
    cols=['#miRNA','total read count','read_count','miRDeep2 score']
    print m.sort('diff',ascending=False)[cols]
    #m['size'] = np.select([m.sc < 2, m.sc < 3, m.sc < 4], [20,40,50], 80)
    m.plot(x='total read count',y='read_count', kind='scatter',s=60,alpha=0.6)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    return

def compareRuns(path1,path2):
    """Compare 2 mirdeep runs"""
    df = getResults(path1)
    r1 = filterExprResults(df,meanreads=200)
    df = getResults(path2)
    r2 = filterExprResults(df,meanreads=200)
    x=pd.merge(r1,r2,on='#miRNA',suffixes=['1','2'])
    fig,ax=plt.subplots(1,1)
    x.plot('read_count1','read_count2',kind='scatter',logx=True,logy=True,
            alpha=0.8,s=40,ax=ax)
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", color='red')
    plt.show()
    return

def test(path):
    print path
    checkQuantifierResults(path)
    #testQuantifier(path)
    #df = getResults(path)
    #df = filterExprResults(df,score=0,freq=30,meanreads=0)
    return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-r", "--run", dest="run", action='store_true',
                           help="run predictions")
    parser.add_option("-i", "--input", dest="input",
                           help="input path or file")
    parser.add_option("-c", "--config", dest="config",
                            help="config file")
    parser.add_option("-a", "--analyse", dest="analyse",
                           help="analyse results of mirdeep2")
    parser.add_option("-t", "--test", dest="test", action='store_true',
                           help="testing")
    opts, remainder = parser.parse_args()
    pd.set_option('display.width', 800)
    base.seabornsetup()

    if opts.run == True:
        #all other options are stored in config file
        if opts.config == None:
            print 'No config file provided.'
            base.writeDefaultConfig('mirdeep2.conf',defaults=mirdeep2options)
            return
        cp = base.parseConfig(opts.config)
        if opts.input != None:
            conf.input = os.path.abspath(opts.input)
        options = cp._sections['base']
        runMultiple(**options)
    elif opts.analyse != None:
        analyseResults(opts.analyse)
    elif opts.test == True:
        test(opts.input)

if __name__ == '__main__':
    main()
