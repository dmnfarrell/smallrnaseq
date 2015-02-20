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
import srnabench as srb
import mirdeep2 as mdp

def first(x):
    return x.iloc[0]

def readLengthDist(df):

    df['length'] = df.seq.str.len()
    bins=np.linspace(1,df.length.max(),df.length.max())
    fig,ax=plt.subplots(1,1,figsize=(10,6))
    #df.hist('length',bins=bins,ax=ax,normed=True)
    x=np.histogram(df.length,bins=bins)
    ax.bar(x[1][:-1],x[0], align='center')
    plt.title('read length distribution')
    plt.xlabel('length')
    plt.ylabel('frequency')
    plt.tight_layout()
    plt.savefig('readlengths_dist.png',dpi=150)
    plt.show()
    return

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

def summariseFastq(f):

    ext = os.path.splitext(f)[1]
    if ext=='.fastq':
        ffile = HTSeq.FastqReader(f, "solexa")
    elif ext == '.fa':
        ffile = HTSeq.FastaReader(f)
    else:
        return
    sequences = [(s.name,s.seq) for s in ffile]
    df = pd.DataFrame(sequences,columns=['id','seq'])
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

def trimAdapters(infile, adapters=[], outfile='cut.fastq'):
    """Trim adapters using cutadapt"""

    #if os.path.exists(outfile):
    #    return
    if len(adapters)==0:
        print 'no adapters!'
        return
    adptstr = ' -a '.join(adapters)
    cmd = 'cutadapt -m 18 -O 5 -q 20 --discard-untrimmed -a %s %s -o %s' %(adptstr,infile,outfile)
    print cmd
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    #print result
    return

def removeKnownRNAs(path, adapters=[], outpath='RNAremoved'):
    """Map to annotated RNAs and put remaining reads in output dir"""

    index = 'bosTau6-tRNAs'
    params = '-v 0 --best'
    files = glob.glob(os.path.join(path,'*.fastq'))
    #files = ['test.fastq']
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for f in files:
        print f
        label = os.path.splitext(os.path.basename(f))[0]
        trimAdapters(f, adapters)
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
    print files
    for f in files:
        print f
        label = os.path.splitext(os.path.basename(f))[0]
        cut = os.path.join(outpath,label+'_cut.fastq')
        cfile = os.path.join(outpath,label+'.fa')
        print cut,cfile
        if not os.path.exists(cfile):
            if not os.path.exists(cut) and adapters!=None:
                trimAdapters(f, adapters, cut)
            elif adapters==None:
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

            else: perc = 0.0
            print '%s: %.4f of total reads aligned' %(index,perc)
            x.append(perc)
        res.append(x)

    df = pd.DataFrame(res,columns=colnames)
    df.to_csv('ncrna_mapped.csv',float_format='%.3f')
    return

def plotRNAmapped(labels=None):
    """Plot RNA map results"""

    df = pd.read_csv('ncrna_mapped.csv',index_col=0)
    df=df.sort('total')
    df = df.drop('total',1)
    df['unmapped'] = 1-df.sum(1)
    df['name'] = df['name'].apply(lambda r: ' '.join(r.split('_')[2:5]))
    df = df.set_index('name')
    print df
    if labels != None:
        df = df.rename(columns=labels)
    plt.figure(figsize=(8,8))
    x=df.mean()
    x.sort()
    explode = [0.00 for i in range(len(x))]
    x.plot(kind='pie',colormap='Spectral',autopct='%.1f%%',startangle=0,
                    labels=None,legend=True,pctdistance=1.1,explode=explode,fontsize=16)
    plt.title('mean percentage small RNAs mapped by category')
    plt.tight_layout()
    plt.savefig('ncrna_means.png',dpi=150)

    df=df.reindex(columns=x.index)
    l = df.plot(kind='bar',stacked=True,cmap='Spectral',figsize=(12,6))
    plt.ylabel('percent mapped')
    plt.legend(ncol=4)
    plt.tight_layout()
    plt.savefig('ncrna_bysample.png',dpi=150)

    '''plt.figure(figsize=(12,8))
    ax=base.sns.boxplot(mp,color='Paired')
    plt.ylabel('percent mapped')
    plt.title('mean percentage small RNAs mapped by category')
    plt.tight_layout()
    plt.savefig('ncrna_dists.png')'''
    #plt.show()
    return

def compareMethods():
    """Compare 2 methods for subsets of samples"""

    path1 = 'results_mirdeep_rnafiltered'
    path2 = 'results_srnabench_isomirs'

    #compare means of filtered knowns
    df = mdp.getResults(path1)
    df = df[df.novel==False]
    mk = mdp.filterExprResults(df,meanreads=200,freq=0.8)
    k,n,i = srb.getResults(path2)
    sk = k[(k['mean read count']>=10) & (k['freq']>=0.8)]
    print sk.columns
    #sk = k[(k['total']>=500)]
    x = pd.merge(mk,sk,left_on='#miRNA',right_on='name',how='inner',suffixes=['1','2'])
    diff = x[(abs(x.total1-x.total2)/x.total2>.3)]
    print diff[['#miRNA','total1','total2']]
    #print np.corrcoef(np.log10(x.total1),np.log10(x.total2))
    fig = plt.figure(figsize=(12,6))
    ax=fig.add_subplot(121)
    base.venndiagram([mk['#miRNA'], sk['name']],['mirdeep2','srnabench'],ax=ax)

    ax=fig.add_subplot(122)
    x.plot('total1','total2',kind='scatter',ax=ax, logx=True,logy=True,alpha=0.8,s=40)
    ax.plot([0, 1], [0, 1], transform=ax.transAxes,color='red',alpha=0.7)
    ax.set_xlabel('mirdeep2')
    ax.set_ylabel('srnabench')
    ax.set_title('total read count comparison')
    plt.tight_layout()
    plt.savefig('mirdeep_vs_srnabench.png',dpi=150)
    plt.show()
    return

def KStest(df):
    """KS test to determine rc freq distributions of replicates
    threshold count is the one which corresponds to the first minimum.
    i.e.  when the distributions of reads amongst replicates begin to be similar
    see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2822534/  """

    def getMin(x):
        r=0
        x=list(x)
        m=x[0]
        for i in range(len(x)):
            if x[i]<m and x[i+1]>x[i] and abs(m-x[i])>0.02:
                r=readcuts[i]
                break
        return r

    from scipy.stats import ks_2samp
    readcuts = np.arange(0,400,10)
    result=[]
    ids=[('s11','s22'),('s23','s28'),('s27','s29'),('s26','s03'),('s36','s37'),
        #,('s53','s54'),('s38','s39'),('s41','s42'),
        ('s44','s45'),('s46','s47')]
    testdata=[]
    for i in ids:
        stats=[]
        for rc in readcuts:
            d = mdp.filterExprResults(df,score=-500,freq=0,meanreads=rc,totalreads=0)
            testdata.append(d)
            s1,s2=i
            x=d[s1+'(norm)']
            y=d[s2+'(norm)']
            val = ks_2samp(x, y)
            print len(x),len(y),rc, val
            stats.append(val[0])
        result.append(stats)
    result = pd.DataFrame(np.array(result).T)
    result=result.set_index(readcuts)
    result.plot(lw=2,colormap='Set1',legend=False)
    tmins = result.apply(getMin)
    print tmins
    mean=tmins.mean()
    std=tmins.std()

    plt.axvspan(mean-std/2,mean+std/2,color='g',alpha=0.3)
    plt.xlabel('read count threshold')
    plt.ylabel('KS')
    plt.savefig('KS_test.png')

    '''f,axs=plt.subplots(4,3)
    grid=axs.flat
    i=0
    for d in testdata[:12]:
        ax=grid[i]
        print len(d)
        bins = 10 ** np.linspace(np.log10(10), np.log10(1e6),80)
        ax.hist(d['s22(norm)'],bins=bins,alpha=0.9)
        #ax.hist(d['s11(norm)'],bins=bins,alpha=0.9)
        ax.set_xscale('log')
        ax.set_xticklabels([])
        i+=1'''
    plt.show()
    return

def mirnaDiscoveryTest(sourcefile):
    """Test miRNAs found per read file size"""

    path = 'benchmarking'
    #sourcefile = 'small.fastq'
    #sizes = np.arange(5e5,7.e6,5e5)
    #sizes = np.insert(sizes,0,[5e4,1e5,2e5])
    #base.createRandomFastqFiles(sourcefile, path, sizes)

    #mirdeep2
    #mdp.runMultiple(path, 'fa')
    outpath = 'benchmarking/mirdeep'
    a = mdp.getResults(outpath)
    #a=pd.concat([k,n])
    mapping = mdp.getFilesMapping(outpath)
    results1 = []
    found = []
    for i,r in mapping.iterrows():
        f = r.filename
        i = float(re.findall("\d+.\d+",f)[0])
        df = a[a[r.id]>0]
        new = df[-df['#miRNA'].isin(found)]
        found.extend(new['#miRNA'].values)
        rcm = new['read_count'].mean()
        results1.append((i,len(found),len(new),rcm))
    r1 = pd.DataFrame(results1,columns=['reads','total','new','rcm'])

    #srnabench
    files = glob.glob(os.path.join(path,'*.fa'))
    outpath = 'benchmarking/srnabench'
    results2 = []
    found = []
    for f in sorted(files):
        outdir = srb.run(f,outpath,overwrite=False)
        k,n = srb.getResults(outdir,plot=False)
        i = float(re.findall("\d+.\d+",f)[0])
        new = k[-k['name'].isin(found)]
        found.extend(new.name.values)
        rcm = new['RPM (lib)'].mean()
        results2.append((i,len(found),len(new),rcm))

    r2 = pd.DataFrame(results2,columns=['reads','total','new','rcm'])
    final = pd.merge(r1,r2,on='reads')
    print final
    fig=plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    #final.plot('reads','new_x',kind='scatter',style='.-',color='blue',ax=ax)

    l1=ax.plot(final.reads,final.total_x,color='blue',lw=2,ls='-')
    #ax2 = ax.twinx()
    l2=ax.plot(final.reads,final.total_y,color='green',lw=2,ls='-')
    ax.set_xlabel('reads (million)')
    ax.set_ylabel('total miRNA found')
    #ax2.set_ylabel('mirna found (srnabench)')
    ax.legend(l1+l2,['mirdeep2','srnabench'],loc=2)
    ax3 = plt.axes([.5, .25, .3, .3])
    ax3.semilogy(final.reads,final.rcm_x,color='blue',ls='-')
    ax3.semilogy(final.reads,final.rcm_y,color='green',ls='-')
    ax3.set_title('mean abundance (new miRNAs)')
    plt.tight_layout()
    fig.savefig('benchmark_discovery.png',dpi=150)
    plt.show()
    return

def novelConservation():
    import ensembl
    df = pd.read_csv('novel_mirdeep.csv')
    ensembl.getmiRNAOrthologs(df)
    ensembl.summarise(df)
    return

def getExpCondMap(path):
    """Get mirdeep results labels mapped to our exp condition table"""

    condmap = pd.read_csv('mapdata_labels.csv')
    idmap = mdp.getFilesMapping(path)
    def getIDMap(x):
        r=idmap[idmap['filename'].str.contains(x)]
        if len(r)>0:
            return r.id.squeeze()
        else:
            return np.nan
    condmap['id'] = condmap.filename.apply(lambda x: getIDMap(x),1)
    condmap = condmap.dropna()
    return condmap

def getSubGroups(condmap, column, catnames=None, subcats=None):
    """Get sub groups of conditions based on any factor column"""

    c = condmap
    if subcats != None:
        c = condmap[condmap[column].isin(subcats)]
    c = c.sort(column)
    x = pd.Categorical(c[column],categories=catnames)
    #labels for DE analysis (used in edgeR)
    c['factor'] = x.labels
    c['label'] = c.apply(lambda x: x.id+'_'+str(x.factor),1)
    c = c.sort('factor')
    #print c
    return c

def DEbyTime(df,condmap):
    """DE over time groups"""

    df = df.set_index('#miRNA')
    infection = ['CTRL','MAP']
    elisa = ['N','P']
    #comparedgroups = [0, 6, 43, 46, 49]
    comparedgroups = [['START','EARLY'], ['EARLY','LATE']]
    groupnames = ['START','EARLY','LATE']
    i=0
    res = []
    for inf in elisa:
        for tps in comparedgroups:
            c = condmap[condmap.elisa==inf]
            c = getSubGroups(c, 'timegroup', groupnames, tps)
            cols = c.id
            data = df[cols]
            data.columns = c.label
            data.to_csv('decounts.csv')
            de = base.runEdgeR('decounts.csv', 1.5)
            print de
            clabel = '_'.join(tps)+' '+inf
            i+=1
            de['comp'] = clabel
            res.append(de)

    #get heatmap
    allde = pd.concat(res)
    #plt.show()
    return allde

def DEbyInfection(df,condmap):
    """DE over time groups"""

    df = df.set_index('#miRNA')
    times = ['START','EARLY','LATE']
    i=0
    res = []
    for tg in times:
        c = condmap[condmap.timegroup==tg]
        c = getSubGroups(c, 'elisa', ['N','P'])
        cols = c.id
        data = df[cols]
        data.columns = c.label
        data.to_csv('decounts.csv')
        de = base.runEdgeR('decounts.csv', 1.5)
        print de
        de['comp'] = 'CTRLVsMAP'+'_'+tg
        res.append(de)
        i+=1
    allde = pd.concat(res)
    return allde

def DE():
    path = 'results_mirdeep_rnafiltered'
    df = mdp.getResults(path)
    df = mdp.filterExprResults(df[df.novel==False],score=0,freq=.9,meanreads=200)
    condmap = getExpCondMap(path)
    det = DEbyTime(df,condmap)
    dei = DEbyInfection(df,condmap)
    allde = pd.concat([det,dei])
    print allde
    allde.to_csv('DE_all.csv')
    #DEheatmap(det)
    return

def plotFactors(path):
    """Plots for various views of mirdeep results using seaborn facetgrids"""

    df = mdp.getResults(path)
    df = mdp.filterExprResults(df,meanreads=200)
    df = df[df.novel==False]
    condmap = getExpCondMap(path)
    df=df.set_index('#miRNA')
    cols = [i for i in df.columns if (i.startswith('s') and len(i)<=3)]
    normcols = [i+'(norm)' for i in cols]
    df = df[normcols]
    names = pd.read_csv('DE_all.csv').name
    #names = df.index[:60]
    #names = ['bta-miR-150','bta-miR-151-3p','bta-miR-186','bta-miR-205',
    #        'bta-miR-92b','bta-miR-29a','bta-miR-101','bta-miR-423-5p','bta-miR-486']
    df=df[df.index.isin(names)]
    tporder = [0, 6, 43, 46, 49]
    tgorder = ['START','EARLY','LATE']
    #get numeric values for time points
    x = pd.Categorical(condmap['month'],categories=tporder)
    condmap['time'] = x.codes
    #prepare data by melting it
    t=df.T
    t.index = cols
    t = t.merge(condmap,left_index=True,right_on='id')
    #t=t[t.timegroup=='LATE']
    tm = pd.melt(t,id_vars=list(condmap.columns),
                   var_name='miRNA',value_name='read count')
    g = base.sns.factorplot('time','read count','elisa', tm, col='miRNA', kind="point",
                            col_wrap=2,size=4,aspect=0.9,legend_out=True,sharey=False)
    g.set(xticklabels=tporder)
    plt.savefig('de_lineplots_timepoints.png')
    g = base.sns.factorplot('timegroup','read count','elisa', tm, col='miRNA', kind="bar",x_order=tgorder,
                            col_wrap=4,size=4,aspect=0.9,legend_out=True,sharey=False)
    plt.savefig('de_timegroups.png')

    g = base.sns.factorplot('miRNA','read count', 'elisa', tm, kind="box",
                            size=8,aspect=1.5,legend_out=False,x_order=names)
    g.despine(offset=10, trim=True)
    g.set(yscale='log')
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)
    plt.tight_layout()
    plt.savefig('counts_byinfection.png')

    g = base.sns.factorplot('timegroup','read count','animal', data=tm, col='miRNA', kind="bar",
                            col_wrap=4,size=4,aspect=0.9,legend_out=True,sharey=False,x_order=tgorder,
                            palette='Spectral')
    plt.savefig('counts_byanimal.png')

    x = dict(list(t.groupby('elisa')))
    import scipy
    ts,p = scipy.stats.ttest_ind(x['N'][names], x['P'][names])
    pv = pd.Series(p, index=names).order()
    print pv[pv<0.05]
    #p=t[['pool','infection']]
    #fig,ax = plt.subplots(1,1)
    #p.groupby(['pool'])
    #p.plot(kind='bar',by='pool')

    #plt.show()
    return

def analyseisomirs():
    """check isomir or other expr profiles by category"""

    fmap = pd.read_csv('srnabench_colnames.csv')
    condmap = pd.read_csv('mapdata_labels.csv')
    fmap['filename'] = fmap.filename.str[:-8]
    c = fmap.merge(condmap,on='filename')
    iso = pd.read_csv('srnabench_isomirs_all.csv')
    i=0
    iso = iso[iso.name=='bta-miR-486'][:10]
    #iso = iso[iso.read=='TCCTGTACTGATCTGCCCCGA']
    df = iso.set_index('read')
    t = df.T
    t = t.merge(c,left_index=True,right_on='col')
    print c.columns
    tm = pd.melt(t,id_vars=list(c.columns),
                   var_name='read',value_name='total')
    print tm
    g = base.sns.factorplot('animal','total', data=tm, col='read', kind="bar", hue='pool',
                            col_wrap=4,size=4,aspect=0.9,legend_out=True,sharey=False)
    plt.show()
    return

def test():
    base.seabornsetup()
    #path = '/opt/mirnaseq/data/vegh_13'
    path = '/opt/mirnaseq/data/combined'
    files = ['/opt/mirnaseq/analysis/test.fastq']
    #files = ['/opt/mirnaseq/data/vegh_13/SRR576286.fastq']
    #bidx =  ['mirbase-mature','Rfam_btau','bosTau6-tRNAs','noncodev4_btau','bos_taurus_alt']
    bidx =  ['mirdeep_found','Rfam_btau','bosTau6-tRNAs','bostau-snRNA','noncodev4_btau','bos_taurus_alt']
    #adapters for our data
    adapters = ['TGGAATTCTCGGGTGCCAAGG','GCATTGTGGTTCAGTGGTAGAATTCTCGC']
    #adapters for Vegh data
    #adapters = ['TAGCTTATCAGACTGATGTTGA','AGATCGGAAGAGCACACGTCTGAACTCC']
    labels = {'bosTau6-tRNAs':'tRNA (GtRNAdb)', 'Rfam_btau':'rRNA (RFAM)',
                'noncodev4_btau':'NONCODE v4', 'bos_taurus_alt':'UMD3.1',
                'mirdeep_found':'miRNA (miRDeep2)'}
    #mapRNAs(path=path, indexes=bidx, adapters=adapters)
    #plotRNAmapped(labels)
    #summariseReads(path)
    #removeKnownRNAs(path, adapters)
    #compareMethods()
    infile = '/opt/mirnaseq/data/combined/miRNA_lib_Pool2_Sample_2_combined.fastq'
    #mirnaDiscoveryTest(infile)
    #novelConservation()
    #DE()
    plotFactors('results_mirdeep_rnafiltered')
    #compareIsomirsRef()
    #analyseisomirs()
    #checknovel()
    return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--summarise", dest="summarise",
                           help="summarise fastq reads")
    parser.add_option("-i", "--input", dest="input",
                           help="input path or file")
    parser.add_option("-t", "--test", dest="test", action='store_true',
                           help="testing")
    opts, remainder = parser.parse_args()
    pd.set_option('display.width', 800)
    base.seabornsetup()
    if opts.summarise != None:
        if os.path.isdir(opts.summarise):
            summariseReads(opts.summarise)
        else:
            df = summariseFastq(opts.summarise)
            print '%s reads' %len(df)
            readLengthDist(df)
    elif opts.test == True:
        test()

if __name__ == '__main__':
    main()
