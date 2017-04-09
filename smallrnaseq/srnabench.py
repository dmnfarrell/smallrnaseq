#!/usr/bin/env python

"""
    Module for sRNAbench wrappers and utilities
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

import sys, os, string, types, re
import shutil, glob, collections
import itertools
import subprocess
import pylab as plt
import numpy as np
import pandas as pd
from . import base, config

modulepath = os.path.dirname(__file__)
path = os.path.dirname(os.path.abspath(__file__)) #path to module
datadir = os.path.join(path, 'data')
mirbase = pd.read_csv(os.path.join(datadir,'miRBase_all.csv'))
mirbase = mirbase.drop_duplicates('mature1')

srbpath = '/local/sRNAbench'
srnabenchcols = ['name','mean read count','mean_norm','total','freq','perc']
srnabenchoptions = {'base': [('input',''),('outpath','srnabench_runs'),
                    ('adapter','TGGAATTCTCGGGTGCCAAGG'),('filetype','fastq'),
                    ('bowtieindex',''),('ref',''), ('species','bta'), ('predict','false'),
                    ('mature',''), ('hairpin',''), ('other',''), ('isomir','false'),
                    ('overwrite',1), ('matureMM',0), ('p',3)]}
isoclasses = {'lv5pT':'5p trimmed',
                'lv5pE':'5p extended',
                'lv5p':'5p length variant',
                'lv3pT':'3p trimmed',
                'lv3pE':'3p extended',
                'lv3p':'3p length variant',
                'mv': 'multiple length variants'}

def get_short_label(label):
    x=label.split('_')
    return x[2]+'_'+x[4]

def run(infile, outpath='srnabench_runs', overwrite=True, adapter='', species='bta',
            ref='bos_taurus_alt', predict='false', p=2, **kwargs):
    """Run sRNAbench for a fastq file"""

    label = os.path.splitext(os.path.basename(infile))[0]
    resfile = os.path.join(outpath, 'mature_sense.grouped')
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outdir = os.path.join(outpath, label)
    print ('running %s' %infile)

    if os.path.exists(outdir):
        if overwrite == False:
            return outdir
        else:
            shutil.rmtree(outdir)

    cmd = ('java -jar %s/sRNAbench.jar dbPath=%s input=%s microRNA=%s'
           ' species=%s output=%s predict=%s plotMiR=true matureMM=0 isoMiR=true' #hierarchical=false'
           ' p=%s' %(srbpath,srbpath,infile,species,ref,outdir,predict,p))
    if adapter != '':
        cmd += ' adapter=%s' %adapter
    #else:
    #    cmd += ' guessAdapter=true'
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print (result)
    return outdir

def run_all(path, outpath='srnabench_runs', filetype='fastq', **kwargs):
    """Run all fastq files in folder"""

    if os.path.isdir(path):
        files = glob.glob(os.path.join(path,'*.'+filetype))
    else:
        files = [path]
    print ('running sRNAbench for %s files' %len(files))
    for f in files:
        res = run(f, outpath, **kwargs)
        if res == None:
            print ('skipped %s' %f)
    return

def read_results_file(path, infile='mature_sense.grouped', filter=True):
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

def plot_results(k):
    fig,ax=plt.subplots(figsize=(8,6))
    ax.set_title('sRNAbench top 10')
    k.set_index('name')['read count'][:10].plot(kind='barh',colormap='Set2',ax=ax,log=True)
    plt.tight_layout()
    fig.savefig('srnabench_summary_known.png')
    return

def normalise_cols(df):

    def n(s):
        return s/s.sum()*1e6
    cols,normcols = get_column_names(df)
    x = df[cols].apply(n, axis=0)
    x = np.around(x,0)
    x.columns = [i+'_norm' for i in cols]
    df = df.merge(x,left_index=True, right_index=True)
    df['mean_norm'] = df[normcols].apply(lambda r: r[r.nonzero()[0]].mean(),1)
    return df

def get_column_names(df):
    cols = [i for i in df.columns if (i.startswith('s') and len(i)<=3)]
    normcols = [i+'_norm' for i in cols]
    return cols, normcols

def filter_expr_results(df, freq=0.5, meanreads=0, totalreads=50):

    c,normcols = get_column_names(df)
    df = df[df.freq>freq]
    df = df[df['total']>=totalreads]
    df = df[df['mean_norm']>=meanreads]
    return df

def parseisoinfo(r):
    '''parse srnabench hierarchical scheme'''

    s = r['isoClass'].split('|')
    lv = s[0]
    if len(s)>1:
        pos = int(s[-1].split('#')[-1])
    else:
        pos = 0
    return pd.Series(dict(pos=pos,variant=lv))

def get_results(path):
    """Fetch results for all dirs and aggregate read counts"""

    k = []
    n = []
    m = []
    outdirs = [os.path.join(path,i) for i in os.listdir(path)]
    outdirs.sort()
    cols = []
    fmap = {} #filename mapping
    c=1
    for o in outdirs:
        if not os.path.isdir(o):
            continue
        r = read_results_file(o, 'mature_sense.grouped')
        x = get_novel(o)
        if x is not None:
            #print len(x), o
            n.append(x)
        if r is not None:
            k.append(r)
        id='s'+str(c)
        cols.append(id)
        fmap[id] = os.path.basename(o)
        c+=1
        iso = get_isomirs(o)
        if iso is not None:
            m.append(iso)

    fmap = pd.DataFrame(fmap.items(),columns=['id','filename'])
    fmap.to_csv(os.path.join(path,'srnabench_colnames.csv'),index=False)

    #combine known mirnas into one file
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
    k = normalise_cols(k)
    k = k.sort('mean_norm',ascending=False)
    k.to_csv('srnabench_known_all.csv')

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
        m['length'] = m.read.str.len()
        m = normalise_cols(m)
        x = m.apply(parseisoinfo,1)
        m = m.merge(x,left_index=True,right_index=True)
        m = m.sort(['total'],ascending=False)
    else:
        m = None
    #novel
    if len(n)>0:
        n = pd.concat(n)
        #print n[n.columns[:8]].sort(['chrom','5pSeq'])
        n = n.pivot_table(index=['5pSeq','chrom'], columns='path', values='5pRC')
        #print n[n.columns[:3]]
    else:
        n = None
    return k,n,m

def get_novel(path):
    """Parse novel.txt file if available"""

    f = os.path.join(path,'novel.txt')
    if not os.path.exists(f):
        return
    df = pd.read_csv(f, delim_whitespace=True, usecols=range(16), converters={'chrom':str})
    df['path'] = os.path.basename(path)
    return df

def get_isomirs(path):
    """Get isomiR results"""

    f = os.path.join(path,'miRBase_isoAnnotation.txt')
    if not os.path.exists(f):
        return
    df = pd.read_csv(f, sep='\t')
    df['path'] = os.path.basename(path)
    return df

def analyse_results(k,n,outpath=None):
    """Summarise multiple results"""

    if outpath != None:
        os.chdir(outpath)
    #add mirbase info
    k = k.merge(mirbase,left_on='name',right_on='mature1')
    ky1 = 'unique reads'
    ky2 = 'read count' #'RC'
    cols = ['name','freq','mean read count','mean_norm','total','perc','mirbase_id']
    print
    print ('found:')
    idcols,normcols = get_column_names(k)
    final = filter_expr_results(k,freq=.8,meanreads=200)
    print (final[cols])
    print ('-------------------------------')
    print ('%s total' %len(k))
    print ('%s with >=10 mean reads' %len(k[k['mean read count']>=10]))
    print ('%s found in 1 sample only' %len(k[k['freq']==1]))
    print ('top 10 account for %2.2f' %k['perc'][:10].sum())

    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    k.set_index('name')['total'][:10].plot(kind='barh',colormap='Spectral',ax=ax,log=True)
    plt.tight_layout()
    fig.savefig('srnabench_top_known.png')
    fig = plot_read_count_dists(final)
    fig.savefig('srnabench_known_counts.png')
    fig,ax = plt.subplots(figsize=(10,6))
    k[idcols].sum().plot(kind='bar',ax=ax)
    fig.savefig('srnabench_total_persample.png')
    print
    k.to_csv('srnabench_known_all.csv',index=False)
    return k

def get_top_isomirs(iso):
    """Get top isomir per mirRNA"""

    g = iso.groupby('name', as_index=False)
    #first add col for each isomirs percentage in its group
    totals = g.agg({'total':np.sum}).set_index('name')
    iso['perc'] = iso.apply(lambda r: r.total/totals.ix[r['name']].total,1)
    t=[]
    for i,x in g:
        r = base.first(x)
        s = x.total.sum()
        perc = r.total/s
        t.append((r['name'],r.read,r.mean_norm,r.total,s,perc,np.size(x.total),r.variant,r.pos))
    cols=['name','read','counts','mean_norm','total','domisoperc','isomirs','variant','pos']
    top = pd.DataFrame(t,columns=cols)
    top = top.sort('counts',ascending=False)
    return top

def analyse_isomirs(iso,outpath=None):
    """Analyse isomiR results in detail"""

    if iso is None:
        return
    if outpath != None:
        os.chdir(outpath)
    subcols = ['name','read','isoClass','NucVar','total','freq']
    iso = iso.sort('total', ascending=False)
    #filter very low abundance reads
    iso = iso[(iso.total>10) & (iso.freq>0.5)]
    top = get_top_isomirs(iso)
    top.to_csv('srnabench_isomirs_dominant.csv',index=False)
    print ('top isomiRs:')
    print (top[:20])
    print ('%s/%s with only 1 isomir' %(len(top[top.domisoperc==1]),len(top)))
    print ('different dominant isomir:', len(top[top.variant!='exact'])/float(len(top)))
    print ('mean dom isomir perc:', top.domisoperc.mean())
    print
    #stats
    fig,ax = plt.subplots(1,1)
    top.plot('isomirs','total',kind='scatter',logy=True,logx=True,alpha=0.8,s=50,ax=ax)
    ax.set_title('no. isomiRs per miRNA vs total adundance')
    ax.set_xlabel('no. isomiRs')
    ax.set_ylabel('total reads')
    fig.savefig('srnabench_isomir_counts.png')
    fig,ax = plt.subplots(1,1)
    #top.hist('domisoperc',bins=20,ax=ax)
    try:
        base.sns.distplot(top.domisoperc,bins=15,ax=ax,kde_kws={"lw": 2})
    except:
        pass
    fig.suptitle('distribution of dominant isomiR share of reads')
    fig.savefig('srnabench_isomir_domperc.png')

    x = iso[iso.name.isin(iso.name[:28])]
    bins=range(15,30,1)
    ax = x.hist('length',bins=bins,ax=ax,by='name',sharex=True,alpha=0.9)
    ax[-1,-1].set_xlabel('length')
    fig.suptitle('isomiR length distributions')
    fig.savefig('srnabench_isomir_lengths.png')
    plt.close('all')

    c=iso.variant.value_counts()
    #c=c[c>10]
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    c.plot(kind='pie',colormap='Spectral',ax=ax, labels=None,legend=True,
             startangle=0,pctdistance=1.1,autopct='%.1f%%',fontsize=10)
    ax.set_title('isomiR class distribution')
    plt.tight_layout()
    fig.savefig('srnabench_isomir_classes.png')
    fig,axs = plt.subplots(2,2,figsize=(8,6))
    grid=axs.flat
    bins=np.arange(-6,6,1)
    i=0
    for v in ['lv3p','lv5p','nta#A','nta#T']:
        ax=grid[i]
        iso[iso.variant==v].hist('pos',ax=ax,bins=bins,grid=False)
        ax.set_title(v)
        ax.set_xticks(bins+0.5)
        ax.set_xticklabels(bins)
        i+=1
    ax.set_xlabel('length difference')
    fig.suptitle('isomiR class length distributions')
    plt.tight_layout()
    fig.savefig('srnabench_isomir_variantlengths.png')
    #plt.show()
    iso.to_csv('srnabench_isomirs_all.csv',index=False)
    return top

def get_file_ids(path):
    """Get file-id mapping"""

    idmap = pd.read_csv(os.path.join(path, 'srnabench_colnames.csv'))
    return idmap

def combine_labels(labels, filename):
    """Combine sample labels with col names file so we can match sample ids"""

    comb = pd.read_csv(filename)
    comb['name'] = comb.filename.apply (lambda x: x.split('.')[0])
    labels = labels.merge(comb,on='name')
    labels['id'] = labels.id+'_norm'
    return labels

def main():
    try:
        base.seabornsetup()
    except:
        pass
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
            config.write_default_config('srnabench.conf', defaults=srnabenchoptions)
        cp = config.parse_config(opts.config)
        options = cp._sections['base']
        print (options)
        if opts.input == None:
            opts.input = options['path']
        run_all(**options)
    elif opts.analyse != None:
        k,n,iso = get_results(opts.analyse)
        analyse_results(k,n)
        analyse_isomirs(iso)

if __name__ == '__main__':
    main()
