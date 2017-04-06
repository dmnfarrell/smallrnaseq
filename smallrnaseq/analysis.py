#!/usr/bin/env python

"""
    smallrnaseq additional analysis routines
    Created June 2014
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
import sys, os, string, types, re, csv
import shutil, glob, collections
import itertools
import subprocess
import matplotlib
import pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
try:
    import HTSeq
except:
    'HTSeq not present'

from . import mirdeep2 as mdp
from . import base, utils

def read_length_dist(df):

    df['length'] = df.seq.str.len()
    bins = np.linspace(1,df.length.max(),df.length.max())
    x = np.histogram(df.length,bins=bins)
    return x

def summarise_reads(path):
    """Count reads in all files in path"""

    resultfile = os.path.join(path, 'read_stats.csv')
    files = glob.glob(os.path.join(path,'*.fastq'))
    vals=[]
    rl=[]
    for f in files:
        label = os.path.splitext(os.path.basename(f))[0]
        s = utils.fastq_to_dataframe(f)
        l = len(s)
        vals.append([label,l])
        print (label, l)

    df = pd.DataFrame(vals,columns=['path','total reads'])
    df.to_csv(resultfile)
    df.plot(x='path',y='total reads',kind='barh')
    plt.tight_layout()
    plt.savefig(os.path.join(path,'total_reads.png'))
    #df = pd.concat()
    return df

def remove_known_rnas(path, adapters=[], outpath='RNAremoved'):
    """Map to annotated RNAs and put remaining reads in output dir"""

    index = 'bosTau6-tRNAs'
    params = '-v 0 --best'
    files = glob.glob(os.path.join(path,'*.fastq'))
    #files = ['test.fastq']
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for f in files:
        print (f)
        label = os.path.splitext(os.path.basename(f))[0]
        trimAdapters(f, adapters)
        fastq2fasta('cut.fastq')
        rem = os.path.join(outpath, label+'.fa')
        #samfile = os.path.join(outpath, '%s_%s_mapped.sam' %(label,index))
        samfile = 'out.sam'
        base.bowtie_align('cut.fa', index, outfile=samfile, params=params, remaining=rem)

    return

def ks_test(df, ids):
    """KS test to determine rc freq distributions of replicates
    threshold count is the one which corresponds to the first minimum.
    i.e.  when the distributions of reads amongst replicates begin to be similar
    see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2822534/
    ids are pairs of replicate column ids for mirna results
    """

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
    readcuts = np.arange(0,300,2)
    result=[]
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
            #print len(x),len(y),rc, val
            stats.append(val[0])
        result.append(stats)
    result = pd.DataFrame(np.array(result).T)
    result=result.set_index(readcuts)
    result.plot(lw=2,colormap='Set1',legend=False)
    tmins = result.apply(getMin)
    print (tmins)
    mean=tmins.mean()
    std=tmins.std()

    plt.axvspan(mean-std/2,mean+std/2,color='g',alpha=0.3)
    plt.xlabel('read count threshold')
    plt.ylabel('KS')
    plt.savefig('KS_test.png',dpi=150)

    return

def mirna_discovery_test(sourcefile):
    """Test miRNAs found per read file size. Assumes runs have been
     done and placed in benchmarking dir"""

    path = 'benchmarking'
    #sizes = np.arange(5e5,7.e6,5e5)
    #sizes = np.insert(sizes,0,[5e4,1e5,2e5])
    #print sizes
    #base.createRandomFastqFiles(sourcefile, path, sizes)

    #mirdeep
    outpath = 'benchmarking/mirdeep'
    a = mdp.getResults(outpath)
    k1 = a[a.novel==False]
    k1=k1.set_index('#miRNA')
    k1 = k1[k1.s16>5]
    mp1 = mdp.get_file_ids(outpath)
    cols,ncols = mdp.get_column_names(k1)
    #srnabench
    outpath = 'benchmarking/srnabench'
    #k2,n,iso = srb.get_results(outpath)
    k2 = k2.set_index('name')
    k2 = k2[k2.s16>5]
    mp2 = srb.getFileIDs(outpath).sort('filename')

    def getFound(k,mapping,rc=5):
        tot = float(len(k))
        f=[]
        results=[]
        for i,r in mapping.iterrows():
            fn = r.filename
            i = float(re.findall("\d+.\d+",fn)[0])
            df = k[k[r.id]>rc]
            n = df[-df.index.isin(f)]
            f.extend(n.index)
            results.append((i,len(f),len(f),n.total.max()))
        r = pd.DataFrame(results,columns=['reads','total','frac','rcm'])
        return r

    r1 = getFound(k1, mp1)
    r2 = getFound(k2, mp2)

    fig=plt.figure(figsize=(8,6))
    ax=fig.add_subplot(111)
    l1=ax.plot(r1.reads,r1.frac,'bo-',color='blue',lw=2)
    l2=ax.plot(r2.reads,r2.frac,'bo-',color='red',lw=2)
    ax.set_xlabel('reads (million)')
    ax.set_ylabel('total miRNA found')
    ax.legend(l1+l2,['mirdeep2','srnabench'],loc=2)
    ax3 = plt.axes([.5, .25, .3, .3])
    ax3.semilogy(r1.reads,r1.rcm,color='blue',ls='-')
    ax3.semilogy(r2.reads,r2.rcm,color='red',ls='-')
    ax3.set_title('mean abundance (new miRNAs)')
    plt.tight_layout()
    fig.savefig('benchmark_discovery.png')

    '''fig,ax=plt.subplots(1,1)
    x=k1[cols]
    from matplotlib.colors import LogNorm
    x=x.replace(0,.1)
    ax.pcolor(x,cmap='Blues',
              norm=LogNorm(vmin=x.min().min(), vmax=x.max().max()))
    plt.yticks(np.arange(0.5, len(x.index), 1), x.index)
    plt.xticks(np.arange(0.5, len(x.columns), 1), x.columns)'''

    plt.show()
    return

def do_pca(X, c=3):
    """Do PCA"""

    from sklearn import preprocessing
    from sklearn.decomposition.pca import PCA, RandomizedPCA
    #do PCA
    #S = standardize_data(X)
    S = pd.DataFrame(preprocessing.scale(X),columns = X.columns)
    pca = PCA(n_components=c)
    pca.fit(S)
    print (pca.explained_variance_ratio_)
    #print pca.components_
    w = pd.DataFrame(pca.components_,columns=S.columns)#,index=['PC1','PC2'])
    #print w.T.max(1).sort_values()
    pX = pca.fit_transform(S)
    pX = pd.DataFrame(pX,index=X.index)
    return pX

def plot_pca(pX, palette='Spectral', labels=False, ax=None, colors=None):
    """Plot PCA result, input should be a dataframe"""

    if ax==None:
        fig,ax=plt.subplots(1,1,figsize=(6,6))
    cats = pX.index.unique()
    colors = sns.mpl_palette(palette, len(cats)+1)
    print (len(cats), len(colors))
    for c, i in zip(colors, cats):
        #print (i, len(pX.ix[i]))
        #if not i in pX.index: continue
        ax.scatter(pX.ix[i, 0], pX.ix[i, 1], color=c, s=90, label=i,
                   lw=.8, edgecolor='black', alpha=0.8)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    if labels == True:
        for i, point in pX.iterrows():
            ax.text(point[0]+.3, point[1]+.3, str(i),fontsize=(9))
    ax.legend(fontsize=10,bbox_to_anchor=(1.5, 1.05))
    sns.despine()
    plt.tight_layout()
    return

def do_mds(X):
    """Do MDS"""

    from sklearn import manifold
    seed = np.random.RandomState(seed=3)
    mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9, random_state=seed,
                        n_jobs=1)
    pX = mds.fit(X.values).embedding_
    pX = pd.DataFrame(pX,index=X.index)
    return pX

def get_aligned_reads_lengths(path, label, refs):
    """Get read length distributions mapped to sam files representing
      different ref annotations. Returns dataframe of read length freqs
      with column label for db mapped to"""

    cols=[]
    for ref in refs:
        f = os.path.join(path,'%s_%s.sam' %(label,ref))
        reads = base.get_aligned(f)
        x = read_length_dist(reads)
        x = pd.Series(x[0],index=x[1][:-1],name=ref)
        x = x/x.sum()
        #print x
        cols.append(x)
    x = pd.DataFrame(cols).T
    x.index = x.index.astype(int)
    x.index.name = 'length'
    return x

def classify(X, y, cl, name=''):
    """Classification using gene features"""

    from sklearn.metrics import classification_report, accuracy_score
    np.random.seed()
    ind = np.random.permutation(len(X))

    from sklearn.cross_validation import train_test_split
    Xtrain, Xtest, ytrain, ytest  = train_test_split(X, y, test_size=0.4)
    #print X
    cl.fit(Xtrain, ytrain)
    ypred = cl.predict(Xtest)

    print (classification_report(ytest, ypred))
    #print accuracy_score(ytest, ypred)
    from sklearn import cross_validation
    yl = pd.Categorical(y).labels
    sc = cross_validation.cross_val_score(cl, X, yl, scoring='roc_auc', cv=5)
    print("AUC: %0.2f (+/- %0.2f)" % (sc.mean(), sc.std() * 2))
    return cl

def read_length_distributions(path, refs):
    """Get read lengths for all aligned files in path mapped
       to ref db names, assumes the fasta files and sam files are present in the
       target folder
    """

    files = glob.glob(path+'/*.fa')
    fnames = [os.path.splitext(os.path.basename(f))[0] for f in files]
    fnames = [i for i in fnames if not i.endswith('_r')]
    print (fnames)

    res = []
    for n in fnames:
        x = get_aligned_reads_lengths(path, n, refs)
        x = x.ix[10:40]
        x['file'] = n
        res.append(x)
    res = pd.concat(res)
    res = res.reset_index()
    #print res[:60]
    #plot histograms
    m = pd.melt(res, id_vars=['file','length'], value_vars=refs,
                var_name='ref', value_name='freq')
    g = sns.factorplot(y='freq',x='length', data=m, row='ref', kind="bar",
                        size=2,aspect=4,sharex=False,palette='Blues_d')
    plt.savefig(os.path.join(path,'read_lengths.png'))
    return res

def get_trna_fragments(samfile, fastafile, truecounts, bedfile=None):
    """Get trf5/3/i fragments from a set reads aligned to a trna sequences.
    This finds the locations of all reads inside the aligned parent and
    classifies the fragments using a scheme similar to Telonis at al."""

    refs = utils.fasta_to_dataframe(fastafile)
    #get rna structure for reference seqs?
    #refs['structure'] = refs.apply(lambda r: utils.rnafold(r.sequence)[0],1)
    a = utils.get_aligned_reads(samfile, truecounts)
    #a = a.merge(truecounts, on='seq')
    print ('%s total sequences with %s counts' %(len(a),a.reads.sum()))
    if len(a) == 0:
        return

    def get_type(x, refs):
        #classify by position in parent seqs
        n = x['name']
        seq = refs.ix[n].sequence
        if x.start<=2:
            return 'trf5'
        elif len(seq)-x.end == 0 and x.seq.endswith('CCA'):
            return 'trf3'
        else:
            return 'itrf'

    def get_loop(x):
        if x.start>11 and x.start<17:
            return 'D'
        elif x.start>30 and x.start<40:
            return 'A'

    #remove sequence redundancy by grouping into unique fragments then get classes
    #first sort reads by sequence and name so we get consistent ids for different samples..
    a = a.sort_values(['seq','name'])
    f = a.groupby('seq').agg({'name':base.first, 'reads':np.max, 'read':np.size,
                                 'start':base.first,'end':base.first})
    f=f.rename(columns={'read':'loci'})
    f = f.reset_index()
    f['anticodon'] = f.apply(lambda x: x['name'].split('-')[1],1)
    f['aa'] = f.anticodon.str[:3]
    f['perc'] = (f.reads/f.reads.sum()*100).round(3)
    f['length'] = f.seq.str.len()

    f['trf'] = f.apply(lambda x: get_type(x, refs), 1)
    f['loop'] = f.apply(lambda x: get_loop(x), 1)
    f['id'] = f.apply(lambda x:'%s_%s.%s.%s' % (x['name'],x.length,x.start,x.end),1)

    print ('%s unique fragments' %len(f))
    print ('%s total counts in fragments' %f.reads.sum())
    return f

if __name__ == '__main__':
    main()
