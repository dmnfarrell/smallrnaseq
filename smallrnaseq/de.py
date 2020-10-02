#!/usr/bin/env python

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
Methods for differential expression analysis
Created Dec 2016
Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, string, types, re, csv
import shutil, glob
import itertools
import subprocess
import numpy as np
import pandas as pd
from . import base, utils

def get_columns_by_label(labels, samplecol, filters=[], querystr=None):
    """Get sample columns according to a condition from a set of labels
    Args:
        labels: dataframe matching sample labels to conditions/factors
        samplecol: name of column holding sample/file names
        filters: tuples containing column/.value pairs to filter on
        querystr: optional string instead of tuple filters
        (see pandas.DataFrame.query documentation)
    """

    if querystr == None:
        q=[]
        for f in filters:
            print (f)
            dtype = labels[f[0]].dtype.kind
            if dtype in 'bifc':
                s = "%s==%s" %(f[0],f[1])
            else:
                s = "%s=='%s'" %(f[0],f[1])
            q.append(s)
        querystr = ' & '.join(q)
    print (querystr)
    x = labels.query(querystr)
    cols = x[samplecol]
    return list(cols)

def get_factor_samples(df, labels, factors, filters=[], samplecol='filename', index=None):
    """Get subsets of samples according to factor/levels specified in another mapping file
       Used for doing differential expr with edgeR.
       Args:
            labels: dataframe matching sample labels to conditions/factors
            factors: conditions to compare
            filters: additional filters for samples e.g. a time point
            samplecol: name of column holding sample/file names
       Returns:
            dataframe of data columns with header labels for edgeR script
            a dict mapping condition to column names
    """

    if index != None:
        df=df.set_index(index)
    l=0
    res = None
    samples = []
    for f in factors:
        cond = f[1]
        f = filters + [f]
        cols = get_columns_by_label(labels, samplecol, f)
        print (cols)
        cols = list(set(cols) & set(df.columns))
        samples.extend(cols)
        x = df[cols]
        print ('%s samples, %s genes' %(len(cols),len(x)))
        if len(x.columns)==0:
            #no data found warning
            print ('WARNING: no data for %s' %f)
            continue
        print()
        x.columns = ['s'+str(cols.index(i))+'_'+str(l) for i in x.columns]
        l+=1
        if res is None:
            res = x
        else:
            res = res.join(x)
    res=res.dropna()
    return res, samples

def run_edgeR(countsfile=None, data=None, cutoff=1.5):
    """Run edgeR from R script"""

    if data is not None:
        countsfile = 'de_counts.csv'
        data.to_csv(countsfile)
    path = os.path.dirname(os.path.abspath(__file__)) #path to module
    descript = os.path.join(path, 'DEanalysis.R')
    cmd = 'Rscript %s %s' %(descript, countsfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print (result)
    #read result back in
    de = pd.read_csv('edger_output.csv')
    de.rename(columns={'Unnamed: 0':'name'}, inplace=True)
    de = de[(de.FDR<0.05) & ((de.logFC>cutoff) | (de.logFC<-cutoff))]
    return de

def run_limma(countsfile=None, data=None, cutoff=1.5):
    """Run limma de from R script"""

    if data is not None:
        countsfile = 'de_counts.csv'
        data.to_csv(countsfile)
    path = os.path.dirname(os.path.abspath(__file__)) #path to module
    descript = os.path.join(path, 'Limma.R')
    cmd = 'Rscript %s %s' %(descript, countsfile)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    #read result back in
    de = pd.read_csv('limma_output.csv')
    de.rename(columns={'Unnamed: 0':'name'}, inplace=True)
    #md_plot(data, de)
    de = de[(de['adj.P.Val']<0.05) & ((de.logFC>cutoff) | (de.logFC<-cutoff))]
    de = de.sort_values('logFC',ascending=False)
    return de

def md_plot(data, de, title=''):
    """MD plot"""

    data = data.reset_index()
    data['mean log expr'] = data.mean(1).apply(np.log)
    df = data.merge(de,on='name')
    #print (df[:10])
    a = df[df['adj.P.Val']<0.05]
    b = df[-df.name.isin(a.name)]
    c = a[a.logFC<0]
    ax=a.plot('mean log expr','logFC',kind='scatter',figsize=(10,10),color='red',s=60)
    c.plot('mean log expr','logFC',kind='scatter',ax=ax,color='g',s=60)
    b.plot('mean log expr','logFC',kind='scatter',ax=ax,color='black')
    ax.set_title(title, fontsize=20)
    return ax

def runEdgeRGLM(countsfile, cutoff=1.5):
    """Run edgeR from R script"""

    cmd = 'Rscript ~/python/sandbox/mirnaseq/GLMDEanalysis.R %s' %countsfile
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print (result)
    #read result back in
    #de = pd.read_csv('de_output.csv')
    #de.rename(columns={'Unnamed: 0':'name'}, inplace=True)
    #de = de[(de.FDR<0.05) & ((de.logFC>cutoff) | (de.logFC<-cutoff))]
    return

def rpyEdgeR(data, groups, sizes, genes):
    """Run edgeR analysis - from http://bcbio.wordpress.com/ """

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    robjects.r('''library(edgeR)''')
    params = {'group' : groups, 'lib.size' : sizes}
    print (params)
    d = robjects.r.DGEList(counts=data, **params)
    print (d)
    robjects.r.calcNormFactors(d)
    robjects.r.estimateCommonDisp(d)
    robjects.r.estimateTagwiseDisp(d)
    robjects.r.exactTest(d)

    #ms = robjects.r.deDGE(dgelist, doPoisson=True)
    #tags = robjects.r.topTags(ms, pair=groups, n=len(genes))
    indexes = [int(t) - 1 for t in tags.rownames()]
    pvals = list(tags.r['adj.P.Val'][0])
    return

def melt_samples(df, labels, names, samplecol='filename', index=None):
    """Melt sample data by factor labels so we can plot with seaborn"""

    if index is not None:
        df=df.set_index(index)
    df = df.loc[names]
    t=df.T
    t.index = df.columns
    t = t.merge(labels,left_index=True,right_on=samplecol)
    m = pd.melt(t,id_vars=list(labels.columns),
                 var_name='name',value_name='read count')
    return m

def cluster_map(data, names):
    """Cluster map of genes"""

    import seaborn as sns
    import pylab as plt
    data = data.loc[names]
    X = np.log(data).fillna(0)
    X = X.apply(lambda x: x-x.mean(), 1)
    cg = sns.clustermap(X,cmap='RdYlBu_r',figsize=(8,10),lw=.5,linecolor='gray')
    mt=plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    mt=plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    return cg
