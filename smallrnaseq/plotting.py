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
    Plotting methods for smallrnaseq package.
    Created June 2014
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, string, types
import itertools
import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style("ticks", {'axes.facecolor': '#F7F7F7',
                        'axes.grid': False,'legend.frameon':True})
sns.set_context("notebook", font_scale=1.3)

from . import base, utils, analysis

def venn_diagram(names,labels,ax=None,**kwargs):
    """Plot venn diagrams"""

    from matplotlib_venn import venn2,venn3
    f=None
    if ax==None:
        f=plt.figure(figsize=(4,4))
        ax=f.add_subplot(111)
    if len(names)==2:
        n1,n2=names
        v = venn2([set(n1), set(n2)], set_labels=labels, **kwargs)
    elif len(names)==3:
        n1,n2,n3=names
        v = venn3([set(n1), set(n2), set(n3)], set_labels=labels, **kwargs)
    ax.axis('off')
    #f.patch.set_visible(False)
    ax.set_axis_off()
    return

def heatmap(df,fname=None,cmap='seismic',log=False):
    """Plot a heat map"""

    from matplotlib.colors import LogNorm
    f=plt.figure(figsize=(8,8))
    ax=f.add_subplot(111)
    norm=None
    df=df.replace(0,.1)
    if log==True:
        norm=LogNorm(vmin=df.min().min(), vmax=df.max().max())
    hm = ax.pcolor(df,cmap=cmap,norm=norm)
    plt.colorbar(hm,ax=ax,shrink=0.6,norm=norm)
    plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    #ax.axvline(4, color='gray'); ax.axvline(8, color='gray')
    plt.tight_layout()
    if fname != None:
        f.savefig(fname+'.png')
    return ax

def plot_read_lengths(filename, df=None):
    """View read length distributions"""

    df = utils.fastq_to_dataframe(filename, size=5e5)
    x = analysis.read_length_dist(df)
    fig,ax=plt.subplots(1,1,figsize=(10,4))
    ax.bar(x[1][:-1],x[0], align='center')
    return fig

def plot_sample_variation(df):

    fig,axs=plt.subplots(2,1,figsize=(6,6))
    axs=axs.flat
    cols,ncols = mirdeep2.get_column_names(m)
    x = m.ix[2][cols]
    x.plot(kind='bar',ax=axs[0])
    x2 = m.ix[2][ncols]
    x2.plot(kind='bar',ax=axs[1])
    sns.despine(trim=True,offset=10)
    plt.tight_layout()
    return fig

def plot_by_label(X, palette='Set1'):
    """Color scatter plot by dataframe index label"""

    import seaborn as sns
    cats = X.index.unique()
    colors = sns.mpl_palette(palette, len(cats))
    #sns.palplot(colors)
    f,ax = plt.subplots(figsize=(6,6))
    for c, i in zip(colors, cats):
        #print X.ix[i,0]
        ax.scatter(X.ix[i, 0], X.ix[i, 1], color=c, s=100, label=i,
                   lw=1, edgecolor='black')
    ax.legend(fontsize=10)
    sns.despine()
    return

def plot_fractions(df, label=None):
    """Process results of multiple mappings to get fractions
    of each annotations mapped
    label: plot this sample only"""

    fig,ax = plt.subplots(figsize=(8,8))
    df = df.set_index('label')
    df = df._get_numeric_data()
    if len(df) == 1:
        label = df.index[0]
    if label != None:
        ax = df.T.plot(y=label,kind='pie',colormap='Spectral',autopct='%.1f%%',
                      startangle=0, labels=None,legend=True,pctdistance=1.1,
                      fontsize=10, ax=ax)
    else:
        ax = df.plot(kind='barh',stacked=True,linewidth=0,cmap='Spectral',ax=ax)
        #ax.legend(ncol=2)
        ax.set_position([0.2,0.1,0.6,0.8])
        ax.legend(loc="best",bbox_to_anchor=(1.0, .9))
    plt.title('fractions mapped')
    #plt.tight_layout()
    return fig

def plot_sample_counts(counts):

    fig,ax = plt.subplots(figsize=(10,6))
    scols,ncols = base.get_column_names(counts)
    counts[scols].sum().plot(kind='bar',ax=ax)
    plt.title('total counts per sample (unnormalised)')
    plt.tight_layout()
    return fig

def plot_read_count_dists(counts, h=8, n=50):
    """Boxplots of read count distributions """

    scols,ncols = base.get_column_names(counts)
    df = counts.sort_values(by='mean_norm',ascending=False)[:n]
    df = df.set_index('name')[ncols]
    t = df.T
    w = int(h*(len(df)/60.0))+4
    fig, ax = plt.subplots(figsize=(w,h))
    if len(scols) > 1:
        sns.stripplot(data=t,linewidth=1.0,palette='coolwarm_r')
        ax.xaxis.grid(True)
    else:
        df.plot(kind='bar',ax=ax)
    sns.despine(offset=10,trim=True)
    ax.set_yscale('log')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylabel('read count')
    #print (df.index)
    #plt.tight_layout()
    fig.subplots_adjust(bottom=0.2,top=0.9)
    return fig

def expression_clustermap(counts, freq=0.8):

    scols,ncols = base.get_column_names(counts)
    X = counts.set_index('name')[ncols]
    X = np.log(X)
    v = X.std(1).sort_values(ascending=False)
    X = X[X.isnull().sum(1)/len(X.columns)<0.2]
    X = X.fillna(0)
    cg = sns.clustermap(X,cmap='YlGnBu',figsize=(12,12),lw=0,linecolor='gray')
    mt = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=9)
    mt = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    return cg

def plot_PCA(X, cmap='Spectral', colors=None, dims=(0,1), ax=None, annotate=None, legend=True, **kwargs):
    '''
    plot PCA from matrix and label
        :X:  dataframe with index as categories
        :dims:  dimensions to plot
        :return: None
    '''

    from sklearn import preprocessing
    from sklearn.decomposition.pca import PCA
    X = X._get_numeric_data()
    S = pd.DataFrame(preprocessing.scale(X),columns = X.columns)
    pca = PCA(n_components=4)
    pca.fit(S)
    out = 'explained variance %s' %pca.explained_variance_ratio_
    print (out)
    #print pca.components_
    w = pd.DataFrame(pca.components_,columns=S.columns)
    #print (w.T.max(1).sort_values())
    pX = pca.fit_transform(S)
    pX = pd.DataFrame(pX,index=X.index)

    ### graph
    if ax is None:
        fig,ax=plt.subplots(1,1,figsize=(8,8))

    cats = pX.index.unique()
    if colors is None:
        colors = sns.mpl_palette(cmap, len(cats))

    y1,y2=dims
    offset = 7
    for c, i in zip(colors, cats):
       ax.scatter(pX.loc[i, y1], pX.loc[i, y2], color=c, label=i, edgecolor='black', **kwargs)

    if annotate is not None:
       pX['lab#el'] = annotate
       i=0
       for idx,r in pX.iterrows():
          x=r[y1]; y=r[y2]
          l=annotate[i]
          ax.annotate(l, (x,y), xycoords='data',xytext=(2,5),textcoords='offset points',fontsize=12)
          i+=1

    ax.set_xlabel("X[%s]" %y1)
    ax.set_ylabel("X[%s]" %y2)
    if legend == True:
       ax.set_position([0.1,0.1,0.5,0.8])
       ax.legend(loc="best",bbox_to_anchor=(1.0, .9))
    ax.set_title("PCA")
    return pX
