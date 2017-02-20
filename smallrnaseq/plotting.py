#!/usr/bin/env python

"""
    Plotting methods for smallrnaseq package
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
import sys, os, string, types
import itertools
import matplotlib
import pylab as plt
import numpy as np
import pandas as pd
from . import utils, analysis

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
    return v

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

    df = utils.fastq_to_dataframe(filename)
    x = analysis.read_length_dist(df)
    fig,ax=plt.subplots(1,1,figsize=(10,4))
    ax.bar(x[1][:-1],x[0], align='center')
    return ax

def plot_sample_variation(df):

    fig,axs=plt.subplots(2,1,figsize=(6,6))
    axs=axs.flat
    cols,ncols = mirdeep2.get_column_names(m)
    x = m.ix[2][cols]
    x.plot(kind='bar',ax=axs[0])
    x2 = m.ix[2][ncols]
    x2.plot(kind='bar',ax=axs[1])
    sns.despine()
    plt.tight_layout()
    return

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

def plot_fractions(df, label=None, path=None):
    """Process results of multiple mappings to get fractions
    of each annotations mapped
    label: plot this sample only"""

    if len(df.columns) == 1:
        label = df.columns[0]
    if label != None:
        explode = [0.05 for i in range(len(df))]
        axs = df.plot(y=label,kind='pie',colormap='Spectral',autopct='%.1f%%',
                      startangle=0,figsize=(6,6),
                      labels=None,legend=True,pctdistance=1.1,
                      explode=explode,fontsize=10)
    else:
        df = df.set_index('label')
        df = df._get_numeric_data()
        l = df.plot(kind='bar',stacked=True,cmap='Spectral',figsize=(12,6))
        plt.legend(ncol=4)
        plt.title('rna fractions mapped')

    plt.tight_layout()
    if path == None:
        path='.'
    plt.savefig(os.path.join(path,'fractions_mapped.png'))
    return