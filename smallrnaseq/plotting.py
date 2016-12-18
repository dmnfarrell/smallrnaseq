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