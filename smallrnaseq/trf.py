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
    smallrnaseq tRNA fragments analysis.
    Created June 2017
    Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, string, types, re, csv
import itertools
import subprocess
import numpy as np
import pandas as pd
from . import base, utils

def get_anticodon(x):
    s = x['first'].split('.')[1].split('-')[1]
    return s

def get_trna_families(ref_fasta):

    trnas = utils.fasta_to_dataframe(ref_fasta).reset_index()
    #print trnas[:3]
    g = trnas.groupby('sequence').agg({'name':[np.size,base.first]})
    g.columns = g.columns.get_level_values(1)
    g['ac'] = g.apply(lambda x: get_anticodon(x) , 1)
    g = g.reset_index()
    g['id'] = g.groupby('ac').cumcount()+1
    g['family'] = g.apply(lambda x: x.ac+'-'+str(x.id)+'-'+str(x['size']), 1)
    print (len(g))
    refname = os.path.splitext(ref_fasta)[0]
    utils.dataframe_to_fasta(g,'%s-fam.fa' %refname,idkey='family',seqkey='sequence')
    return

def tdr_mapper(samfile, collapsed, ref_trnas, threshold=20):
    """Get trf5/3/i fragments from a set reads aligned to a trna sequences.
    This finds the locations of primary trfs inside each aligned parent 'family' trna and
    classifies the fragments using a scheme similar to tdrmapper"""

    refs = utils.fasta_to_dataframe(ref_trnas)
    #print samfile, collapsed
    a = utils.get_aligned_reads(samfile, collapsed)
    a = a[a.reads>threshold]
    total = float(a.drop_duplicates('seq').reads.sum())
    print ('%s total sequences with %s counts' %(len(a),a.reads.sum()))
    if len(a) == 0:
        return

    def overlap(start, end, x1, x2):
        if ((start<x1) & (end>x1)) or ((start>x1) & (start<x2)):
            return True

    def pos_coverage(r, p):
        x = [r.reads if (i>=r.start and i<=r.end) else 0 for i in p]
        return pd.Series(x,index=p)

    #find primary trna by getting coverage over each trna, so group by gene
    grps = a.groupby('name')
    f = []
    for name,df in grps:
        parent = refs.ix[name]
        tlen = len(parent.sequence)
        p = range(1,tlen)
        m = df.apply( lambda x: pos_coverage(x,p), 1 )
        cov = m.sum()/df.reads.sum()
        pr = cov[cov>=.5]
        if len(pr) == 0:
            continue
        start, end = pr.index[0],pr.index[-1]
        seq = parent.sequence[start-1:end]
        l = len(seq)
        reads = df[(df['start']>=start-1) & (df.end<=end+1)].reads.sum()
        #relative abundance
        readcov = round(reads/float(df.reads.sum()),2)

        if l<41 and l>=28:
            frtype = 'tRH'
        elif l>14 and l<28:
            frtype = 'tRF'
        else:
            continue
        region = ''
        if start == 1:
            region = '5'
        elif tlen-end < 6:
            region = '3'
        else:
            if overlap(start, end, 13,22) == True:
                region = 'D'
            if overlap(start, end, 31,39) == True:
                region += 'A'
            t1=cov.index[-23]; t2=cov.index[-15]
            if overlap(start, end, t1,t2) == True:
                region += 'T'
        if region == '':
            continue
        f.append( [name, frtype, region, start, end, seq, reads, readcov] )

    f = pd.DataFrame(f, columns=['family','frtype','region','start','end','seq','reads','coverage'])
    f['anticodon'] = f.apply(lambda x: x.family.split('-')[0], 1)
    f['aa'] = f.anticodon.str[:3]
    f['length'] = f.seq.str.len()
    f['id'] = f.apply(lambda x: x.family+'-'+x.frtype+'-'+x.region, 1)
    f['abundance'] = (f.reads/total*100).round(4)
    f = f[f.coverage>=.6]
    f = f[f.length>15]
    f = f.sort_values('reads',ascending=False)#.set_index('id')
    s = f.groupby('seq').first()

    print ('%s primary tdrs, %s unique sequences' %(len(f), len(s)))
    return f
