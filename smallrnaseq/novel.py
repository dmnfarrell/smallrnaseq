#!/usr/bin/env python

"""
    Novel miRNA prediction
    Created Feb 2017
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
import shutil, glob, collections
import itertools
from itertools import islice
import numpy as np
import pandas as pd
from . import base, utils

def get_triplets(seq, struct):
    """triplet elements"""

    tr=['(((', '((.', '(..', '(.(', '.((', '.(.', '..(' , '...']
    nuc = ['A','G','T','C']
    d={}
    for i in nuc:
        for j in tr:
            d[i+j]=0
    struct=struct.replace(')','(')
    l = len(seq)-len(seq)%3
    for i in range(0,l,3):
        n = seq[i+1]+struct[i:i+3]
        if n in d:
            d[n]+=1
    return d

def get_biggest_stem(bg):
    biggest_stem = (-1, 'x')
    for s in bg.stem_iterator():
        if bg.stem_length(s) > biggest_stem[0]:
            biggest_stem = (bg.stem_length(s), s)
    return biggest_stem

def get_stem_pairs(bg):
    pairs=[]
    for s in bg.stem_iterator():
        for p in bg.stem_bp_iterator(s):
            pairs.append( (bg.seq[p[0]-1],bg.seq[p[1]-1]) )
    return pairs

def get_stem_matches(bg):
    pairs = get_stem_pairs(bg)
    wc = {'G':'C','C':'G','T':'A','A':'T'}
    matches = [True if i[1]==wc[i[0]] else False for i in pairs]
    return matches

def build_rna_features(seq, mature=None):
    """Get features for mirna sequence"""

    from Bio.SeqUtils import GC
    import forgi.graph.bulge_graph as cgb

    struct,sc = utils.rnafold(seq)

    feats = {}
    #feats['reads'] = reads.reads.sum()
    feats['length'] = len(seq)
    feats['mfe'] = round(sc/len(seq),3)
    #print seq
    #print struct

    bg = cgb.BulgeGraph()
    bg.from_dotbracket(struct)
    bg.seq = seq
    try:
        h0seq = bg.get_define_seq_str('h0')[0]
        feats['loops'] = len(list(bg.hloop_iterator()))
        feats['loop_length'] = bg.get_bulge_dimensions('h0')[0]
    except:
        h0seq=''
        feats['loops'] = 0
        feats['loop_length'] = 0
    feats['loop_gc']= GC(h0seq)
    feats['stem_length'] = len(get_stem_pairs(bg))
    feats['longest_stem'] = get_biggest_stem(bg)[0]
    bulges = [bg.get_bulge_dimensions(i) for i in bg.iloop_iterator()]
    feats['bulges'] = len(bulges)
    try:
        feats['longest_bulge'] =  max(max(zip(*bulges)))
    except:
        feats['longest_bulge'] = 0
    bulgematches = [True if i[0]==i[1] else False for i in bulges]
    feats['bulges_symmetric'] = bulgematches.count(True)
    feats['bulges_asymmetric'] = bulgematches.count(False)
    sm = get_stem_matches(bg)
    feats['stem_mismatches'] = sm.count(False)

    if mature == None:
        #mature should be given - place holder for code to guess it later
        start = np.random.randint(1,len(sm)-20)
    else:
        start = utils.find_subseq(seq, mature)
        end = start+len(mature)
        #print start, end
    feats['mature_mismatches'] = sm[start:end].count(False)

    tr = get_triplets(seq, struct)
    feats.update(tr)
    return feats

def get_star(seq, mature, struct=None):
    """Estimate the star sequence from a given mature and precursor."""

    import forgi.graph.bulge_graph as cgb
    start = utils.find_subseq(seq, mature)+1
    end = start + len(mature)
    if struct == None:
        struct,sc = utils.rnafold(seq)
    bg = cgb.BulgeGraph()
    bg.from_dotbracket(struct)
    bg.seq = seq
    stempairs = []
    for s in bg.sorted_stem_iterator():
        stempairs.extend( list(bg.stem_bp_iterator(s)) )
    m = zip(*stempairs)
    stem1 = list(m[0])
    stem2 = list(m[1])
    matidx = range(start, end)
    #print stem1
    #print start, end
    #is mature on 5' or 3' end?
    if start < max(stem1):
        print ('5p')
        matidx = [i for i in matidx if i in stem1]
        staridx = [i[1] for i in stempairs if i[0] in matidx]
        gaps = [abs(t-s) for s, t in zip(staridx, staridx[1:])]
        for i in range(len(gaps)):
            if gaps[i]>3 and i>=len(gaps)-5:
                staridx = staridx[:i+1]
        offset = len(matidx)-len(staridx)+2
        starseq = seq[staridx[-1]:staridx[0]+offset]
    else:
        print ('3p')
        matidx = [i for i in matidx if i in stem2]
        staridx = [i[0] for i in stempairs if i[1] in matidx]
        offset = len(matidx)-len(staridx)+2
        #print matidx
        #print staridx
        starseq = seq[staridx[0]+offset:staridx[-1]]

    #print matidx
    #print staridx
    #print mature
    #print starseq
    return starseq

def get_positives(species='hsa'):
    """Get known mirbase hairpins for training precursor classifier. """

    reload(base)
    mirs = base.get_mirbase(species)
    feats=[]
    for i,row in mirs.iterrows():
        f = build_rna_features(row.precursor, row.mature1_seq)
        f['seq'] = row.precursor
        f['mature'] = row.mature1_seq
        f['star'] = row.mature2_seq
        feats.append(f)

    known = pd.DataFrame(feats)
    known.to_csv('known_mirna_features.csv')
    return known

def get_negatives():
    """negative pseudo mirna set"""

    cds = utils.fasta_to_dataframe('../genomes/human/Homo_sapiens.GRCh38.cds.all.fa')
    cds = cds.drop_duplicates('sequence')
    cds = cds[cds.sequence.str.len()>50]

    def split_seqs(r):
        #maxlen = int(np.random.normal(81,17))
        maxlen = int(np.random.gamma(9.5,9))
        #print len(r.sequence), maxlen
        s = [r.sequence[ind:ind+maxlen] for ind in range(0, len(r.sequence), maxlen)]
        return pd.Series(s)

    seqs = cds[:2000].apply(split_seqs,1).stack().reset_index(drop=True)
    seqs = seqs[seqs.str.len()>50]
    seqs = seqs[-seqs.str.contains('N')]
    result=[]
    for seq in seqs:
        ms = int(np.random.randint(2,5))
        f = build_rna_features(seq, mature=seq[ms:ms+22])
        f['seq'] = seq
        result.append(f)
    result = pd.DataFrame(result)
    result = result[(result.loops==1) & (result.mfe*result.length<=-15) & (result.stem_length>18)]
    return result

def precursor_classifier(known, neg):
    """Train a mirna precursor classifier"""

    import sklearn
    from sklearn.ensemble import (RandomForestClassifier, ExtraTreesClassifier, RandomForestRegressor)
    from sklearn.model_selection import train_test_split,cross_val_score

    print (len(known), len(neg))
    known['target'] = 1
    neg['target'] = 0
    data = pd.concat([known,neg]).reset_index(drop=True)
    data = data.sample(frac=1)

    y = data.target
    data = data.drop('target',1)
    X = data.select_dtypes(['float','int'])
    #print data[:5]
    #X = sklearn.preprocessing.scale(X)

    #test
    #X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.4)
    #rf = RandomForestClassifier()
    rf = RandomForestRegressor()
    #rf.fit(X_train, y_train)
    #y_score = rf.predict(X_test)
    scores = cross_val_score(rf, X, y, cv=5, scoring='roc_auc')
    print (scores)
    #print sklearn.metrics.classification_report(y_test, y_score)

    rf = RandomForestClassifier()
    #rf = RandomForestRegressor()
    rf.fit(X,y)
    names = data.columns
    importances = rf.feature_importances_
    indices = np.argsort(importances)[::-1]
    #print indices
    print ('feature ranking:')
    for f in range(X.shape[1])[:10]:
        print("%d. %s (%f)" % (f + 1, names[indices[f]], importances[indices[f]]))
    return rf

def score_features(data, rf):
    """Score a set of features"""

    X = data.select_dtypes(['float','int'])
    #data['score'] = rf.predict(X)
    return rf.predict(X)

def build_cluster_trees(alnmt, cluster_distance=100, min_size=2, key='read_id'):
    """Build cluster tree of reads from a dataframe of locations e.g from
        a set of aligned reads from a sam file.
    Args:
        cluster_distance: Distance in basepairs for two reads to be in the same cluster;
       for instance 20 would group all reads with 20bp of each other
        min_size: Number of reads necessary for a group to be considered a cluster;
       2 returns all groups with 2 or more overlapping reads
    Returns:
        dict of ClusterTrees per chromosome
    """

    import collections
    from bx.intervals.cluster import ClusterTree
    cluster_trees = collections.defaultdict(lambda:
            ClusterTree(cluster_distance, min_size))
    for i, row in alnmt.iterrows():
        chrom = row['name']
        #print chrom, row.read_id, row.start, row.end
        cluster_trees[chrom].insert(row.start, row.end, row[key])
    return dict(cluster_trees)

def get_read_clusters(reads):
    """get clusters of reads from a dataframe with alignment fields
      i.e. from samfile info"""

    df = reads
    df = df[df.length<=26]
    #get clusters of reads and store by read_id
    clustertrees = build_cluster_trees(df, 10, 5)
    df.set_index('read_id',inplace=True)

    groups = []
    i=1
    for chrom, cltree in clustertrees.items():
        #print chrom
        for start, end, ids in cltree.getregions():
            #print start, end, ids
            c = df.ix[ids]
            c['cl_start'] = start
            c['cl_end'] = end
            c['cluster'] = i
            groups.append(c)
            i+=1
    df = pd.concat(groups)
    return df

def find_precursor(ref_fasta, cluster, cluster2=None, step=5, score_cutoff=1):
    """Find most likely precursor from a genomic sequence and
       one or two mapped read clusters.
       Args:
           ref_fasta: genomic reference sequence
           cluster: reads in a cluster, a dataframe
           cluster2: a pair to the first cluster, optional
           step: increment for extending precursors
       Returns:
           the top precursor
    """

    x = cluster.iloc[0]
    maturereads = reads1 = cluster.reads.sum()
    mature = x.seq
    star = None
    starreads = 0
    #dtermine mature/star if two clusters present
    if cluster2 is not None:
        reads2 = cluster2.reads.sum()
        y = cluster2.iloc[0]
        if reads2>reads1:
            mature = y.seq
            maturereads = reads2
            star = x.seq
            starreads = reads1

    print (mature, star, maturereads, starreads)
    #check mature for non templated additions?

    chrom = x['name']
    strand = x.strand
    loop = 15
    N = []
    #generate candidate precursors
    for i in range(1,45,step):
        #5' side
        start5 = x.start - i
        end5 = x.start + 2 * len(x.seq)-1 + loop + i
        coords = [chrom,start5,end5,strand]
        prseq = utils.sequence_from_coords(ref_fasta, coords)
        N.append({'precursor':prseq, 'chrom':chrom,'start':start5,'end':end5,'mature':x.seq,'strand':strand})
        #3' side
        start3 = x.start - (loop + len(x.seq) + i)
        end3 = x.end + i
        coords = [chrom,start3,end3,strand]
        prseq = utils.sequence_from_coords(ref_fasta, coords)
        N.append({'precursor':prseq, 'chrom':chrom,'start':start3,'end':end3,'mature':x.seq,'strand':strand})

    N = pd.DataFrame(N)
    print (len(N))
    N['mature_reads'] = maturereads
    N['star_reads'] = starreads
    if cluster2 is not None:
        N['star'] = cluster2.iloc[0].seq
        #'check star for non templated additions?

    f = N.apply(lambda x: pd.Series(build_rna_features(x.precursor, x.mature)), 1)

    N['score'] = score_features(f, rf)
    N['mfe'] = f.mfe
    #filter by feature
    N = N[(f.loops==1) & (f.stem_length>18) & (f.mfe*f.length<-15)]
    N = N[N.score>=score_cutoff]
    #print N
    N = N.sort_values('mfe')
    if len(N)>0:
        found = N.iloc[0]
        return found
    else:
        return

def find_mirnas(samfile, ref_fasta, truecounts):
    """Find novel miRNAs in reference mapped reads. Assumes we have already
        mapped to known miRNAs. """

    reads = utils.get_aligned_reads(samfile, truecounts)
    df = get_read_clusters(reads)

    print ('%s unique reads in clusters' %len(df))
    clusts = df.groupby(['name','cluster','cl_start','cl_end','strand'])\
                            .agg({'reads':np.sum,'length':np.max})\
                            .reset_index()\
                            .rename(columns={'cl_start':'start','cl_end':'end'})
    print ('%s read clusters' %len(clusts))

    #find pairs of read clusters - likely to be mature/star sequences
    clustpairs = build_cluster_trees(clusts, 120, min_size=2, key='cluster')

    def get_pairs(r):
        #get ids for corresponding pairs
        id = r.cluster
        pairs = clustpairs[r['name']].getregions()
        if len(pairs)>0 and id in pairs[0][2]:
            x = pairs[0][2]
            return x[0]
        return

    def get_coords(r):
        return r['chrom']+':'+str(r.start)+'..'+str(r.end)+':'+r.strand

    clusts['pair'] = clusts.apply(get_pairs, 1)
    #print clusts

    n1 = []
    for i,r in clusts.groupby('pair'):
        a,b = list(r.cluster)
        c1 = df[df.cluster==a]
        c2 = df[df.cluster==b]
        print ('paired clusters found')
        p = find_precursor(ref_fasta, c1, c2, step=7)
        if p is None:
            print ('no precursor predicted')
            continue
        n1.append(p)
    n1 = pd.DataFrame(n1)
    print
    n2 = []
    #guess precursors for single clusters
    singleclusts = clusts[clusts.pair.isnull()]
    for i,r in islice(singleclusts.iterrows(),5,20):
        c = df[df.cluster==r.cluster]
        p = find_precursor(ref_fasta, c, step=7)
        #print p
        if p is None:
            print ('no precursor predicted')
            continue
        #estimate star sequence
        p['star'] = get_star(p.precursor, p.mature)
        n2.append(p)
    n2 = pd.DataFrame(n2)

    novel = pd.concat([n1,n2])
    #get seed seq and mirbase matches
    novel['seed'] = novel.apply(lambda x: x.mature[2:8], 1)
    #get coords column
    novel['coords'] = novel.apply(get_coords,1)

    print ('found %s novel mirnas' %len(novel))
    return novel
