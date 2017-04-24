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

path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(path, 'data')
CLASSIFIER = None

def get_triplets(seq, struct):
    """Get triplet elements used by Xue et al."""

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
    """Stem pairs of rna hairpin"""

    pairs=[]
    for s in bg.stem_iterator():
        for p in bg.stem_bp_iterator(s):
            pairs.append( (bg.seq[p[0]-1],bg.seq[p[1]-1]) )
    return pairs

def get_stem_matches(bg):
    """Find rna stem mismatches"""

    pairs = get_stem_pairs(bg)
    wc = {'G':'C','C':'G','T':'A','A':'T'}
    matches = [True if i[1]==wc[i[0]] else False for i in pairs]
    return matches

def build_rna_features(seq, struct=None, sc=None, mature=None):
    """Get features for mirna sequence"""

    from Bio.SeqUtils import GC
    import forgi.graph.bulge_graph as cgb

    if struct == None:
        struct,sc = utils.rnafold(seq)

    feats = {}
    #feats['reads'] = reads.reads.sum()
    feats['length'] = len(seq)
    feats['mfe'] = round(sc/len(seq),3)
    #print seq
    #print struct

    bg = utils.get_bg(seq, struct)
    try:
        h0seq = bg.get_define_seq_str('h0')[0]
        feats['loops'] = len(list(bg.hloop_iterator()))
        feats['loop_length'] = bg.get_bulge_dimensions('h0')[0]
    except:
        h0seq=''
        feats['loops'] = 0
        feats['loop_length'] = 0
    feats['gc'] = GC(seq)
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
        end = start+22
    else:
        start = utils.find_subseq(seq, mature)
        end = start+len(mature)
        #print start, end
    feats['mature_mismatches'] = sm[start:end].count(False)
    feats['struct'] = struct
    tr = get_triplets(seq, struct)
    feats.update(tr)
    return feats

def find_star_sequence(seq, mature, struct=None):
    """Estimate the star sequence from a given mature and precursor."""

    bg = utils.get_bg(seq, struct)
    start = utils.find_subseq(seq, mature)+1
    end = start + len(mature)
    #print (start, end)
    #find end of mature to estimate start of star
    for p in bg.adjacent_stem_pairs_iterator():
        #print (p)
        for s in p:
            #print (s)
            x = list(bg.stem_bp_iterator(s))
            #print (x)
            for i in x:
                if start == i[0]:
                    #print ('found start', i )
                    star = seq[i[1]-len(mature):i[1]+2]
                    return star
                elif start == i[1]:
                    #print ('found start', i )
                    star = seq[i[0]-len(mature)+1:i[0]]
                    return star
    print ('star not found')
    return

def check_hairpin(seq, struct):
    """Remove free ends of hairpin"""

    bg = utils.get_bg(seq, struct)
    loops = list(bg.hloop_iterator())
    if len(loops) != 1:
        return seq, struct
    p = list(bg.stem_bp_iterator('s0'))[0]
    s = p[0]-1; e=p[1]
    return seq[s:e], struct[s:e]

def check_mature(seq, struct, mature):
    """Check if the mature sequence is not in hairpin loop and inside stem"""

    bg = utils.get_bg(seq, struct)
    start = utils.find_subseq(seq, mature)+1
    end = start + len(mature)
    loops = list(bg.hloop_iterator())
    if len(loops)==0:
        return 'too many loops'
    l = list(bg.define_range_iterator('h0'))[0]
    if (start<l[0] and end>l[0]) or (start<l[1] and end>l[1]):
        return 'mature in loop'
    p = list(bg.stem_bp_iterator('s0'))[0]
    #print (p[0], p[1], start, end)
    if start == 0:
        return 'mature not found'
    if end>p[1]+2:
        return '3p mature outside'
    return 'ok'

def get_positives(species='hsa'):
    """Get known mirbase hairpins for training precursor classifier. """

    reload(base)
    mirs = base.get_mirbase(species)
    feats=[]
    for i,row in mirs.iterrows():
        f = build_rna_features(row.precursor, mature=row.mature1_seq)
        f['seq'] = row.precursor
        f['mature'] = row.mature1_seq
        f['star'] = row.mature2_seq
        feats.append(f)

    result = pd.DataFrame(feats)
    result.to_csv('known_mirna_features.csv', index=False)
    return result

def get_negatives(fasta_file=None, samples=2000):
    """Create a negative pseudo mirna set for training classifier"""

    if fasta_file == None:
        cds = utils.fasta_to_dataframe('../genomes/human/Homo_sapiens.GRCh38.cds.all.fa')
    cds = cds.drop_duplicates('sequence')
    cds = cds[cds.sequence.str.len()>50]

    def split_seqs(r):
        #maxlen = int(np.random.normal(81,17))
        maxlen = int(np.random.gamma(9.5,9))
        #print len(r.sequence), maxlen
        s = [r.sequence[ind:ind+maxlen] for ind in range(0, len(r.sequence), maxlen)]
        return pd.Series(s)

    seqs = cds[:5000].apply(split_seqs,1).stack().reset_index(drop=True)
    seqs = seqs[seqs.str.len()>50]
    seqs = seqs[-seqs.str.contains('N')]
    result=[]
    i=1
    for seq in seqs:
        ms = int(np.random.randint(2,5))
        f = build_rna_features(seq, mature=seq[ms:ms+22])
        f['seq'] = seq
        if f['loops'] > 2:# or f['stem_length']>18:
            continue
        result.append(f)
        i+=1
        if i > samples:
            break
    result = pd.DataFrame(result)
    #result = result[(result.mfe*result.length<=-15)]
    result.to_csv('negative_mirna_features.csv', index=False)
    return result

def get_training_data(known=None, neg=None):
    """Get training data for classifier
        Args:
            known: known precursor data, a dataframe
            neg: negatives
        Returns:
            a dataframe with all features and a set of true/false values
    """

    if known is None:
        known = pd.read_csv(os.path.join(datadir, 'training_positives.csv'))
        #known = get_positives()
    if neg is None:
        neg = pd.read_csv(os.path.join(datadir, 'training_negatives.csv'))
    #print (len(known), len(neg))
    known['target'] = 1
    neg['target'] = 0
    data = pd.concat([known,neg]).reset_index(drop=True)
    data = data.sample(frac=1)
    y = data.target
    data = data.drop('target',1)
    X = data.select_dtypes(['float','int'])
    #X = sklearn.preprocessing.scale(X)
    return X, y

def precursor_classifier(known=None, neg=None, kind='regressor'):
    """Get a miRNA precursor classifier using given training data.
       Args:
        X: numpy array/dataframe with features
        y: true/false values matching feature rows, 1's and 0's
        kind: use 'classifier' or 'regressor' random forest
       Returns:
        random forest classifier fitted to X,y
    """

    X, y = get_training_data(known, neg)
    from sklearn.ensemble import (RandomForestClassifier, RandomForestRegressor)
    if kind == 'classifier':
        rf = RandomForestClassifier(n_estimators=100)
    else:
        rf = RandomForestRegressor(n_estimators=100)
    #print ('fitting..')
    rf.fit(X,y)
    return rf

def test_classifier(known=None, neg=None):

    from sklearn.ensemble import (RandomForestClassifier, RandomForestRegressor)
    X, y = get_training_data(known, neg)
    rf = RandomForestClassifier(n_estimators=50)
    rf.fit(X,y)
    from sklearn.model_selection import train_test_split,cross_val_score
    #X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.4)
    #rf = RandomForestRegressor()
    scores = cross_val_score(rf, X, y, cv=5, scoring='roc_auc')
    print (scores)
    #print sklearn.metrics.classification_report(y_test, y_score)
    names = X.columns
    importances = rf.feature_importances_
    indices = np.argsort(importances)[::-1]
    print ('feature ranking:')
    for f in range(X.shape[1])[:10]:
        print("%d. %s (%f)" % (f + 1, names[indices[f]], importances[indices[f]]))
    '''a = neg[:2000].drop('target',1)
    a['score'] = score_features(a, rf)
    b = known[:2000].drop('target',1)
    b['score'] = score_features(b, rf)
    x = a.score.value_counts().sort_index()
    y = b.score.value_counts().sort_index()
    res = pd.DataFrame({'neg':x,'pos':y})
    res.plot(kind='bar')
    '''
    return

def score_features(data, rf):
    """Score a set of rna features"""

    X = data.select_dtypes(['float','int'])
    #data['score'] = rf.predict(X)
    return rf.predict(X)

def build_cluster_trees(reads, cluster_distance=2, min_size=2):
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
    for i, row in reads.iterrows():
        chrom = row['name']
        #print chrom, row.read_id, row.start, row.end
        cluster_trees[chrom].insert(row.start, row.end, row.name)
    return dict(cluster_trees)

def get_read_clusters(reads, cluster_distance=0, min_size=3, key='align_id'):
    """Assign reads to clusters from a dataframe from an alignment
      Args:
        reads: pandas dataframe of reads with start, end, read_id fields
        key: key to use for assigning reads to clusters, should be unique
      Returns:
        the dataframe with cluster numbers assigned
    """

    reads = reads.copy()
    if key in reads.columns:
        reads.set_index(key, inplace=True)
    #build clustertrees per chromosome of reads and store by index values
    clustertrees = build_cluster_trees(reads, cluster_distance, min_size)

    groups = []
    i=1
    for chrom, cltree in clustertrees.items():
        #print (chrom)
        for start, end, ids in cltree.getregions():
            #print (start, end, ids)
            c = reads.ix[ids].copy()
            c['cl_start'] = start
            c['cl_end'] = end
            c['cluster'] = i
            #remove wrongly added opposite strand reads
            c = c.groupby(['strand']).filter(lambda x: len(x) > 1)
            groups.append(c)
            i+=1
    df = pd.concat(groups)
    return df

def generate_precursors(ref_fasta, coords, mature=None, step=5):
    """Create a set of possible precursor sequences from flanking sequence given
       genomic coordinates and the reference genome sequence.
       Args:
        ref_fasta: reference genome fasta file
        cooords: coordinates of estimated mature sequence to start from
        seq: mature sequence
       Returns:
        dataframe with precursor sequences and their coordinates
    """

    chrom,start,end,strand=coords
    loop = 15
    N = []
    if mature != None:
        seqlen = len(mature)
    else:
        seqlen = 22
    #generate candidate precursors
    for i in range(2,40,step):
        #5' side
        start5 = start - i
        end5 = start + 2 * seqlen-1 + loop + i
        coords = [chrom,start5,end5,strand]
        prseq = utils.sequence_from_coords(ref_fasta, coords)
        if prseq == None:
            continue
        struct,sc = utils.rnafold(prseq)
        #prseq, struct = check_hairpin(prseq, struct)
        mstatus = check_mature(prseq, struct, mature)
        #print (mstatus)
        #print (i)
        #print (prseq)
        #print (struct)
        N.append({'precursor':prseq,'struct':struct,'score':sc,
                  'chrom':chrom,'start':start5,'end':end5,
                  'mature':mature,'strand':strand, 'mature_check': mstatus})
    for i in range(2,40,step):
        #3' side
        start3 = start - (loop + seqlen + i)
        end3 = end + i + seqlen+1
        coords = [chrom,start3,end3,strand]
        prseq = utils.sequence_from_coords(ref_fasta, coords)
        if prseq == None:
            continue
        struct,sc = utils.rnafold(prseq)
        #prseq, struct = check_hairpin(prseq, struct)
        mstatus = check_mature(prseq, struct, mature)
        N.append({'precursor':prseq,'struct':struct,'score':sc,
                  'chrom':chrom,'start':start3,'end':end3,
                  'mature':mature,'strand':strand,'mature_check': mstatus})
    N = pd.DataFrame(N)
    return N

def score_precursors(N):
    """Filter/score a set of precursor sequences, requires a dataframe
      with a column of precursor sequences and structures made by the
      generate_precursors method. """

    global CLASSIFIER
    if CLASSIFIER == None:
        rf = precursor_classifier()
    else:
        rf = CLASSIFIER
    f = N.apply( lambda x: pd.Series(build_rna_features(x.precursor,
                                                        x.struct, x.score, x.mature)), 1 )
    #N['struct'] = f.struct
    #check mature out of stem range also
    #check mature for non templated additions?

    N['score'] = score_features(f, rf)
    N['mfe'] = f.mfe
    #filter by features
    N = N[(f.loops==1) & (f.stem_length>18) & (f.mfe*f.length<-15) & \
          (f.longest_bulge<10) & (f.gc<75) & (f.bulges_asymmetric<7)]

    N = N.sort_values('mfe')
    #print (N)
    return N

def get_consensus_read(ref, df):
    """Get consensus read from cluster of aligned reads"""

    strand = df.iloc[0].strand
    if strand == '+':
        s='start'
    else:
        s='end'
    g = df.groupby([s,'length']).agg({'reads':np.sum})\
        .sort_values(by='reads',ascending=False)\
        .reset_index()
    #print g
    x = g.iloc[0]
    cons = df[(df[s]==x[s]) & (df.length==x.length)].sort_values(by='reads',ascending=False)
    #print cons
    mature = None
    for i,r in cons.iterrows():
        if r.seq in ref:
            mature = r.seq
            break
    if mature==None:
        mature = cons.iloc[0].seq
    return mature

def find_precursor(ref_fasta, m, o=None, step=5, score_cutoff=.7):
    """Find the most likely precursor from a genomic sequence and
       one or two mapped read clusters.
       Args:
           ref_fasta: genomic reference sequence
           cluster: reads in a cluster, a dataframe
           cluster2: a pair to the first cluster, optional
           step: increment for extending precursors
           score_cutoff: if using non-classifier, optional
       Returns:
           the top precursor
    """

    x=m.iloc[0]
    rcoords = (x['name'], x.start-10, x.end+10, x.strand)
    refseq = utils.sequence_from_coords(ref_fasta, rcoords)
    #get a consensus mature sequence
    mature = get_consensus_read(refseq, m)

    coords = (x['name'], x.start, x.start, x.strand)
    N = generate_precursors(ref_fasta, coords, mature=mature, step=step)
    if len(N)==0:
        return
    N = score_precursors(N)
    #print ('%s candidates' %len(N))
    N = N[N.score>=score_cutoff]
    if len(N)>0:
        P = N.iloc[0].copy()
        #print(P)
    else:
        return

    #print (o)
    maturecounts = m.reads.sum()
    star = find_star_sequence(P.precursor, mature, P.struct)
    starcounts = 0
    if o is not None and star != None:
        #check reads from local cluster that are inside star seq?
        s = utils.find_subseq(P.precursor, star)
        #if P.strand == '+':
        ss = P.start+s; se = ss+len(star)
            #sreads = o[(o.start>=ss-2) & (o.end<=se+3)]
        #else:
            #se = P.end-s; ss = se-len(star)
        sreads = o[(o.start>=ss-2) & (o.end<=se+3)]
        starcounts = sreads.reads.sum()
        #print (ss, se)
        #print (sreads)
        #print display(HTML(forna_url(P.precursor, mature, star)))

    print (P.start, P.end, P.strand)
    print (mature, maturecounts, star, starcounts)
    P['mature_reads'] = maturecounts
    P['star_reads'] = starcounts
    P['star'] = star
    P['cluster'] = x.cluster
    #print (mature, star, maturecounts, starcounts)
    #print ('')
    return P

def find_mirnas(reads, ref_fasta, score_cutoff=.8, read_cutoff=50, species='',
                max_length=25, min_length=18, min_size=3):
    """Find novel miRNAs in reference mapped reads. Assumes we have already
        mapped to known miRNAs.
        Args:
            reads: unique aligned reads with counts in a dataframe
            ref_fasta: reference genome fasta file
            score_cutoff: max score to keep for precursors
            species: three letter mirbase code for species, optional
        Returns:
            tuple of dataframes: novel mirnas and the read clusters
    """

    global CLASSIFIER
    if CLASSIFIER == None:
        print ('getting default classifier')
        CLASSIFIER = precursor_classifier(kind='regressor')

    reads = reads[(reads.length<=max_length) & (reads.length>=min_length)]
    #assign reads to clusters
    rcl = get_read_clusters(reads, 10, min_size)
    clusts = rcl.groupby(['name','cluster','cl_start','cl_end','strand'])\
                            .agg({'reads':np.sum,'length':np.max})\
                            .reset_index()\
                            .rename(columns={'cl_start':'start','cl_end':'end'})
    print ('%s read clusters in %s reads' %(len(clusts),len(rcl)))
    clusts['clust_size'] = clusts.end-clusts.start
    clusts = clusts[clusts.reads>=read_cutoff]
    print ('%s clusters above reads cutoff' %len(clusts))

    X = []
    N = []
    for i,c in clusts.iterrows():
        df = rcl[rcl.cluster==c.cluster]
        if df.reads.sum() < read_cutoff:
            continue
        df = df.sort_values('reads',ascending=False)
        #small clusters should belong to a single mature
        if c.clust_size<28:
            df['mature'] = True
            X.append(df)
            p = find_precursor(ref_fasta, df, score_cutoff=score_cutoff)
            #print p
            if p is not None:
                N.append(p)
                #print (p.mature)
            continue

        #for larger clusters choose most likely mature reads
        anchor = df.iloc[0]
        st = anchor.start
        end = anchor.end
        m = df.loc[(abs(df.start-st)<=3) & (abs(df.end-end)<=5)].copy()
        m['mature'] = True
        X.append(m)
        #remainder of reads assigned as non-mature
        o = df.loc[-df.index.isin(m.index)].copy()
        o['mature'] = False
        X.append(o)

        p = find_precursor(ref_fasta, m, o, score_cutoff=score_cutoff)

        #check for nearby clusters inside this precursor
        #nc = clusts[(abs(clusts.start-c.start)<100) & (clusts.cluster!=c.cluster)]

        if p is not None:
            N.append(p)
            print (p.mature)

    new = pd.DataFrame(N)
    new['seed'] = new.apply(lambda x: x.mature[1:7], 1)
    #get coords column
    new['coords'] = new.apply(get_coords_string,1)
    new = new.reset_index(drop=True)
    assign_names(new, species)
    kp = base.get_mirbase(species)
    #check if the mature are already in known precursors
    new['known_id'] = new.apply( lambda x: find_from_known_prec(x, kp), 1)
    new = new.sort_values(by='mature_reads', ascending=False)

    print ('found %s novel mirnas' %len(new))
    print ('%s with known mature sequences' %len(new[-new.known_id.isnull()]))
    #also return all the reads found in clusters
    found = pd.concat(X)
    return new, found

def get_coords_string(r):
    """coords string from fields"""
    if 'chrom' not in r:
        r['chrom']=r['name']
    return r['chrom']+':'+str(r.start)+'..'+str(r.end)+':'+r.strand

def find_from_known_prec(x, prec):
    for i,r in prec.iterrows():
        if r.precursor.find(x.mature[:-1]) != -1:
            return r.mirbase_id

def assign_names(df, species=''):
    """Assign name to novel mirna, precursor/mature ids should allow consistent
       identification across datasets"""

    df['precursor_id'] = df.apply( lambda x: species+'_novel_'+x.chrom+'_'+str(x.start),1 )
    df['mature_id'] = df.apply( lambda x: species+'_'+encode_name(x.mature), 1 )
    return

def encode_name(s):
    """hash a sequence into a short string"""

    import hashlib
    h = hashlib.md5(s.encode())
    s = h.digest().encode('base64')[:8]
    s.replace('/','x')
    return s

def forna_url(precursor, mature, star=None, struct=None):
    """Create link to view mirna structure in forna web app"""

    #print x
    pd.options.display.max_colwidth = 500
    #seq = x.precursor
    #mature = x[col1]
    mstart = utils.find_subseq(precursor, mature)+1
    mend = mstart+len(mature)-1
    if struct==None:
        struct = utils.rnafold(precursor)[0]
    colors = '%s-%s:lightgreen\\n' %(mstart, mend)
    if star != None:
        #print x[maturecol], x[col2]
        #star = x[col2]
        if star != None and not pd.isnull(star):
            sstart = utils.find_subseq(precursor, star)+1
            send = sstart+len(star)-1
            colors += '%s-%s:pink' %(sstart,send)
    url='http://nibiru.tbi.univie.ac.at/forna/forna.html'\
        '?id=url/name&sequence=%s&structure=%s&colors=%s' %(precursor,struct,colors)
    html = '<a href="%s" target="_blank"> view structure </a>' %url
    return html

def string_to_html(s):
    """Convert lines of strings for html rendering"""

    html=''
    x = s.split('\n')
    for line in x:
        line = line.replace(' ','&nbsp') #preserve spaces
        html += line+'<br>'
    return html

def create_report(df, reads, species=None, outfile='report.html'):
    """Novel miRNA predictions html report"""

    pd.options.display.max_colwidth = 500
    css = get_css()
    h = '<html><head><meta charset="utf-8">  <title>novel miRNA</title>'
    h += '<style media="screen" type="text/css"> %s </style>' %css
    h += '</head>'
    h += '<body>'
    h += '<div class="header">'
    h += '<h3>novel miRNA predictions</h3>'
    h += '</div>'
    h += '<div class="sidebar">'
    links = df[['mature_id','mature_reads','chrom']].copy()
    links['mature_id'] = links.mature_id.apply(lambda x: ('<a href=#%s > %s </a>' %(x,x)))
    h += links.to_html(escape=False, classes='sidebar', index=False)
    h += '</div>'

    df = df.copy()
    df = df.set_index('mature_id')

    ens_sp = pd.read_csv(os.path.join(datadir, 'ensembl_names.csv'), index_col=0)
    if species in ens_sp.index:
        ensname = ens_sp.ix[species]['scientific name']
        df['coords'] = df.coords.apply(
            lambda x: ('<a href=http://www.ensembl.org/%s/Location/View?r=%s target="_blank">%s </a>' %(ensname,x,x)),1)
    df['link'] = df.apply(lambda x: forna_url(x.precursor, x.mature, x.star, x.struct), 1)

    h += '<div class="content">'
    for i,r in df.iterrows():
        #print (r.name, r.mature, r.star)
        h += '<div class="box">'
        h += '<a name=%s></a>' %i
        h += r.to_frame().to_html(escape=False)
        x = reads[(reads.cluster==r.cluster)]
        #print (x)
        #print (x.reads.sum())
        s = utils.print_read_stack(x, r.precursor, by='reads')
        if s==None:
            h+='</div>'
            continue
        h += '<p>'
        h += string_to_html(s)
        h += '</p>'
        h += '</div>'

    h+='</div>'
    h += '</body>'
    f=open(outfile,'w')
    f.write(h)
    f.close()
    return h

def get_css():
    """Get css style for embedding in html page"""

    fname = os.path.join(datadir, 'styles.css')
    with open(fname) as f:
        content = f.readlines()
        content = ''.join(content)
    return content
