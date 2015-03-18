#!/usr/bin/env python

"""Module for methods using ensembl with pycogent api
   Created September 2014
   Copyright (C) Damien Farrell
"""

import numpy as np
import pandas as pd
import cogent
from cogent.db.ensembl import HostAccount, Genome, Compara
import base

account = None
release = '78'
species = ['cow','human','mouse','rat','chimp','gorilla','orangutan',
           'macaque','dog','pig','cat','olive baboon','sheep']
#add recent ensembl species
from cogent.db.ensembl import Species
Species.amendSpecies('Papio anubis', 'olive baboon')
Species.amendSpecies('Ovis aries', 'sheep')
Species.amendSpecies('Erinaceus europaeus', 'hedgehog')
Species.amendSpecies('Mustela putorius furo', 'ferret')


def getOrthologs(refgenome, ensid=None, symbol=None):
    '''if ensid!=None:
        mygene = cow.getGeneByStableId(StableId=ensid)
    else:
        mygene = cow.getGenesMatching(symbol=symbol)'''
    #print mygene
    orthologs = comp.getRelatedGenes(gene_region=mygene,
                    Relationship='ortholog_one2one')
    print orthologs.Members
    #seqs = orthologs.getSeqCollection(feature_types='gene')
    #print seqs.Names
    #syntenicregions = comp.getSyntenicRegions(region=mygene,
    #                align_method='EPO', align_clade='eutherian')
    return

def getGenesFromLocation(ref, coords, pad=0):
    """Get genes from a set of coords.
       pad will add n bases to either side to expand area"""

    genome = Genome(Species=ref, Release=release, account=account)
    print genome
    chrom,start,end,strand = coords
    genes = list(genome.getFeatures(CoordName=chrom, Start=start-pad, End=end+pad, feature_types='gene'))
    return genes

def getAlignmentTree(fname):
    """Build a neighbour joining tree"""

    from cogent.phylo import distance, nj
    from cogent.evolve.models import HKY85,F81
    al = cogent.LoadSeqs(fname,format='fasta')
    d = distance.EstimateDistances(al, submodel= F81())
    d.run(show_progress=False)
    mytree = nj.nj(d.getPairwiseDistances())
    mytree = mytree.balanced()
    print mytree.asciiArt()
    print
    '''from cogent.draw import dendrogram
    p = dendrogram.SquareDendrogram(mytree)
    p.drawToPDF('tree-scaled.pdf', 500, 400, stroke_width=2.0,
                shade_param = 'r', max_value = 1.0,)'''
    return

def getSyntenicAlignment(comp, ref, coords, fname='ensembl.aln.fa'):
    """Get ensembl genomic alignment for a location using compara and return seqs"""

    print
    chrom,start,end,strand = coords
    regions = comp.getSyntenicRegions(Species=ref, CoordName=chrom,
                                       Start=start, End=end, align_method="EPO_LOW_COVERAGE",
                                       align_clade="39 eutherian", Strand=strand)
    regions = [r for r in regions]
    if len(regions)==0:
        print 'no alignments'
        return None, None
    #usually only 1 alignment
    A = regions[0].getAlignment()
    A.writeToFile(fname)
    print cogent.LoadSeqs(fname)
    #if len(A.Seqs)>3:
    #    getAlignmentTree(fname)
    if regions is not None:
        print '%s syntenic regions found' %len(regions)
    return regions, A

def getmiRNAOrthologs(df, comp=None, ref='cow'):
    """Get all possible orthologs/conservation for miRNAs using ensembl"""

    if comp == None:
        comp = Compara(species, account=None, Release='79')
    results=[]
    for i, r in list(df.iterrows()):
        #base.RNAfold(r['consensus precursor sequence'], r['#miRNA']+'_'+ref)
        mature = r['consensus mature sequence'].replace('u','t').upper()
        star = r['consensus star sequence'].replace('u','t').upper()
        seed = r['seed'].replace('u','t').upper()
        mbmatch = r['mirbase seed match']
        print r['#miRNA'], seed, mature, star, mbmatch
        c,locs,strand = r['precursor coordinate'].split(':')
        start,end = locs.split('..')
        coords = (c,int(start),int(end),strand)
        regions, aln = getSyntenicAlignment(comp, ref, coords, fname=r['#miRNA']+'.aln.fa')

        if aln != None:
            a = base.cogentAlignment2DataFrame(aln.degap())
            a['#miRNA'] = r['#miRNA']
            a['seed'] = seed
            a['ident'] = getIdentities(aln)
            print 'max identity: %s' %a.ident.max()
            a['seedcons'] = getSeqConservation(aln, seed)
            a['mirbase'] = mbmatch
            orthgenes = getGenesinRegion(regions[0])
            targets=[]
            #find possible targets in 3'UTR of containing gene
            for g in orthgenes:
                try:
                    tr = g[0].CanonicalTranscript
                    sd = seed.replace('u','t').upper()
                    l = str(tr.Utr3.complement()).find(sd)
                    targets.append(l)
                except:
                    targets.append(-1)

            a['genes'] = [g[0].Symbol if len(g)>0 else np.nan for g in orthgenes]
            a['targets'] = targets
            a['location'] = getLocations(regions[0])

            #get RNAfold energy for each sequence
            a['energy'] = a.apply(lambda x : base.RNAfold(x.seq)[1],1)
            results.append(a)
            print a
        print '--------------------------------------------------------'
    results = pd.concat(results).reset_index(drop=True)
    results.to_csv('novel_orthologs.csv')
    return

def getLocations(region):
    """Locations with region alignments"""
    l=[]
    #need try except because test for None fails.. why
    for r in region.Members:
        try:
            loc=r.Location
            #remove sp name
            loc = ':'.join(str(loc).split(':')[2:])
            l.append(loc)
        except:
            continue
    return l
    #return [r.Location for r in region.Members if r is not None]

def getSeqConservation(aln, seq):
    """Get position of a sub sequence and return position if it's
        conserved for each aligned species, -1 if not"""

    vals=[]
    for s in aln.Seqs:
        print seq, s, str(s).find(seq)
        vals.append(str(s).find(seq))
    return vals

def getIdentities(aln):
    """Identities for all seqs in a cogent alignment"""

    names=aln.Names
    ref=names[0]
    vals=[None]
    for n in names[1:]:
        new = aln.takeSeqs([ref, n])
        ident = round(len(new.filtered(lambda x: len(set(x)) == 1))/float(len(new)),3)
        vals.append(ident)
    return vals

def getGenesinRegion(region):
    """Get orthologous genes in all species in a syntenic region"""

    print region.ref_location
    print 'checking genes in aligned species'
    orthologs=[]
    for r in region.Members:
        try:
            loc = r.Location
        except:
            continue
        sp = loc.Species
        coords = loc.CoordName,loc.Start,loc.End,loc.Strand
        #print sp, loc
        #genes = getGenesFromLocation(sp, coords, pad=1000)
        genes = list(r.genome.getFeatures(feature_types='gene',region=r))
        #genes = [g.Symbol for g in genes]
        orthologs.append(genes)
    return orthologs

def getESTs(region):
    for r in region.Members:
        loc = r.Location
        sp = loc.Species
        #genome = Genome(Species=sp, Release=release, account=account)
        ests = r.genome.getFeatures(feature_types='est', region=r)
        for est in ests:
            print est
    return

def summarise(df):
    """Summarise candidates from ensembl results"""

    n = pd.read_csv('novel_orthologs.csv')
    n = n[n.ident!=0]
    x = n.groupby('#miRNA').agg({
                    'seq':np.size,
                    'energy':np.min,
                    'ident': np.max,
                    'seedcons': lambda r: len(r[r>-1]),
                    'mirbase': base.first,
                    'genes': base.first,
                    'targets': lambda r: len(r[r>-1])})
    x.columns = x.columns.get_level_values(0)
    x = x.merge(df[['#miRNA','read_count','miRDeep2 score','freq','precursor coordinate',
                    'seed','consensus mature sequence']],
                left_index=True,right_on='#miRNA',how='outer')
    #print x
    x=x.set_index('#miRNA')
    #entries with at least 1 alignment and conserved seed are best candidates
    def isconserved(x):
        return (x.seq>1) & (x.seedcons>=2)
    #x['conserved'] = x.apply(isconserved,1)
    x = x.sort(['seq'],ascending=False)
    #x=x.fillna('')
    x.to_csv('novel_conserved.csv',float_format='%2.2f')
    return x

def test():

    comp = Compara(species, account=account, Release=release)
    #print comp.method_species_links
    coords = [('15', 85273455, 85273507, '+'),('18',12423976,12424060,'+')]
    for c in coords:
        getSyntenicAlignment(comp, 'cow', c)
        #getAlignmentTree('ensembl_aln.fa')
    return

