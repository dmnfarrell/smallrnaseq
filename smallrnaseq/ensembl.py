#!/usr/bin/env python

"""
    Module for methods using ensembl mostly using pycogent API
    Created September 2014
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
import numpy as np
import pandas as pd
from . import base

try:
    import cogent
    from cogent.db.ensembl import HostAccount, Genome, Compara, Species
    account = None
    release = '87'
    #add recent ensembl species
    species = ['cow','human','mouse','rat','chimp','gorilla','orangutan',
               'macaque','dog','pig','cat','olive baboon','sheep']
    Species.amendSpecies('Papio anubis', 'olive baboon')
    Species.amendSpecies('Ovis aries', 'sheep')
    Species.amendSpecies('Erinaceus europaeus', 'hedgehog')
    Species.amendSpecies('Mustela putorius furo', 'ferret')
except Exception as e:
    print (e)


def get_orthologs(refgenome, ensid=None, symbol=None):
    '''if ensid!=None:
        mygene = cow.getGeneByStableId(StableId=ensid)
    else:
        mygene = cow.getGenesMatching(symbol=symbol)'''
    #print mygene
    orthologs = comp.getRelatedGenes(gene_region=mygene,
                    Relationship='ortholog_one2one')
    print (orthologs.Members)
    #seqs = orthologs.getSeqCollection(feature_types='gene')
    #print seqs.Names
    #syntenicregions = comp.getSyntenicRegions(region=mygene,
    #                align_method='EPO', align_clade='eutherian')
    return

def get_genes_from_location(ref, coords, pad=0):
    """Get genes from a set of genome coordinates.
       pad will add n bases to either side to expand area"""

    genome = Genome(Species=ref, Release=release, account=account)
    chrom,start,end,strand = coords
    genes = list(genome.getFeatures(CoordName=chrom, Start=start-pad,
                    End=end+pad, feature_types='gene'))
    return genes

def find_in_gene(g, start, end):
    """Find if coords are inside intron or exon"""

    start=int(start)
    end=int(end)
    tr = g.CanonicalTranscript
    for i in tr.Exons:
        s,e = i.Location.Start,i.Location.End
        if start-5>=s and end<=e+5:
            return 'exon'
    if tr.Introns is None: return
    for i in tr.Introns:
        s,e = i.Location.Start,i.Location.End
        if start-5>=s and end<=e+5:
            return 'intron'
    #return 'both'

def get_alignment_tree(fname):
    """Build a neighbour joining tree"""

    from cogent.phylo import distance, nj
    from cogent.evolve.models import HKY85,F81
    al = cogent.LoadSeqs(fname,format='fasta')
    d = distance.EstimateDistances(al, submodel= F81())
    d.run(show_progress=False)
    mytree = nj.nj(d.getPairwiseDistances())
    mytree = mytree.balanced()
    print (mytree.asciiArt())
    print
    '''from cogent.draw import dendrogram
    p = dendrogram.SquareDendrogram(mytree)
    p.drawToPDF('tree-scaled.pdf', 500, 400, stroke_width=2.0,
                shade_param = 'r', max_value = 1.0,)'''
    return

def get_syntenic_alignment(comp, ref, coords, fname='ensembl.aln.fa'):
    """Get ensembl genomic alignment for a location using compara and return seqs"""

    print
    chrom,start,end,strand = coords
    regions = comp.getSyntenicRegions(Species=ref, CoordName=chrom,
                                       Start=start, End=end, align_method="EPO_LOW_COVERAGE",
                                       align_clade="39 eutherian", Strand=strand)
    regions = [r for r in regions]
    if len(regions)==0:
        print ('no alignments')
        return None, None
    #usually only 1 alignment
    A = regions[0].getAlignment()
    A.writeToFile(fname)
    #print cogent.LoadSeqs(fname)
    #if len(A.Seqs)>3:
    #    getAlignmentTree(fname)
    if regions is not None:
        print ('%s syntenic regions found' %len(regions))
    return regions, A

def get_host_genes(df, ref='cow'):
    """Get all genes containing the given miRNAs using ensembl"""

    results=[]
    comp = Compara(species, account=None, Release='79')
    for i, r in list(df.iterrows()):
        name = r['#miRNA']
        c,locs,strand = r['precursor coordinate'].split(':')
        start,end = locs.split('..')
        coords = c,int(start),int(end),strand
        genes = get_genes_from_location(ref, coords)
        for g in genes:
            if g.BioType != 'miRNA':
                tu = findinGene(g, start, end)
                results.append((name,g.Symbol,g.Location,g.BioType,tu,g.StableId))
                #print name,g.Symbol,tu
    if len(results)>0:
        results = pd.DataFrame(results,columns=['#miRNA','gene','location','biotype',
                                'tu','ensid'])
        return results
    return

#def getGeneProperty(genes, label):
#    x = [g[0].Symbol if len(g)>0 else np.nan for g in orthgenes]
#    return x

def get_mirna_orthologs(df, comp=None, ref='cow'):
    """Get all possible orthologs/conservation for miRNAs using ensembl"""

    if comp == None:
        comp = Compara(species, account=None, Release='79')
    results=[]
    for i, r in list(df.iterrows()):
        #base.RNAfold(r['consensus precursor sequence'], r['#miRNA']+'_'+ref)
        mature = r['consensus mature sequence'].replace('u','t').upper()
        seed = r['seed'].replace('u','t').upper()
        c,locs,strand = r['precursor coordinate'].split(':')
        start,end = locs.split('..')
        print (r['#miRNA'], seed, mature, locs, strand)
        coords = (c,int(start),int(end),strand)
        regions, aln = getSyntenicAlignment(comp, ref, coords, fname=r['#miRNA']+'.aln.fa')
        if aln == None:
            x=pd.DataFrame([r])
            a = getHostGenes(x)
            if a is not None:
                results.append(a)
            continue
        region = regions[0]
        a = base.cogentAlignment2DataFrame(aln.degap())
        a['#miRNA'] = r['#miRNA']
        a['ident'] = getIdentities(aln)
        print ('max identity: %s' %a.ident.max())
        a['seedcons'] = getSeqConservation(aln, seed)
        orthgenes = getGenesinRegion(region)

        a['gene'] = [g[0].Symbol if len(g)>0 else np.nan for g in orthgenes]
        a['gene_loc'] = [g[0].Location if len(g)>0 else np.nan for g in orthgenes]
        locs = getLocations(region)
        a['location'] = [':'.join(str(l).split(':')[2:]) for l in locs]
        a['biotype'] = [g[0].BioType if len(g)>0 else np.nan for g in orthgenes]
        a['ensid'] = [g[0].StableId if len(g)>0 else np.nan for g in orthgenes]
        #find where in gene the miRNA is, usually introns
        trpts = []
        for l,g in zip(locs,orthgenes):
            if len(g)==0:
                trpts.append(np.nan)
            else:
                tr = findinGene(g[0],l.Start,l.End)
                trpts.append(tr)
        a['tu'] = trpts

        #get RNAfold energy for each sequence
        a['energy'] = a.apply(lambda x : base.RNAfold(x.seq)[1],1)
        results.append(a)
        print (a)
        print ('--------------------------------------------------------')
    results = pd.concat(results).reset_index(drop=True)
    #results.drop(columns='seq')
    results.to_csv('novel_orthologs.csv')
    return

def get_locations(region):
    """Locations with region alignments"""

    l=[]
    for r in region.Members:
        try:
            loc = r.Location
            l.append(loc)
        except:
            continue
    return l

def get_seq_conservation(aln, seq):
    """Get position of a sub sequence and return position if it's
        conserved for each aligned species, -1 if not"""

    vals=[]
    for s in aln.Seqs:
        print (seq, s, str(s).find(seq))
        vals.append(str(s).find(seq))
    return vals

def get_identities(aln):
    """Identities for all seqs in a cogent alignment"""

    names=aln.Names
    ref=names[0]
    vals=[None]
    for n in names[1:]:
        new = aln.takeSeqs([ref, n])
        ident = round(len(new.filtered(lambda x: len(set(x)) == 1))/float(len(new)),3)
        vals.append(ident)
    return vals

def get_genes_in_region(region):
    """Get orthologous genes in all species in a syntenic region"""

    print ('checking genes in aligned species')
    orthologs=[]
    for r in region.Members:
        try:
            loc = r.Location
        except:
            continue
        sp = loc.Species
        coords = loc.CoordName,loc.Start,loc.End,loc.Strand
        #print sp, loc
        genes = list(r.genome.getFeatures(feature_types='gene',region=r))
        orthologs.append(genes)
    return orthologs

def get_ests(region):

    for r in region.Members:
        loc = r.Location
        sp = loc.Species
        ests = r.genome.getFeatures(feature_types='est', region=r)
        for est in ests:
            print (est)
    return

def get_sequence_from_location(species, coords):
    """Get sequence from a genomic location in an ensembl species genome."""

    from cogent.db.ensembl import HostAccount, Genome, Compara, Species
    genome = Genome(Species=species, Release='87', account=None)
    chrom,start,end,strand = coords
    #print coords
    r = genome.getRegion(CoordName=str(chrom), Start=start,End=end,Strand=strand)
    return r.Seq

def test():

    comp = Compara(species, account=account, Release=release)
    #print comp.method_species_links
    coords = [('15', 85273455, 85273507, '+'),('18',12423976,12424060,'+')]
    for c in coords:
        getSyntenicAlignment(comp, 'cow', c)
        #getAlignmentTree('ensembl_aln.fa')
    return
