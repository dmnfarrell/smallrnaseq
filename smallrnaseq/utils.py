#!/usr/bin/env python

"""
    Misc miRNA analysis routines
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
try:
    import HTSeq
except:
    'HTSeq not present'

def gzipfile(filename, remove=False):
    """Compress a file with gzip"""

    import gzip
    fin = open(filename, 'rb')
    fout = gzip.open(filename+'.gz', 'wb')
    fout.writelines(fin)
    fout.close()
    fin.close()
    if remove == True:
        os.remove(filename)
    return

def create_html(df,name,path='.'):
    """Create a basic html page for dataframe results"""

    s = ['<script src="sorttable.js"></script>']
    s.append('<link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.5.0/pure-min.css">')
    s.append('<body><h2>'+name+'</h2><div class="pure-div">')
    table = df.to_html(classes=['sortable','pure-table-striped'])
    s.append(table)
    body = '\n'.join(s)
    f = open(os.path.join(path,name)+'.html','w')
    f.write(body)
    return

def run_blastn(database, query):
    """Run blast"""

    out = os.path.splitext(query)[0]
    cmd = 'blastall -d %s -i %s -p blastn -m 7 -e .1 > %s.xml' %(database,query,out)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    gzipfile(out+'.xml', remove=True)
    return

def parse_blast_rec(rec):
    """Parse blast record alignment(s)"""

    #if len(rec.alignments) == 0 : print 'no alignments'
    recs=[]
    qry = rec.query.split()[0]
    for align in rec.alignments:
        hsp = align.hsps[0]
        subj = align.title.split()[1]
        if qry == subj: continue
        recs.append([qry, subj, hsp.score, hsp.expect, hsp.identities,
                    hsp.positives, hsp.align_length])
    return recs

def get_blast_results(handle=None, filename=None, n=80):
    """Get blast results into dataframe"""

    from Bio.Blast import NCBIXML
    import gzip
    if filename!=None:
        #handle = open(filename)
        handle = gzip.open(filename, 'rb')
    blastrecs = NCBIXML.parse(handle)
    rows=[]
    for rec in blastrecs:
        r = parseBlastRec(rec)
        rows.extend(r)
        #print r
    df = pd.DataFrame(rows, columns=['query','subj','score','expect','identity',
                            'positive','align_length'])
    df['perc_ident'] = df.identity/df.align_length*100
    return df

def blastDB(f, database, ident=100):
    """Blast a blastdb and save hits to csv"""

    outname = os.path.splitext(f)[0]
    runBlastN(database, f)
    df = getBlastResults(filename=outname+'.xml.gz')
    df = df[df['perc_ident']>=ident]
    #print df[:10]
    g = df.groupby('query').agg({'subj':first})
    g = g.sort('subj',ascending=False)
    g = g.reset_index()
    #print g[:15]
    print ('found %s hits in db' %len(df))
    print ()
    #outname = os.path.splitext(f)[0]+'_hits.csv'
    #g.to_csv(outname)
    return g

def fastq_to_fasta(infile, rename=True):
    """Fastq to fasta"""

    fastqfile = HTSeq.FastqReader(infile, "solexa")
    outfile = open(os.path.splitext(infile)[0]+'.fa','w')
    i=1
    for s in fastqfile:
        if rename==True:
            s.name=str(i)
        s.write_to_fasta_file(outfile)
        i+=1
    outfile.close()
    return

def dataframe_to_fasta(df, seqkey='seq', idkey='id', outfile='out.fa'):
    """Convert dataframe to fasta"""

    df = df.reset_index() #in case key is the index
    fastafile = open(outfile, "w")
    for i,row in df.iterrows():
        seq = row[seqkey].upper().replace('U','T')
        if idkey in row:
            d = row[idkey]
        else:
            d = ''
        myseq = HTSeq.Sequence(seq, d)
        myseq.write_to_fasta_file(fastafile)
    return

def fasta_to_dataframe(infile,idindex=0):
    """Get fasta proteins into dataframe"""

    keys = ['name','sequence','description']
    fastafile = HTSeq.FastaReader(infile)
    data = [(s.name, s.seq, s.descr) for s in fastafile]
    df = pd.DataFrame(data,columns=(keys))
    df.set_index(['name'],inplace=True)
    return df

def fastq_to_dataframe(f):
    """Convert fastq to dataframe - may use a lot of memory"""

    ext = os.path.splitext(f)[1]
    if ext=='.fastq':
        ffile = HTSeq.FastqReader(f, "solexa")
    elif ext == '.fa':
        ffile = HTSeq.FastaReader(f)
    else:
        return
    sequences = [(s.name,s.seq) for s in ffile]
    df = pd.DataFrame(sequences,columns=['id','seq'])
    return df

def get_subset_fasta(infile, labels=['bta'], outfile='found.fa'):
    """Get a subset of sequences matching a label"""

    fastafile = HTSeq.FastaReader(infile)
    sequences = [(s.name, s.seq, s.descr) for s in fastafile]
    #print sequences[0][2]
    df = pd.DataFrame(sequences, columns=['id','seq','descr'])
    found=[]
    for l in labels:
        f = df[df.id.str.contains(l) | df.descr.str.contains(l)]
        found.append(f)
    df = pd.concat(found)
    print ('found %s sequences' %len(df))
    dataframe_to_fasta(df,outfile=outfile)
    return

def filter_fasta(infile):

    fastafile = HTSeq.FastaReader(infile)
    sequences = [(s.name, s.seq, s.descr) for s in fastafile]
    out = open('filtered.fa', "w")
    for s in sequences:
        if s[1] == 'Sequence unavailable':
            continue
        myseq = HTSeq.Sequence(s[1], s[0])
        myseq.write_to_fasta_file(out)
    return

def create_random_subset(sourcefile=None, sequences=None, size=1e5,
                        outfile='subset.fa'):
    """Generate random subset of reads"""

    if sequences==None:
        fastqfile = HTSeq.FastqReader(sourcefile, "solexa")
        sequences = [s.seq for s in fastqfile]
    randidx = np.random.randint(1,len(sequences),size)
    ffile = open(outfile, "w")
    for r in randidx:
        sequences[r].name = str(r)
        sequences[r].write_to_fasta_file(ffile)
    print ('wrote %s sequences to %s' %(size, outfile))
    return

def create_random_fastq(sourcefile, path, sizes=None):
    """Generate multiple random subsets of reads for testing"""

    fastqfile = HTSeq.FastqReader(sourcefile, "solexa")
    sequences = [s for s in fastqfile]
    print ('source file has %s seqs' %len(sequences))
    if sizes==None:
        sizes = np.arange(5e5,7.e6,5e5)
    for s in sizes:
        label = str(s/1e6)
        name = os.path.join(path,'test_%s.fa' %label)
        create_random_subset(sequences=sequences, size=s, outfile=name)
    return

def get_mifam():
    """Get miRBase family data"""

    cr=list(csv.reader(open('miFam.csv','r')))
    data=[]
    i=0
    for row in cr:
        if row[0]=='ID':
            fam=row[1]
        elif row[0]=='AC' or row[0]=='//':
            continue
        else:
            data.append((row[1],row[2],fam))
        i+=1
    df = pd.DataFrame(data,columns=['id','name','family'])
    return df

def trim_adapters(infile, adapters=[], outfile='cut.fastq'):
    """Trim adapters using cutadapt"""

    #if os.path.exists(outfile):
    #    return
    if len(adapters) == 0:
        print ('no adapters!')
        return
    adptstr = ' -a '.join(adapters)
    cmd = 'cutadapt -m 18 -O 5 -q 20 --discard-untrimmed -a %s %s -o %s' %(adptstr,infile,outfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    #print result
    return

def runEdgeR(countsfile, cutoff=1.5):
    """Run edgeR from R script"""

    cmd = 'Rscript ~/python/sandbox/mirnaseq/DEanalysis.R %s' %countsfile
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print (result)
    #read result back in
    de = pd.read_csv('de_output.csv')
    de.rename(columns={'Unnamed: 0':'name'}, inplace=True)
    de = de[(de.FDR<0.05) & ((de.logFC>cutoff) | (de.logFC<-cutoff))]
    return de

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

def cogentalignment_to_dataframe(A):
    """Pycogent alignment to dataframe"""

    res=[]
    for s in zip(A.Names,A.Seqs):
        res.append((s[0].split(':')[0],str(s[1])))
    df = pd.DataFrame(res,columns=['species','seq'])
    return df

def rnafold(seq, name=None):
    """Run RNAfold for precursor"""

    import RNA
    x = RNA.fold(seq)
    colors = [" 1. 0. .2", " 0. .9 .5"]
    if name != None:
        path=''
        RNA.svg_rna_plot(seq,x[0],os.path.join(path,name+'.svg'))
        #macro = format_cmark_values(range(0,10), rgb=colors[0])
        #RNA.PS_rna_plot_a(seq, x[0], name+'.ps', '', macro)
    return x

def format_cmark_values(values, rgb=" 1. 0. .2"):
    """PS colored marks for rnaplot"""

    minval , maxval = min ( values ) ,max ( values )
    valtab = [" %s %s cfmark"%(i,rgb) for i in values]
    #valtab = ["%s cmark " %i for i in values]
    x = "". join (valtab)
    macro = "/cfmark {setrgbcolor newpath 1 sub coor exch get aload"
    macro += " pop fsize 2 div 0 360 arc fill} bind def"+x
    return macro

def plot_rna_structure(seq, path='', subseqs=[], name='test'):
    """plot RNA structure using Vienna package"""

    import cogent.app.vienna_package as vienna
    colors = [" 1. 0. .2", " 0. .9 .5"]
    seq,struct,e = vienna.get_secondary_structure(seq)
    seqname='test'
    rp = vienna.RNAplot()
    i=0
    x=''
    if len(subseqs) > 0:
        for s in subseqs:
            ind = seq.find(s)+1
            e = ind+len(s)
            x += format_cmark_values(range(ind,e), rgb=colors[i])
            i+=1
        rp.Parameters['--pre'].on('"%s"' %x)
    rp(['>'+seqname,seq,struct])
    filename = os.path.join(path,'%s.png' %name)
    os.system('convert test_ss.ps %s' %filename)
    return filename

def muscle_alignment(filename=None, seqs=None):
    """Align 2 sequences with muscle"""

    if filename == None:
        filename = 'temp.faa'
        SeqIO.write(seqs, filename, "fasta")
    name = os.path.splitext(filename)[0]
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=filename, out=name+'.txt')
    stdout, stderr = cline()
    align = AlignIO.read(name+'.txt', 'fasta')
    return align

def sam_to_bam(filename):
    """Convert sam to bam"""

    import pysam
    infile = pysam.AlignmentFile(filename, "r")
    name = os.path.splitext(filename)[0]+'.bam'
    bamfile = pysam.AlignmentFile(name, "wb", template=infile)
    for s in infile:
        bamfile.write(s)
    pysam.sort("-o", name, name)
    pysam.index(name)
    bamfile = pysam.AlignmentFile(name, "rb")
    return