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
from itertools import islice
import subprocess
import matplotlib
import pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
try:
    import HTSeq
except:
    'HTSeq not present'

def first(x):
    return x.iloc[0]

def move_files(files, path):
    if not os.path.exists(path):
        os.mkdir(path)
    for f in files:
        shutil.move(f, os.path.join(path,os.path.basename(f)))
    return

def remove_files(path, wildcard=''):
    """remove files from a folder"""

    files = glob.glob(os.path.join(path, wildcard))
    for f in files:
        os.remove(f)
    return

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

def make_blastdb(fastafile, title=None):
    """Make local blast db"""

    if title==None: title=os.path.splitext(fastafile)[0]
    cmd = 'makeblastdb -in %s -dbtype nucl -out %s' % (fastafile, title)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    return

def run_blastn(database, query, params='-e .1 -G 10'):
    """Run blast"""

    out = 'blast_result.csv'
    cmd = 'blastall -d %s -i %s -p blastn -m 8 %s > %s' %(database,query,params,out)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    #zipfile(out+'.xml', remove=True)
    return

def local_blast(fastafile, database, ident=100, params='-e 10 -a 2', results='all'):
    """Blast a blastdb and save hits to a dataframe.
       Args:
            fastafile: file with queries
            database: name of blast db
            ident: cutoff percentage identity
            params: custom blast parameters (see blastall -h)
            results: return 'all' or 'best'
        Returns: dataframe of hits"""

    outname = 'blast_result.csv'
    run_blastn(database, fastafile, params)
    cols = ['query', 'subj', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore']
    res = pd.read_csv(outname, names=cols, sep='\t')
    print ('found %s hits in db' %len(res))
    res = res[res['pident']>=ident]
    if results == 'best':
        res = res.groupby('query').first().reset_index()
    queries = fasta_to_dataframe(fastafile).reset_index()
    res = res.merge(queries, left_on='query', right_on='name', how='left')
    print ()
    return res

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

def dataframe_to_fasta(df, outfile='out.fa', seqkey='seq', idkey='id'):
    """Convert dataframe to fasta"""

    if idkey not in df.columns:
        df = df.reset_index()
    fastafile = open(outfile, "w")
    for i,row in df.iterrows():
        try:
            seq = row[seqkey]
        except:
            continue
        seq = seq.upper().replace('U','T')
        if idkey in row:
            d = str(row[idkey])
        else:
            d = row.name
        seq = seq.encode()
        myseq = HTSeq.Sequence(seq, str(d))
        myseq.write_to_fasta_file(fastafile)
    fastafile.close()
    return

def fasta_to_dataframe(infile, idindex=0):
    """Get fasta proteins into dataframe"""

    keys = ['name','sequence','description']
    fastafile = HTSeq.FastaReader(infile)
    data = [(s.name, s.seq.decode(), s.descr) for s in fastafile]
    df = pd.DataFrame(data,columns=(keys))
    df.set_index(['name'],inplace=True)
    return df

def fastq_to_dataframe(f, size=None):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    ext = os.path.splitext(f)[1]
    if ext=='.fastq':
        ffile = HTSeq.FastqReader(f, "solexa")
    elif ext == '.fa':
        ffile = HTSeq.FastaReader(f)
    else:
        return
    if size != None:
        sequences = [(s.name, s.seq, s.descr) for s in islice(fastfile, i, i+size)]
    else:
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

def trim_adapters(infile, adapter, outfile='cut.fastq', method='default'):
    """Trim adapters using cutadapt"""

    if not type(adapter) is str:
        print ('not valid adapter')
        return

    if method == 'default':
        newfile = open( outfile, "w" )
        fastfile = HTSeq.FastqReader(infile, "solexa")
        a = HTSeq.Sequence(adapter)
        for s in fastfile:
            new = s.trim_right_end(a, mismatch_prop = 0.)
            new.write_to_fastq_file( newfile )
        newfile.close()
    elif method == 'cutadapt':
        cmd = 'cutadapt -m 18 -O 5 -q 20 -a %s %s -o %s' %(adapter,infile,outfile)
        print (cmd)
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
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
    try:
        x = RNA.fold(seq)
    except Exception as e:
        print (e)
        return
    return x

def rnaplot(seq, struct=None, path='rnaplots', name='temp'):

    import RNA
    if struct==None:
        struct = RNA.fold(seq)[0]
    filename = os.path.join(path,name+'.ps')
    #RNA.svg_rna_plot(seq,struct,filename)
    colors = [" 1. 0. .2", " 0. .9 .5"]
    macro = format_cmark_values(range(0,10), rgb=colors[0])
    RNA.PS_rna_plot_a(seq, struct, filename, '', macro)
    return filename

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

def bed_to_dataframe(bedfile):
    """Bed file to dataframe"""

    header=['chrom','chromStart','chromEnd','name','score','strand','thickStart',
            'thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    feats = pd.read_csv(bedfile, sep='\t', names=header)
    #feats['chr'] = feats.chrom.str.extract('(\d+)')
    feats['chr'] = feats.chrom.str[3:]
    return feats

def features_to_gtf(df, filename):
    """Take generic dataframe of features and create ensembl gtf file. Note some fields
       will be redundant as they require ensembl specific information"""

    #ensembl gtf header format
    gtfheader=['chrom', 'start', 'end', 'exon_id', 'exon_number', 'exon_version', 'gene_biotype', 'gene_id',
           'gene_name', u'gene_source', u'gene_version', 'id', 'protein_id', 'protein_version',
           'strand', 'transcript_biotype', 'transcript_id', 'transcript_name',
           'transcript_source', 'transcript_version']
    rows=[]
    for i,r in df.iterrows():
        #print r
        row = [r.chr,r.chromStart+1,r.chromEnd,r['name'],1,1,'tRNA',r['name'],r['name'],
               'gtrnadb',1,'','','',r.strand,'tRNA',r['name'],'','gtrnadb',1]
        rows.append(row)
    gtf = pd.DataFrame(rows,columns=gtfheader)

    f=open(filename,'w')
    #f.write('#custom gtf file\n')
    for idx,r in gtf.iterrows():
        c1 = ['chrom', 'start', 'end']
        s1 = '\t'.join([str(r.chrom), 'gtrnadb','exon', str(r.start), str(r.end)])
        s2 = '\t'.join(['.',r.strand,'.'])
        c2 = ['gene_id','gene_version','transcript_id','transcript_version','exon_number',
              'gene_source','gene_biotype','transcript_source','transcript_biotype',
              'exon_id','exon_version']
        s3 = '; '.join(i[0]+' '+'"%s"' %str(i[1]) for i in zip(c2,r[c2]))
        s = '\t'.join([s1,s2,s3])
        f.write(s); f.write('\n')
    return gtf

def sequence_from_coords(fastafile, coords):
    """Fasta sequence from genome coords"""

    from pybedtools import BedTool
    chrom,start,end,strand = coords
    if not os.path.exists(fastafile):
        print ('no such file')
        return
    try:
        if strand == '+':
            seq = str(BedTool.seq(coords, fastafile))
        else: #reverse strand
            seq = str(BedTool.seq(coords, fastafile))
            seq = str(HTSeq.Sequence(seq).get_reverse_complement())
    except Exception as e:
        print (e)
        return
    return seq

def sequence_from_bedfile(fastafile, features=None, bedfile=None, pad5=0, pad3=0):
    """Fasta sequences from set of genomic features in a bed file
        Args:
            fastafile: fasta file with genomic sequence
            features: dataframe of features/coords with bed file col names
            bedfile: optionally provide a bed file instead
            pad5,pad3: flanking sequence at 5' or 3' ends
        Returns:
            a pandas dataframe with name, sequence and coord columns"""

    from pybedtools import BedTool
    if bedfile != None:
        features = utils.bed_to_dataframe(bedfile)
    new = []
    for n,r in features.iterrows():
        if r.strand == '+':
            coords = (r.chr,r.chromStart-pad5,r.chromEnd+pad3)
            seq = str(BedTool.seq(coords, fastafile))
        else: #reverse strand
            coords = (r.chr,r.chromStart-pad3,r.chromEnd+pad5)
            seq = str(BedTool.seq(coords, fastafile))
            seq = HTSeq.Sequence(seq).get_reverse_complement()
        #print n, coords, r['name']
        new.append([r['name'],str(seq),coords])
    new = pd.DataFrame(new, columns=['name','seq','coords'])
    return new

def get_csv_files(path, filename, names, **kwargs):
    """Get multiple csv files in set of sub folders, used for re-loading previously
       created results for different datasets into one pandas dataframe"""

    res = []
    for n in names:
        p = os.path.join(path, n)
        f = os.path.join(p,filename)
        if not os.path.exists(f):
            continue
        #print f
        data = pd.read_csv(f,**kwargs)
        data['study'] = n
        #r = studies[studies.label==n].iloc[0]
        #data['source'] = r['source']
        #data['type'] = r['type']
        res.append(data)
    return pd.concat(res)

def read_collapsed_file(collapsed):
    """Read seqs/counts from a collapsed file using values in the
       fasta id to get original counts"""

    #original read counts are encoded in fasta names
    df = fasta_to_dataframe(collapsed).reset_index()
    #original count stored in id as second value
    df['read_id'], df['reads'] = df.name.str.split('_', 1).str
    df['reads'] = df.reads.astype(int)
    df['read_id'] = df.read_id.astype(int)
    df = df.drop(['name','description'],1)
    df = df.rename(columns={'sequence':'seq'})
    #print (df[:5])
    return df

def get_aligned_reads(samfile, collapsed=None, readcounts=None):
    """Get all aligned reads from a sam file into a pandas dataframe"""

    sam = HTSeq.SAM_Reader(samfile)
    f=[]
    for a in sam:
        if a.aligned == True:
            seq = a.read.seq.decode()
            f.append((seq,a.read.name,a.iv.chrom,a.iv.start,a.iv.end,a.iv.strand))
        #else:
        #    f.append((seq,a.read.name,'_unmapped'))
    counts = pd.DataFrame(f, columns=['seq','read','name','start','end','strand'])
    counts['length'] = counts.seq.str.len()
    counts = counts.drop(['read'],1)
    if collapsed is not None:
        readcounts = read_collapsed_file(collapsed)
    if readcounts is not None:
        counts = counts.merge(readcounts, on='seq')
        counts['align_id'] = counts.index
    return counts

def combine_aligned_reads(path, idx, filenames=None):
    """Combine reads from mapping of multiple samples to a genome/library.
    Args:
        path: folder with sam files and collapsed counts from mapping, should
        contain results corresponding to filenames
        idx: name of index mapped to
        filenames: names of files
    Returns: pandas dataframe with all read counts summed
    """

    a = []
    #files are collapsed files
    if filenames == None:
        filenames = glob.glob(os.path.join(path, '*[!_r].fa'))
    #print (filenames)
    for f in filenames:
        #print(f)
        name = os.path.splitext(os.path.basename(f))[0]
        samfile = os.path.join(path, '%s_%s.sam' %(name,idx))
        if not os.path.exists(samfile):
            print ('no sam file')
            continue
        reads = get_aligned_reads(samfile, collapsed=f)
        a.append(reads)
    if len(a) == 0:
        return
    a = pd.concat(a)
    cols = ['seq', 'name', 'start', 'end', 'strand', 'length']
    s = a.groupby(cols).agg({'reads':np.sum}).reset_index()
    s['read_id'] = s.index.copy()
    s = s.sort_values(by='reads',ascending=False)
    print ('pooled %s files into %s unique reads' %(len(filenames),len(s)))
    return s

def print_read_stacks(reads, fastafile, outfile=None, name=None, by=None):
    """Print multiple read alignments to file or stdout
       Args:
        reads: dataframe of read counts with position info
        fastafile: optional fasta file with references sequences mapped to
    """

    if name != None:
        names = [name]
    else:
        names = reads.name.unique()
    if outfile != None:
        f = open(outfile, 'w')

    refs = fasta_to_dataframe(fastafile)
    s = ''
    for n in names:
        x = reads[reads.name==n]
        refseq = refs.ix[n].sequence
        s += print_read_stack(x, refseq, by=by, cutoff=0, label=n)
    #if f != None:
    #    f.close()
    return s

def print_read_stack(reads, refseq=None, cutoff=0, by=None, label=''):
    """Print local read alignments from a sam file against the mapped sequence.
       Args:
        reads: dataframe of read counts with position info
        refseq: sequence segment or gene that reads are mapped to
        cutoff: don't display read with <=cutoff counts
        by: column to sort reads by, e.g. 'reads', 'start' or 'end'
       Returns:
        string of printed read stacks to display or write to a file
    """

    if refseq != None:
        seqlen = len(refseq)
    else:
        seqlen = reads.end.max()
    f = None
    reads = reads[reads.reads>=cutoff]
    if by is not None:
        reads = reads.sort_values(by, ascending=False)

    l = len(reads)
    if l==0: return
    s = ''
    s += '%s unique reads\n' %l
    s += '-------------------------------------\n'
    s += refseq+'\n'
    for idx,r in reads.iterrows():
        seq = r.seq
        count = r.reads
        pos = find_subseq(refseq, seq)
        if pos == -1: continue
        i = len(seq)+pos
        s += "{:>{w}} ({c})\n".format(seq,w=i,c=count)
    s += '\n'
    return s

def plot_read_stack(reads, refseq=None, by=None, cutoff=0, ax=None):
    """Plot read stack using coverage at each position"""

    reads = reads.copy()
    if by != None:
        reads = reads.sort_values(by)
    else:
        reads = reads.sort_values(by=['start','reads'],ascending=[1,0])
    if refseq != None:
        seqlen = len(refseq)
    else:
        offset = reads.start.min()
        reads.start = reads.start-offset
        reads.end = reads.end-offset
        seqlen = reads.end.max()
    reads = reads[reads.reads>cutoff]
    if len(reads)==0:
        return

    def pos_coverage(r, p):
        x = [r.reads if (i>=r.start and i<=r.end) else 0 for i in p]
        return pd.Series(x,index=p)

    reads = reads.set_index('reads',drop=False)
    p = range(1,seqlen+1)
    m = reads.apply( lambda x: pos_coverage(x,p), 1 )
    m = m.replace(0,1)

    from matplotlib.colors import LogNorm
    if ax == None:
        h = 12*len(m)/80+1
        fig,ax = plt.subplots(1,1,figsize=(12,h))
    ax = sns.heatmap(m,ax=ax,norm=LogNorm(vmin=m.min(), vmax=m.max()),
                cmap='Blues',cbar=False)
    plt.gca().xaxis.grid(True)
    start, end = ax.get_xlim()
    xlbls = np.arange(start, end, 5).astype(int)
    ax.xaxis.set_ticks(xlbls)
    ax.set_xticklabels(xlbls)
    ax.tick_params(axis='y', which='major', labelsize=9)
    for xmin in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmin,ls='--',lw=0.5,color='gray')
    #m.plot(kind='bar',width=.8,figsize=(12,2))
    return ax

def close_match(ref, seq, mm=1):
    """Find string inside another with single mismatch"""

    sr = re.compile('|'.join(seq[:i]+'.'+seq[i+1:] for i in range(len(seq))))
    m = sr.findall(ref)
    if len(m)==0: return -1
    idx = ref.find(m[0])
    return idx

def find_subseq(seq, s):
    """Find a sub sequence inside another
    Returns:
        index inside seq or -1 if not found
    """

    for i in range(16,4,-4):
        c = seq.find(s[:i])
        if c != -1: return c
    #if above fails try to find single mismatches
    return close_match(seq, s)

def featurecounts(samfile, gtffile):
    """Count aligned features with the featureCounts program.
       Returns: a dataframe of counts"""

    params = '-T 5 -t exon -g transcript_id'
    cmd = 'featureCounts %s -a %s -o counts.txt %s' %(params, gtffile, samfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    counts =  pd.read_csv('counts.txt', sep='\t', comment='#')
    counts = counts.rename(columns={samfile:'reads'})
    counts = counts.sort('reads', ascending=False)
    return counts

def get_bg(seq, struct=None):
    """Get bulgegraph object from sequence/struct, requires
    the forgi package to be installed
    Returns:
        a BulgeGraph object
    """

    import forgi.graph.bulge_graph as cgb
    if struct == None:
        struct,sc = utils.rnafold(seq)
    bg = cgb.BulgeGraph()
    bg.from_dotbracket(struct)
    bg.seq = seq
    return bg