#!/usr/bin/env python

"""Module for core utilities
   Created July 2014
   Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, string, time
import types, re, subprocess
import pylab as plt
import numpy as np
import pandas as pd
import configparser
try:
    import HTSeq
except:
    'HTSeq not present'

def writeDefaultConfig(conffile='default.conf', defaults={}):
    """Write a default config file"""

    if not os.path.exists(conffile):
        cp = createConfigParserfromDict(defaults, ['base'])
        cp.write(open(conffile,'w'))
        print ('wrote config file %s' %conffile)
    return conffile

def createConfigParserfromDict(data, sections, **kwargs):
    """Helper method to create a ConfigParser from a dict and/or keywords"""

    cp = configparser.ConfigParser()
    for s in sections:
        cp.add_section(s)
        if not data.has_key(s):
            continue
        for i in data[s]:
            name,val = i
            cp.set(s, name, val)
    #use kwargs to create specific settings in the appropriate section
    for s in cp.sections():
        opts = cp.options(s)
        for k in kwargs:
            if k in opts:
                cp.set(s, k, kwargs[k])
    return cp

def parseConfig(conffile=None):
    """Parse a configparser file"""

    f = open(conffile,'r')
    cp = configparser.ConfigParser()
    try:
        cp.read(conffile)
    except Exception as e:
        print ('failed to read config file! check format')
        print ('Error returned:', e)
        return
    f.close()    
    return cp

def getOptions(cp):
    """Makes sure boolean opts are parsed"""

    options = cp._sections['base']
    for o in options:
        try:
            options[o] = cp.getboolean('base', o)
        except:
            pass
    return options

def seabornsetup():
    global sns
    import seaborn as sns
    sns.set_style("ticks", {'axes.grid': False,'legend.frameon':True})
    sns.set_context("paper", rc={'axes.labelsize':16,'axes.titlesize':15,
                    'lines.color':1.0,'xtick.labelsize':12,
                    'ytick.labelsize': 12, 'legend.fontsize':12, 'title.fontsize':14,
                    'figure.figsize': np.array([8, 8])})
    return

def first(x):
    return x.iloc[0]

def doHeatMap(df,fname=None,cmap='seismic',log=False):
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

def venndiagram(names,labels,ax=None,**kwargs):
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

def createHtml(df,name,path='.'):
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

def getSubsetFasta(infile, labels=['bta'], outfile='found.fa'):
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
    dataframe2Fasta(df,outfile=outfile)
    return

def filterFasta(infile):
    fastafile = HTSeq.FastaReader(infile)
    sequences = [(s.name, s.seq, s.descr) for s in fastafile]
    out = open('filtered.fa', "w")
    for s in sequences:
        if s[1] == 'Sequence unavailable':
            continue
        myseq = HTSeq.Sequence(s[1], s[0])
        myseq.write_to_fasta_file(out)
    return

def fasta2DataFrame(infile,idindex=0):
    """Get fasta proteins into dataframe"""

    keys = ['name','sequence','description']
    fastafile = HTSeq.FastaReader(infile)
    data = [(s.name, s.seq, s.descr) for s in fastafile]
    df = pd.DataFrame(data,columns=(keys))
    df.set_index(['name'],inplace=True)
    return df

def dataframe2Fasta(df, seqkey='seq', idkey='id', outfile='out.fa'):
    """Convert dataframe to fasta"""

    df = df.reset_index() #in case key is the index
    fastafile = open(outfile, "w")
    for i,row in df.iterrows():
        seq = row[seqkey].upper().replace('U','T')
        myseq = HTSeq.Sequence(seq, row[idkey])
        #if 'descr' in df.columns:
        #    myseq.descr = str(row.descr)
        myseq.write_to_fasta_file(fastafile)
    return

def runBlastN(database, query):
    """Run blast"""

    out = os.path.splitext(query)[0]
    cmd = 'blastall -d %s -i %s -p blastn -m 7 -e .1 > %s.xml' %(database,query,out)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    gzipfile(out+'.xml', remove=True)
    return

def parseBlastRec(rec):
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

def getBlastResults(handle=None, filename=None, n=80):
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

def bwaMap(infile, ref=None, outfile=None):
    """Map with bwa"""

    bwaindexes = '/opt/mirnaseq/genomes/bwa_index'
    ref = os.path.join(bwaindexes, ref)
    label = os.path.splitext(os.path.basename(infile))[0]
    outfile = label+'_'+ref+'_bwa.sam'
    cmd1 = 'bwa aln -n 0 -t 2 %s %s > out.sai' %(ref,infile)
    cmd2 = 'bwa samse %s out.sai %s > %s' %(ref,infile,outfile)
    result = subprocess.check_output(cmd1, shell=True, executable='/bin/bash')
    result = subprocess.check_output(cmd2, shell=True, executable='/bin/bash')
    return

def bowtieMap(infile, ref, outfile=None, bowtieindex=None, params='-v 0 --best',
                remaining=None):
    """Map reads using bowtie"""

    label = os.path.splitext(os.path.basename(infile))[0]
    outpath = os.path.dirname(os.path.abspath(infile))
    if outfile == None:
        outfile = label+'_'+ref+'_bowtie.sam'
    if bowtieindex == None:
        bowtieindex = '/opt/mirnaseq/genomes/bowtie_index'
    os.environ["BOWTIE_INDEXES"] = bowtieindex

    if remaining == None:
        remaining = os.path.join(outpath, label+'_r.fastq')
    cmd = 'bowtie -f -p 2 -S %s --un %s %s %s > %s' %(params,remaining,ref,infile,outfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    print (result)
    return remaining

def createRandomSubset(sourcefile=None, sequences=None, size=1e5,
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

def createRandomFastqFiles(sourcefile, path, sizes=None):
    """Generate multiple random subsets of reads for testing"""

    fastqfile = HTSeq.FastqReader(sourcefile, "solexa")
    sequences = [s for s in fastqfile]
    print ('source file has %s seqs' %len(sequences))
    if sizes==None:
        sizes = np.arange(5e5,7.e6,5e5)
    for s in sizes:
        label = str(s/1e6)
        name = os.path.join(path,'test_%s.fa' %label)
        createRandomSubset(sequences=sequences, size=s, outfile=name)
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

def RNAfold(seq, name=None):
    """Run RNAfold for precursor"""

    import RNA
    x = RNA.fold(seq)
    print (name, x)
    if name != None:
        path='RNAplots'
        RNA.svg_rna_plot(seq,x[0],os.path.join(path,name+'.svg'))
    return x

def cogentAlignment2DataFrame(A):
    """Pycogent alignment to dataframe"""

    res=[]
    for s in zip(A.Names,A.Seqs):
        res.append((s[0].split(':')[0],str(s[1])))
    df = pd.DataFrame(res,columns=['species','seq'])
    return df

def formatcmarkValues(values, rgb=" 1. 0. .2"):
    """PS colored marks for rnaplot"""

    minval , maxval = min ( values ) ,max ( values )
    valtab = [" %s %s cfmark"%(i,rgb) for i in values]
    #valtab = ["%s cmark " %i for i in values]
    x = "". join (valtab)
    macro = "/cfmark {setrgbcolor newpath 1 sub coor exch get aload"
    macro += " pop fsize 2 div 0 360 arc fill} bind def"+x
    return macro

def plotRNA(seq, path='', subseqs=[], name='test'):
    """plot miRNA"""

    import cogent.app.vienna_package as vienna
    colors = [" 1. 0. .2", " 0. .9 .5"]
    seq,struct,e = vienna.get_secondary_structure(seq)
    seqname='test'
    r=vienna.RNAplot()
    i=0
    x=''
    if len(subseqs) > 0:
        for s in subseqs:
            ind=seq.find(s)+1
            e=ind+len(s)
            x += formatcmarkValues(range(ind,e), rgb=colors[i])
            i+=1
        r.Parameters['--pre'].on('"%s"' %x)
    r(['>'+seqname,seq,struct])
    filename = os.path.join(path,'%s.ps' %name)
    os.system('convert test_ss.ps %s' %filename)
    #os.system('cp test_ss.ps %s' %filename)
    return filename
