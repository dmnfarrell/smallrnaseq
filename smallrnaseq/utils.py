#!/usr/bin/env python

"""Misc miRNA analysis routines
   Created June 2014
   Copyright (C) Damien Farrell
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
from . import base

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

def rnafold(seq, name=None):
    """Run RNAfold for precursor"""

    import RNA
    x = RNA.fold(seq)
    print (name, x)
    if name != None:
        path='RNAplots'
        RNA.svg_rna_plot(seq,x[0],os.path.join(path,name+'.svg'))
    return x

def cogentalignment_to_dataframe(A):
    """Pycogent alignment to dataframe"""

    res=[]
    for s in zip(A.Names,A.Seqs):
        res.append((s[0].split(':')[0],str(s[1])))
    df = pd.DataFrame(res,columns=['species','seq'])
    return df

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
    """plot RNA structure using vienna package"""

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
            x += format_cmark_values(range(ind,e), rgb=colors[i])
            i+=1
        r.Parameters['--pre'].on('"%s"' %x)
    r(['>'+seqname,seq,struct])
    filename = os.path.join(path,'%s.ps' %name)
    os.system('convert test_ss.ps %s' %filename)
    #os.system('cp test_ss.ps %s' %filename)
    return filename
