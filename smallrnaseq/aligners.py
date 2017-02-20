#!/usr/bin/env python

"""
    Methods for calling short read aligners
    Created Jan 2017
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
import sys, os, string, types, re
import shutil, glob, collections
import itertools
import subprocess
import numpy as np
import pandas as pd
from . import utils

BOWTIE_INDEXES = None
BOWTIE_PARAMS = None
SUBREAD_INDEXES = None
SUBREAD_PARAMS = '-m 2 -M 2'
BWA_INDEXES = None


def bwa_align(infile, ref=None, bowtie_index=None, outfile=None):
    """Align reads with bwa"""

    if bwa_index == None:
        bwa_index = BWA_INDEXES
    ref = os.path.join(bwaindexes, ref)
    label = os.path.splitext(os.path.basename(infile))[0]
    outfile = label+'_'+ref+'_bwa.sam'
    cmd1 = 'bwa aln -n 0 -t 2 %s %s > out.sai' %(ref,infile)
    cmd2 = 'bwa samse %s out.sai %s > %s' %(ref,infile,outfile)
    result = subprocess.check_output(cmd1, shell=True, executable='/bin/bash')
    result = subprocess.check_output(cmd2, shell=True, executable='/bin/bash')
    return

def build_bowtie_index(fastafile, path):
    """Build a bowtie index"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'bowtie-build -f %s %s' %(fastafile, name)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    files = glob.glob(name+'*.ebwt')
    utils.move_files(files, path)
    return

def build_subread_index(fastafile, path):
    """Build an index for subread"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'subread-buildindex -o %s %s' %(name,fastafile)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    exts = ['.00.b.array','.00.b.tab','.files','.reads']
    files = [name+i for i in exts]
    utils.move_files(files, path)
    return

def bowtie_align(infile, ref, outfile=None, remaining=None, verbose=True):
    """Map reads using bowtie"""

    label = os.path.splitext(os.path.basename(infile))[0]
    outpath = os.path.dirname(os.path.abspath(infile))
    if outfile == None:
        outfile = label+'_'+ref+'_bowtie.sam'

    if BOWTIE_INDEXES == None:
        print ('aligners.BOWTIE_INDEXES variable not set')
        return
    os.environ["BOWTIE_INDEXES"] = BOWTIE_INDEXES
    params = BOWTIE_PARAMS
    if remaining == None:
        remaining = os.path.join(outpath, label+'_r.fa')
    cmd = 'bowtie -f -p 2 -S %s --un %s %s %s > %s' %(params,remaining,ref,infile,outfile)
    if verbose == True:
        print (cmd)
    try:
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print (str(e.output))
    if verbose == True:
        print (result)
    return remaining

def subread_align(infile, ref, outfile):
    """Align reads with subread"""

    if SUBREAD_INDEXES == None:
        print ('aligners.SUBREAD_INDEXES variable not set')
        return
    ref = os.path.join(SUBREAD_INDEXES, ref)
    params = '-t 0 --SAMoutput -T 2 %s' %SUBREAD_PARAMS
    from subprocess import Popen, PIPE
    cmd = 'subread-align %s -i %s -r %s -o %s' %(params, ref, infile, outfile)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash',
                                     stderr= subprocess.STDOUT)
    return