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
BOWTIE_PARAMS = '-v 1 --best'
SUBREAD_INDEXES = None
SUBREAD_PARAMS = '-m 2 -M 1'

def get_current_params(aligner):
    if aligner == 'bowtie':
        global BOWTIE_PARAMS
        return BOWTIE_PARAMS

def set_params(aligner, params=None):
    """set aligner parameters"""
    if aligner == 'bowtie':
        global BOWTIE_PARAMS
        BOWTIE_PARAMS = params
    elif aligner == 'subread':
        global SUBREAD_PARAMS
        SUBREAD_PARAMS = params

def build_bowtie_index(fastafile, path):
    """Build a bowtie index"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'bowtie-build -f %s %s' %(fastafile, name)
    try:
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print (str(e.output))
        return
    files = glob.glob(name+'*.ebwt')
    utils.move_files(files, path)
    return

def build_subread_index(fastafile, path):
    """Build an index for subread"""

    name = os.path.splitext(fastafile)[0]
    cmd = 'subread-buildindex -o %s %s' %(name,fastafile)
    try:
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print (str(e.output))
        return
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
        result = subprocess.check_output(cmd, shell=True, executable='/bin/bash',
                                         stderr= subprocess.STDOUT)
        if verbose == True:
            print (result.decode())
    except subprocess.CalledProcessError as e:
        print (str(e.output))
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