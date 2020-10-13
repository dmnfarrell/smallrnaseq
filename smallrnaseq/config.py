#!/usr/bin/env python

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
Module for smallrnaseq configuration file. Used with command line app.
Created Jan 2017
Copyright (C) Damien Farrell
"""

from __future__ import absolute_import, print_function
import sys, os, string, time
import types, re, subprocess, glob, shutil
import pandas as pd
try:
    import configparser
except:
    import ConfigParser as configparser

path = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(path, 'data')
from . import aligners

baseoptions = {'base': [('filenames',''),('path',''),('overwrite',0),
                    ('adapter',''),
                    ('index_path','indexes'),
                    ('libraries',''),
                    ('ref_fasta',''),('features',''),
                    ('output','results'),('add_labels',0),
                    ('aligner','bowtie'),
                    ('mirna',0),('species','hsa'),('pad5',3),('pad3',5),
                    ('verbose', 1),
                    ('cpus',1)],
               'aligner': [('default_params','-v 1 --best'),
                    ('mirna_params',aligners.BOWTIE_MIRBASE_PARAMS)],
               'novel': [('score_cutoff',.7), ('read_cutoff',100),
                        ('strict',0)],
               'de': [('count_file',''),('sample_labels',''),('sep',','),
                    ('sample_col',''),('factors_col',''),
                    ('conditions',''),('logfc_cutoff',1.5),
                    ('de_plot','point')]
                    }

def write_default_config(conffile='default.conf', defaults={}):
    """Write a default config file"""

    if not os.path.exists(conffile):
        cp = create_config_parser_from_dict(defaults, ['base','novel','aligner','de'])
        cp.write(open(conffile,'w'))
        print ('wrote config file %s' %conffile)
    return conffile

def create_config_parser_from_dict(data, sections, **kwargs):
    """Helper method to create a ConfigParser from a dict and/or keywords"""

    cp = configparser.ConfigParser()
    for s in sections:
        cp.add_section(s)
        if not data.has_key(s):
            continue
        for i in data[s]:
            name,val = i
            cp.set(s, name, str(val))

    #use kwargs to create specific settings in the appropriate section
    for s in cp.sections():
        opts = cp.options(s)
        for k in kwargs:
            if k in opts:
                cp.set(s, k, kwargs[k])
    return cp

def parse_config(conffile=None):
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

def get_options(cp):
    """Makes sure boolean opts are parsed"""

    from collections import OrderedDict
    options = OrderedDict()
    #options = cp._sections['base']
    for section in cp.sections():
        options.update( (cp._sections[section]) )
    for o in options:
        for section in cp.sections():
            try:
                options[o] = cp.getboolean(section, o)
            except:
                pass
            try:
                options[o] = cp.getint(section, o)
            except:
                pass
    return options

def print_options(options):
    """Print option key/value pairs"""

    for key in options:
        print (key, ':', options[key])
    print ()

def check_options(opts):
    """Check for missing default options in dict. Meant to handle
       incomplete config files"""

    sections = baseoptions.keys()
    for s in sections:
        defaults = dict(baseoptions[s])
        for i in defaults:
            if i not in opts:
                opts[i] = defaults[i]
    return opts
