#!/usr/bin/env python

"""
    smallrnaseq entry point for command line functions
    Created November 2016
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
from smallrnaseq import base, analysis, utils

def main():
    """Run the application from outside the module - used for
       deploying as frozen app"""

    import sys, os
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config",
                        help="config file to use", metavar="FILE")
    parser.add_option("-r", "--run", dest="run",  action="store_true",
                        default=False, help="run smallrna mapping")
    parser.add_option("-t", "--test", dest="test",  action="store_true",
                        default=False, help="run tests")

    opts, remainder = parser.parse_args()
    if opts.test == True:
        pass
    elif opts.run == True:
        if opts.config == None:
            print ('No config file provided.')
            base.write_default_config('default.conf', defaults=base.baseoptions)
            return
        cp = base.parse_config(opts.config)
        options = base.get_options(cp)
        print (options)
        #base.map_rnas(**options)

if __name__ == '__main__':
    main()
