#!/usr/bin/env python

"""
    Module for miRDeep2 wrappers and utilities
    Created July 2014
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

import sys, os, string, types, re
import shutil, glob, collections
import itertools
import subprocess
import pylab as plt
import numpy as np
import pandas as pd
from . import base, utils

mirdeep2options = {'base': [('input',''),('adapter','TGGAATTCTCGGGTGCCAAGG'),('filetype','fastq'),
                    ('bowtieindex',''),('refgenome',''),('species','hsa'),
                    ('mature',''), ('hairpin',''), ('other',''),('mirbase',os.getcwd()),
                    ('randfold',1), ('overwrite',1)]}
mirdeepcols = ['#miRNA','read_count','mean_norm','miRDeep2 score','chr','seed','precursor',
                'freq','precursor coordinate','mirbase seed match','star read count','rfam alert',
                'consensus mature sequence','consensus star sequence',
                'consensus precursor sequence']

def create_mirbase_files(species,path):
    """Generate species specific mature/hairpin files for input to mirdeep"""

    mature = os.path.join(path, 'mature.fa')
    hairpin = os.path.join(path, 'hairpin.fa')
    names=[]
    for f in [mature,hairpin]:
        fname = os.path.splitext(f)[0]+'_'+species+'.fa'
        utils.get_subset_fasta(f, labels=[species], outfile=fname)
        names.append(fname)
    print ('wrote mirbase files for species %s' %species)
    return names

def create_sample_map(path, ext='fastq'):
    """Create filename mapping to run all samples at once.
       This is required."""

    os.chdir(path)
    files = sorted(glob.glob('*.'+ext))
    print (files)
    fname = 'combined.txt'
    i=1
    rows=[]
    for f in files:
        rows.append((os.path.basename(f), 's%02d'%i))
        i+=1
    res=pd.DataFrame(rows)
    res.to_csv(fname, index=False,sep=' ',header=False)
    return (fname)

def combine_labels(labels, filename):
    """Combine sample labels with mirdeep combined.txt file so we can match sample ids"""

    comb = pd.read_csv(filename,sep=' ',names=['filename','id'])
    comb['name'] = comb.filename.apply (lambda x: x.split('.')[0])
    labels = labels.merge(comb,on='name')
    labels['id'] = labels.id+'(norm)'
    return labels

def run_multiple(**kwargs):
    """Prepare and run mirdeep2"""

    if kwargs['filetype'] == 'fasta':
        ext='fa'
    else:
        ext='fastq'
    path = kwargs['input']
    #create filename/id mapping
    samplemap = create_sample_map(path, ext)
    #get mirbase subset for species if provided
    if kwargs['species'] != '':
        mature, hairpin = create_mirbase_files(kwargs['species'], kwargs['mirbase'])
        kwargs['mature'] = mature
        kwargs['hairpin'] = hairpin
    if kwargs['other'] != '':
        kwargs['other'], h = create_mirbase_files(kwargs['other'], kwargs['mirbase'])
    run(samplemap, **kwargs)
    return

def run(infile, refgenome, bowtieindex, mature='', hairpin='', other='',
        randfold=False, overwrite=True, filetype='fastq', adapter=None,
        clean=True, outpath=None, **kwargs):
    """Run all mirdeep2 steps including adapter trimming.
       Uses a config file even if we only have one sample."""

    label = os.path.splitext(os.path.basename(infile))[0]
    print ('running %s' %label)
    os.environ["BOWTIE_INDEXES"] = os.path.dirname(bowtieindex)
    collapsed = 'collapsedreads.fa'
    if filetype=='fasta': mapparams='-c'
    else: mapparams='-e -h'
    if adapter =='': adapter = 'none'
    if randfold == False: params='-c'
    else: params = ''
    if other=='': other = 'none'

    #if mapping has been done already we can skip it
    if not os.path.exists('mapped.arf') or overwrite == True:
        try:
            os.remove('mapped.arf')
            os.remove(collapsed)
        except:
            pass
        cmd1 = ('mapper.pl %s -d %s -j -l 18 -m -k %s -s %s'
                ' -p %s -t mapped.arf -v' %(infile,mapparams,adapter,collapsed,bowtieindex))
        print (cmd1)
        result = subprocess.check_output(cmd1, shell=True, executable='/bin/bash')
    else:
        print ('arf file found, skipping mapper step')

    #mirdeep core
    cmd2 = ('miRDeep2.pl %s %s mapped.arf'
           ' %s %s %s -z _%s %s'
           ' -d > report.log' %(collapsed,refgenome,mature,other,hairpin,label,params))
    print (cmd2)
    result = subprocess.check_output(cmd2, shell=True, executable='/bin/bash')
    #remove junk
    if clean == True:
        tempfiles = glob.glob('dir_*')
        for f in tempfiles:
            shutil.rmtree(f)
    #move results to dest folder
    if outpath != None:
        pass #move files..
    return

def quantifier(path, mature, precursor, star=None, collapsed='collapsedreads.fa', time='novel'):
    """Run quantifier module using custom known mature/precursors"""

    current = os.getcwd()
    os.chdir(path)
    cmd = 'quantifier.pl -p %s -m %s -r %s -y %s -k -d -g 1 -U' %(precursor,mature,collapsed,time)
    print (cmd)
    result = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    os.chdir(current)
    return

def get_pdf_path(path):
    """get path to pdfs"""
    ppath = glob.glob(os.path.join(path,'pdf*'))[0]
    return ppath

def get_chromosome(x):
    val = x.split('_')[0]
    try:
        return '%02d' %int(val)
    except:
        return val

def get_coords(x):
    """Get start/end from precursor coords string"""

    l = re.split("[\:\..]+",x)
    return pd.Series(l[:4],index=['chr','start','end','strand'])

def get_score_stats(path):
    """Get mirdeep results from the summary file"""

    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    if os.path.splitext(resfile)[1] != '.csv':
        return
    df = pd.read_csv(resfile, sep='\t',header=0,nrows=22)
    df = df.dropna()
    df['known'] = df['known miRNAs detected by miRDeep2'].apply(lambda r: int(r.split()[0]))
    df['FDR'] = df['novel miRNAs, estimated false positives'].apply(lambda r: int(r.split()[0]))
    return df

def plot_score_stats(df):

    f,axs=plt.subplots(2,2,figsize=(10,6))
    grid=axs.flat
    df.plot('miRDeep2 score','known',marker='o',ax=grid[0],legend=False)
    grid[0].set_title('known miRNAs')
    df.plot('miRDeep2 score','novel miRNAs reported by miRDeep2',marker='o',ax=grid[1],legend=False)
    grid[1].set_title('novel miRNAs')
    df.plot('miRDeep2 score','estimated signal-to-noise',marker='o',ax=grid[2],legend=False)
    grid[2].set_title('signal-to-noise')
    df.plot('miRDeep2 score','FDR',marker='o',ax=grid[3],legend=False)
    grid[3].set_title('FDR')
    plt.tight_layout()
    f.savefig('mirdeep_score_stats.png')
    return

def read_results_file(infile):
    """Get mirdeep results from the summary file"""

    if os.path.splitext(infile)[1] != '.csv':
        return
    scol = 'miRDeep2 score'
    df = pd.read_csv(infile, sep='\t',header=23)
    df = df.dropna()
    #remove junk
    idx = df[df['provisional id']=='tag id'].index[0]
    df = df.drop(idx)
    df['novel'] = np.where(df['miRBase miRNA']=='-',True,False)
    colstorename = {'example miRBase miRNA with the same seed':'mirbase seed match',
         'significant randfold p-value': 'randfold'}
    df = df.rename(columns=colstorename)
    df = df.convert_objects(convert_numeric=True)
    #df['chr'] = df['provisional id'].apply(get_chromosome)
    df['seed'] = df['consensus mature sequence'].apply(lambda x: x[1:8])
    coords = df['precursor coordinate'].apply(get_coords)
    df = df.join(coords) #pd.concat([df,coords])
    return df

def get_results(path):
    """Process known and novel results from mirdeep run.
       Combines with expression data and removes redundant entries"""

    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    df = read_results_file(resfile)

    #use quantifier module to get novel expression results from predicted precursors
    #if not done already
    novelmature = os.path.join(path,'novel_mature.fa')
    novelprecursor = os.path.join(path,'novel_precursor.fa')
    res = os.path.join(path, 'expression_novel.html')
    reads = 'collapsedreads.fa'
    if not os.path.exists(res):
        novel = df[df.novel==True]
        mkey = 'consensus mature sequence'
        pkey = 'consensus precursor sequence'
        utils.dataframe_to_fasta(novel, mkey, 'provisional id', outfile=novelmature)
        utils.dataframe_to_fasta(novel, pkey, 'provisional id', outfile=novelprecursor)
        quantifier(path, os.path.abspath(novelmature), os.path.abspath(novelprecursor))

    #get expression results and merge with prediction results to get other info
    files = glob.glob(os.path.join(path,'miRNAs_expressed_all_samples*.csv'))
    res=[]

    for f in files:
        if 'novel' in f:
            key='provisional id'
        else:
            key='miRBase miRNA'
        q = pd.read_csv(f,sep='\t')
        samples = float(len(q.filter(regex="norm").columns))
        #print 'samples: %s' %samples
        q['freq'] = q.filter(regex="norm").apply(lambda r: len(r.nonzero()[0])/samples,1)
        #apply 5p id so we can merge with results file and keep star seqs
        q['id'] = q['#miRNA'].apply(lambda x: x[:-2]+'5p' if str(x).endswith('3p') else x)
        #loses information on multiple precursors for a mature seq
        q = q.merge(df,left_on=['id'],right_on=key).drop_duplicates('#miRNA')
        res.append(q)
    res = pd.concat(res)

    #get mean normalised count
    res['mean_norm'] = res.filter(regex="norm").apply(lambda r: r[r.nonzero()[0]].mean(),1)
    res = res.sort_values(by=['read_count'], ascending=False)
    res = res.drop_duplicates('#miRNA')
    res = res.reset_index(drop=True)
    #res['std'] = res.filter(regex="norm").std(1)
    #res['cv'] = res['std']/res['mean_norm']
    return res

def get_column_names(df):
    """Extract column names for multiple samples"""

    cols = [i for i in df.columns if (i.startswith('s') and len(i)<=3)]
    normcols = [i+'(norm)' for i in cols]
    return cols, normcols

def filter_expr_results(df, cols=None, score=0, freq=0.5, mean_norm=0, total_reads=0):
    """Additional filters for abundances/no. samples"""

    if cols is None:
        cols = [i for i in df.columns if (i.startswith('s') and len(i)<=3)]
    normcols = [i+'(norm)' for i in cols]
    df = df[(df['miRDeep2 score']>=score)]
    df = df[df.freq>=freq]
    df = df[df['read_count']>=total_reads]
    df = df[df['mean_norm']>=mean_norm]
    df = df[df['randfold'].isin(['yes','-'])]
    #n = n[n['rfam alert']=='-']
    #df = df.reset_index(drop=True)
    return df

def analyse_results(path, outpath=None, **kwargs):
    """General analysis of mirdeep results"""

    if outpath != None:
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        os.chdir(outpath)
    df = get_results(path)
    idcols,normcols = get_column_names(df)
    known = df[df.novel==False]
    novel = df[df.novel==True]
    idmap = get_file_ids(path)

    #cutoffs and freqs need to be configurable..
    k = filter_expr_results(known,score=0,freq=.5,total_reads=50)
    n = filter_expr_results(novel,score=4,freq=.8,total_reads=50)
    cols = mirdeepcols
    core = pd.concat([k,n])
    utils.dataframe_to_fasta(core, 'consensus mature sequence', '#miRNA', 'mirdeep_core.fa')

    k[cols].to_csv('known_mirdeep.csv')
    n[cols].to_csv('novel_mirdeep.csv')
    utils.create_html(n[cols],'novel_mirdeep')
    k['perc'] = k['read_count']/k['read_count'].sum()

    print (k[cols[:8]])
    print
    print (n[cols[:9]])

    print ('mirdeep summary')
    print ('-------------------------------')
    print (('%s (%s novel) identified' %(len(df),len(novel))))
    print ('quantifier results after filtering:')
    print ('%s/%s known' %(len(k),len(known)))
    print ('%s/%s novel' %(len(n),len(novel)))
    print ('top 10 known account for %2.2f' %k['perc'][:10].sum())
    print ('%s are 3p strand' %len(k[k['#miRNA'].str.contains('3p')]))
    #print df[(df.mean_norm>300) & (df.freq<0.5)][cols[:7]]
    print

    k=k.set_index('#miRNA')
    n=n.set_index('#miRNA')

    fig, ax = plt.subplots(figsize=(8,6))
    k['read_count'][:10].plot(kind='barh',colormap='Spectral',ax=ax,log=True)
    plt.title('miRDeep2 top 10')
    plt.tight_layout()
    fig.savefig('mirdeep_top_known.png',dpi=150)
    fig, ax = plt.subplots(figsize=(8,8))
    df.plot('freq','mean_norm',kind='scatter',ax=ax,logy=True,alpha=0.8)
    fig.savefig('mirdeep_freqsvcounts.png')
    fig=plt.figure()

    fig = plot_read_count_dists(n,h=5)
    fig.savefig('mirdeep_novel_counts.png')
    fig = plot_read_count_dists(k)
    fig.savefig('mirdeep_known_counts.png')

    #perSampleDists(k)

    fig,ax = plt.subplots(figsize=(10,6))
    core[idcols].sum().plot(kind='bar',ax=ax)
    plt.title('total miRNA counts per sample (unnormalised)')
    plt.tight_layout()
    fig.savefig('mirdeep_total_persample.png')

    #found per chromosome
    fig = plt.figure(figsize=(8,6))
    x = core.sort('chr').groupby('chr').size()
    x.plot(kind='bar')
    plt.title('miRNA per chromosome')
    fig.savefig('mirdeep_chromosome_dist.png')

    ss = get_score_stats(path)
    plot_score_stats(ss)
    #plt.show()
    #plt.close()
    return df,k,n

def plot_read_count_dists(df,h=8):
    """Boxplots of read count distributions per miRNA - use seaborn?"""

    w=int(h*(len(df)/60.0))+4
    fig, ax = plt.subplots(figsize=(w,h))
    cols,normcols = get_column_names(df)
    df = df[normcols]
    t=df.T
    t.index = cols
    try:
        import sns
        sns.boxplot(t,linewidth=1.0,saturation=0.2,palette='coolwarm_r')
        sns.despine(trim=True)
    except:
        t.plot(kind='box',color='black',grid=False,whis=1.0,ax=ax)
    ax.set_yscale('log')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.ylabel('read count')
    plt.tight_layout()
    return fig

def get_file_ids(path):
    """Get file<->mirdeep2 id mapping"""

    idmap = pd.read_csv(os.path.join(path, 'combined.txt'),sep=' ',
                header=None,names=['filename','id'])
    return idmap

def get_label_map(path, labels):
    """Get results labels mapped to labels with the filenames"""

    condmap = pd.read_csv(labels)
    idmap = get_file_ids(path)
    def matchname(x):
        r = idmap[idmap['filename'].str.contains(x)]
        if len(r)>0:
            return r.id.squeeze()
        else:
            return np.nan
    condmap['id'] = condmap.filename.apply(lambda x: matchname(x),1)
    condmap = condmap.dropna()
    return condmap

def test_quantifier(path):
    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    df = read_results_file(resfile)
    novelmature = os.path.join(path, 'novel_mature.fa')
    novelstar = os.path.join(path, 'novel_star.fa')
    novelprecursor = os.path.join(path, 'novel_precursor.fa')
    reads = '../results_mirdeep_combined/collapsedreads.fa'
    #base.dataframe2Fasta(df[df.novel==True], 'consensus star sequence', 'provisional id',
    #            outfile=novelstar)
    quantifier(path, os.path.abspath(novelmature), os.path.abspath(novelprecursor),
                    os.path.abspath(novelstar), reads)
    return

def check_quantifier_results(path):
    """Check quantifier vs results file in case of miscounts"""

    resfile = glob.glob(os.path.join(path,'result*.csv'))[0]
    df = read_results_file(resfile)
    files = glob.glob(os.path.join(path,'miRNAs_expressed_all_samples*.csv'))
    q = pd.read_csv(files[0],sep='\t')
    key='provisional id'
    m=q.merge(df,left_on='#miRNA',right_on=key).drop_duplicates('#miRNA')
    m.sc = m['miRDeep2 score']
    m['err'] = abs(m['read_count']-m['total read count'])
    cols=['#miRNA','total read count','read_count','miRDeep2 score']
    print (m[m.err>400].sort('total read count',ascending=False)[cols])
    m['size'] = np.select([m.sc < 2, m.sc < 3, m.sc < 4], [20,40,50], 80)
    f,ax=plt.subplots(1,1)
    plt.xscale('log')
    plt.yscale('log')
    m.plot(x='total read count',y='read_count', kind='scatter',s=60,alpha=0.6,ax=ax)
    #ax.plot([0, 1], [0, 1], transform=ax.transAxes,color='red',alpha=0.7)
    plt.show()
    return

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-r", "--run", dest="run", action='store_true',
                           help="run predictions")
    parser.add_option("-i", "--input", dest="input",
                           help="input path or file")
    parser.add_option("-c", "--config", dest="config",
                            help="config file")
    parser.add_option("-a", "--analyse", dest="analyse",
                           help="analyse results of mirdeep2")

    opts, remainder = parser.parse_args()
    pd.set_option('display.width', 800)

    if opts.run == True:
        #all other options are stored in config file
        if opts.config == None:
            print ('No config file provided.')
            base.write_default_config('mirdeep2.conf', defaults=mirdeep2options)
            return
        cp = base.parse_config(opts.config)
        if opts.input != None:
            conf.input = os.path.abspath(opts.input)
        options = base.get_options(cp)
        run_multiple(**options)
    elif opts.analyse != None:
        analyse_results(opts.analyse)
    elif opts.test == True:
        test(opts.input)

if __name__ == '__main__':
    main()
