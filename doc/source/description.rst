Introduction
============

**smallrnaseq** is a Python package for processing of small RNA seq data.
This is used to elucidate the small RNA contents of deep sequencing reads.
For non Python users there is a command line interface that is quite simple to use.
Typically a lot of disparate tools are integrated to create pipelines for this kind
of analysis. This can be quite cumbersome. Our objective is to perform the
analyses using Python packages as far as possible allowing almost all requirements
to be installed using the pip tool. An aligner such as bowtie are still needed.

Python users may find the various modules useful in creating their own more flexible workflows.
A functional approach is used with a flat hierarchy of modules with limited use of classes.

The primary source for this documentation is on the github wiki pages. This site
is mostly intended for documentating the module API. For dataexplore see below.

Command Line Interface
===========================

Installing the package provides the command *smallrnaseq* in your path.
This allows users is a command line interface to the library without the need
for any Python coding at all. It provides a set of pre-defined functions with
parameters specified in a text configuration file. This is documented on the github wiki.

Screencast
=====

A screencast tutorial for using the command line interface: https://www.youtube.com/watch?v=m24cuLyTqg0

Citation
========

If you use this software in your work please cite the following article:

**Farrell, D. (2017). smallrnaseq : short non coding RNA-seq analysis with Python.
Bioarxiv. https://doi.org/10.1101/110585**

Installation
============

Linux
-----

On most linux operating systems installations of Python should include the pip tool.
If not use your distributions package manager to install pip first. Then the simple
call below should install all dependencies. However if this fails see the linux section
below for commands to run for installing any pre-requisites that might not be on your
system. We hope to provide a snap package that will simplify installation on linux.

``pip install smallrnaseq``

For python 3 installs
+++++++++++++++++++++

You may need to use the command pip3 instead if python 2 is also on your system, like in Ubuntu.
When installing packages with apt you likely need to specify python 3. e.g. python3-numpy instead
of python-numpy.

For python 2.7 ONLY
+++++++++++++++++++

You might also need the future package. Run `pip install future` to install.

If pip fails you can run the following commands first to fix likely missing packages.
These are mainly needed for HTSeq to install. Then run pip again.

Ubuntu
++++++

    sudo apt install python-dev samtools BEDtools liblzma-dev libbz2-dev zlib1g-dev python-scipy
    sudo pip install smallrnaseq
    sudo apt install bowtie

Fedora
++++++

    sudo dnf install zlib-devel bzip2-devel xz-devel samtools swig redhat-rpm-config python-devel
    sudo pip install cython pysam
    sudo pip install smallrnaseq
    sudo dnf install bowtie

Mac OSX
-------

You will need to make sure you have Python. Anaconda is recommended as it provides most of the
package requirements built in. Follow these steps in order:

Download and run the Mac OS X installer from https://www.continuum.io/downloads.
The installer will automatically configure your system to use the Anaconda Python.
Close the terminal and start a new one.

You should then add the bioconda channels::

    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda

Type this command to install the remaining requirements::

    conda install pybedtools bx-python HTseq

smallrnaseq can then be installed using pip::

    pip install smallrnaseq

bowtie on OSX
+++++++++++++
You can install bowtie with conda too but it may hang or give an error on the latest
version of OSX (Sierra). Running ``conda install bowtie=1.1.2``
will install the older version which should work.

Windows
-------

In theory this package will work on Windows but has not been tested. If you are a windows user
it is recommended to use linux running inside virtualbox.
See http://www.makeuseof.com/tag/how-to-use-virtualbox/

Vienna package
--------------

This is needed if you want to do novel miRNA prediction. It has to be installed separately
on all systems. Go to https://www.tbi.univie.ac.at/RNA/#binary_packages and download the
binary for your system.

Required dependencies

* numpy
* pandas
* matplotlib
* seaborn (requires scipy)
* HTSeq
* scikit-learn