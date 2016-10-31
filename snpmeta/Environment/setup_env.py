#!/usr/bin/env python
"""Functions to prepare the working environment of SNPMeta"""

import os
import sys
from Bio import SeqIO


def get_dir_contents(directory, noblast):
    """Returns the list of files to work on. Depends on if --no-blast was given
    on the command line."""
    #   If --no-blast was supplied...
    if noblast:
        #   If --no-blast was specified, operate on files ending in 'xml'
        suffix = 'xml'
    else:
        #   If not, then we want the FASTA to BLAST them
        suffix = 'fasta'
    #   The full path of the directory supplied
    #   We need this since the script looks relative to where it is located
    #   on the disk.
    full_path = os.path.abspath(directory)
    #   This is the full directory listing of the current directory, as a list
    file_list = os.listdir(full_path)
    #   Trim the file list by those ending in the proper suffix
    file_list = [x for x in file_list if x.endswith(suffix)]
    return file_list


def build_targets(fasta, directory, noblast):
    """Returns a list of annotation targets. This will return either a list of
    SeqRecord objects, or a list of XML reports to use for annotation."""
    #   Are we annotating a directory or a single file?
    if directory:
        #   Get the list of files we are going to annotate
        record_list = get_dir_contents(directory, noblast)
        #   Drop a message if the directory is empty
        #   An empty list evaluates to Boolean False
        if not record_list:
            sys.stderr.write(
                'No SNPs or BLAST reports in specified directory!\n')
            exit(1)
    elif fasta:
        #   Get the generator object for the records in the FASTA file
        #   We use BioPython's SeqIO.parse() function, which returns
        #   a generator that spits out one FASTA record at a time
        record_list = list(SeqIO.parse(fasta, 'fasta'))
    return record_list
