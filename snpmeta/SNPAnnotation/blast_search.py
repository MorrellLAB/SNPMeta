#!/usr/bin/env python
"""A class to store information related to BLAST searching."""

import sys
import tempfile
import subprocess
import time
from urllib.error import HTTPError
from urllib.error import URLError
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastxCommandline


class BlastSearch(object):
    """A class to store information related to BLAST searching. Will store
    database names and paths, run BLAST, and return results."""

    def __init__(self, prog, evalue, maxhits, db, entrez):
        self.prog = prog
        self.evalue = evalue
        self.maxhits = maxhits
        self.db = db
        if not self.db:
            self.web = True
        else:
            self.web = False
        self.entrez = entrez

        self.commandline = None
        self.blastout = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='SNPMeta_BlastOut_',
            suffix='.xml',
            delete=False)
        self.blastin = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='SNPMeta_BlastIn_',
            suffix='.fasta',
            delete=False)
        return

    def build_commandline(self, query):
        """Build the command line based on the arguments that were provided.
        If local database is provided, will create a command line based on the
        Ncbi____Commandline function, based on which program was specified."""
        #   We build a dictionary of command lines, which we will use to select
        #   the command to run.
        command_dict = {
            'blastn': NcbiblastnCommandline(
                query=self.blastin.name,
                out=self.blastout.name,
                db=self.db,
                evalue=self.evalue,
                outfmt=5,
                max_target_seqs=self.maxhits),
            'tblastx': NcbitblastxCommandline(
                query=self.blastin.name,
                out=self.blastout.name,
                db=self.db,
                evalue=self.evalue,
                outfmt=5,
                max_target_seqs=self.maxhits)
            }
        if not self.web:
            #   Write the contents of the query sequence into the temp FASTA
            #   file. Unfortunately, command line BLAST only accepts input
            #   files and not sequences
            SeqIO.write(query, self.blastin, 'fasta')
            self.commandline = command_dict[self.prog]
            self.blastin.close()
        return

    def run_blast(self, query):
        """Runs BLAST using the given query sequence. Reads the flags set in
        self.web and self.entrez to execute the appropriate command, and returns
        a handle to the XML results."""
        #   Set the BLAST output handle to None
        blast_handle = None
        #   If we are running over the web
        if self.web:
            #   Sometimes, we can fail to execute BLAST over the web, especially
            #   in times of high traffic. We wrap the web queries in a
            #   try/except block, and retry if we get an error, up to a max of
            #   three tries.
            tries = 0
            success = False
            while not (success or tries >= 3):
                tries += 1
                try:
                    #   And then, if we are limited by an Entrez query
                    if self.entrez:
                        blast_handle = NCBIWWW.qblast(
                            self.prog,
                            'nr',
                            query.seq,
                            format_type='XML',
                            descriptions=self.maxhits,
                            alignments=self.maxhits,
                            hitlist_size=self.maxhits,
                            expect=self.evalue,
                            entrez_query=self.entrez)
                    else:
                        blast_handle = NCBIWWW.qblast(
                            self.prog,
                            'nr',
                            query.seq,
                            format_type='XML',
                            descriptions=self.maxhits,
                            alignments=self.maxhits,
                            hitlist_size=self.maxhits,
                            expect=self.evalue)
                    success = True
                except HTTPError as e:
                    #   If we catch an HTTPError, print the error code, and wait
                    #   5 seconds to try again.
                    sys.stderr.write(
                        'Caught HTTP error ' +
                        str(e.code) +
                        ' with reason ' +
                        str(e.reason) +
                        '. Retrying in 5 seconds ...\n')
                    time.sleep(5)
                except (URLError, ConnectionRefusedError) as e:
                    #   Sometimes our connection is refused. We handle this by
                    #   retrying after 5 seconds. The reason for this error is
                    #   usually that the service that listens on the requested
                    #   port can't handle the request.
                    sys.stderr.write(
                        'Caught URL error, connection refused. Retrying in 5 '
                        'seconds ...\n')
                    time.sleep(5)
                finally:
                    #   Print a helpful message here about how the web BLAST
                    #   went.
                    if not success and tries >= 3:
                        sys.stderr.write(
                            'Failed to get BLAST results after three tries!\n'
                            'Moving on to next SNP ... \n')
        else:
            #   We are not running over the web, and will execute a local BLAST
            #   command.
            subprocess.call(str(self.commandline).split())
            blast_handle = open(self.blastout.name, 'r')
        return blast_handle
