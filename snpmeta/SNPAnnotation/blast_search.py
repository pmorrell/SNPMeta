#!/usr/bin/env python
"""A class to store information related to BLAST searching."""

import tempfile
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline


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
        return

    def build_commandline(self, query):
        """Build the command line based on the arguments that were provided.
        If local database is provided, will create a command line based on the
        Ncbi____Commandline function, based on which program was specified."""
        #   We build a dictionary of command lines, which we will use to select
        #   the command to run.
        command_dict = {
            'blastn': NcbiblastnCommandline(
                query=query,
                out=self.blastout.name,
                db=self.db,
                evalue=self.evalue,
                outfmt=5,
                num_descriptions=self.maxhits,
                num_alignments=self.maxhits)
            'tblastx': NcbitblastxCommandline(
                query=query,
                out=self.blastout.name,
                db=self.db,
                evalue=self.evalue,
                outfmt=5,
                num_descriptions=self.maxhits,
                num_alignments=self.maxhits)
            }
        if not self.web:
            self.commandline = command_dict[self.prog]
        return

    def run_blast(self, query):
        """Runs BLAST using the given query sequence. Reads the flags set in
        self.web and self.entrez to execute the appropriate command, and returns
        a handle to the XML results."""
        #   If we are running over the web
        if self.web:
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
        else:
            #   We are not running over the web, and will execute a local BLAST
            #   command.
            subprocess.call(str(blast_command).split())
            blast_handle = open(self.blastout.name, 'r')
        return blast_handle
