#!/usr/bin/env python
"""A class to store information related to BLAST searching."""


class BlastSearch(object):
    """A class to store information related to BLAST searching. Will store
    database names and paths, run BLAST, and return results."""

    def __init__(self, prog, evalue, maxhits, db, entrez):
        self.prog = prog
        self.evalue = evalue
        self.maxhits = maxhits
        self.db = db
        self.entrez = entrez

        self.commandline = None
        return

    def build_commandline(self):
        """Build the command line based on the arguments that were provided.
        If local database is provided, will create a command line based on the
        Ncbi____Commandline function, based on which program was specified. If
        a local database was not provided, it will set a flag to run BLAST over
        the network against NCIB's servers."""
        pass

    def run_blast(self, query):
        """Runs BLAST using the given query sequence. Reads the value out of
        self.commandline and runs BLAST accordingly. If it is set to web-blastn,
        it will run BLASTN over the network, if it is set to web-tblastx, it
        will run TBLASTX over the network, else it will execute the command that
        is present. Returns a handle to the XML results."""
        pass
