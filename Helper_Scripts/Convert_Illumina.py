#!/usr/bin/env python

"""
A SNPMeta companion script to convert from the Illumina format to FASTA,
acceptable as input for SNPMeta.
"""
###############################################################################
#
#   January 25, 2013
#   Thomas Kono
#   Saint Paul, MN
#   
#   Dependencies:
#       None
#
###############################################################################

#   To read arguments off the command line
import sys
#   To match regular expressions
import re
#   To check files
import os

#   The lookup table for IUPAC ambiguitites. BioPython certainly does this
#   but their implementation is cumbersome. This makes for REALLY SIMPLE lookup
#   and takes up little memory.
IUPAC = {
        "A": {
            "T": "W",
            "C": "M",
            "G": "R",
            },
        "T": {
            "A": "W",
            "C": "Y",
            "G": "K",
            },
        "C": {
            "A": "M",
            "T": "Y",
            "G": "S",
            },
        "G": {
            "A": "R",
            "T": "K",
            "C": "S",
            },
        }

#   Check the arguments given
if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " [FILE]"
    print """Converts an Illumina-formatted contextual sequence file to FASTA,
suitable for input to SNPMeta. [FILE] must contain two fields, separated
by a tab. The first field is the name of the SNP, and the second is the 
sequence."""
    exit(1)
try:
    #   If we can't get stats on it, then we can't read it for one reason or
    #   another.
    tmp_stats = os.stat(sys.argv[1])
except OSError:
    print "Error! The file provided is not readable, or does not exist!"
    exit(1)

#   Convert the file!
#       First, read it in
with open(sys.argv[1], 'r') as f:
    for line in f:
        #   Strip the newlines, split on tabs
        tmp = line.strip().split('\t')
        #   save the name and the sequence
        name = tmp[0]
        ill_seq = tmp[1]
        #   Search for the ambiguity, it is listed as [N/N]
        search = re.compile('\[[ATCG]\/[ATCG]\]', re.IGNORECASE)
        ill_match = re.search(search, ill_seq)
        ill_amb = ill_match.group(0)
        #   Then, we get the IUPAC ambiguity to drop in
        #   We slice like this because it is in the form
        #       [N/N]
        #   We convert to uppercase, since we can find both lowercase and
        #   uppercase, but out lookup is only for upppercase
        iupac_ambiguity = IUPAC[ill_amb[1].upper()][ill_amb[3].upper()]
        #   Drop it in instead
        fasta_seq = re.sub(search, iupac_ambiguity, ill_seq, re.IGNORECASE)
        #   Now print out the finished sequences
        print '>' + name
        print fasta_seq
