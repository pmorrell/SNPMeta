#!/usr/bin/env python
"""Contains a class to store and print SNP annotation information."""

import datetime


class SNPAnnotation(object):
    """A class to store and emit the SNP annotation data."""

    def __init__(self, sequence):
        """Initalize the class with the annotation time, and the SNP name.
        These instance variables correspond to table headers in the tabular
        output format."""
        self.SNPName = sequence.id
        self.FiveFlank = '-'
        self.ThreeFlank = '-'
        self.Observed = '-'
        self.Organism = '-'
        self.GenBankID = '-'
        self.ProteinID = '-'
        self.GeneShortName = '-'
        self.Position = '-'
        self.ThreeUTR = '-'
        self.FiveUTR = '-'
        self.Silent = '-'
        self.AA1 = '-'
        self.AA2 = '-'
        self.Grantham = '-'
        self.CDSPos = '-'
        self.Codon1 = '-'
        self.Codon2 = '-'
        self.Ambiguity = '-'
        self.Ancestral = '-'
        self.Product = '-'
        self.Notes = '-'
        self.AltGene = '-'
        self.AltOrg = '-'
        self.ContextSeq = '-'
        self.AlignScore = '-'
        self.DateTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        #   These attributes are merely in support of gathering the above
        #   information.
        self.query_seq = sequence.seq
        self.context = '-'
        return


    def calculate_clen(self, s, l, gbs, illumina):
        """Calculate the contextual sequence length, given the full SNP
        sequence, the estimated length, whether it came from GBS tags, and/or
        whether it came from an Illumina array."""
        pass

    def print_dbsnp(self, out):
        """Print the SNP annotation information in dbSNP report format."""
        pass

    def print_tabular(self, out):
        """Print the SNP annotation information in tabular format."""
        pass
