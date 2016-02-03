#!/usr/bin/env python
"""Contains a class to store and print SNP annotation information."""

import datetime


class SNPAnnotation(object):
    """A class to store and emit the SNP annotation data."""

    def __init__(self, sequence):
        """Initalize the class with the annotation time, and the SNP name.
        These instance variables correspond to table headers in the tabular
        output format."""
        self.snp_name = sequence.id
        self.five_flank = '-'
        self.three_flank = '-'
        self.observed = '-'
        self.organism = '-'
        self.genbank_id = '-'
        self.protein_id = '-'
        self.gene_short_name = '-'
        self.position = '-'
        self.three_utr = '-'
        self.five_utr = '-'
        self.silent = '-'
        self.aa_1 = '-'
        self.aa_2 = '-'
        self.grantham = '-'
        self.cds_pos = '-'
        self.codon_1 = '-'
        self.codon_2 = '-'
        self.ambiguity = '-'
        self.ancestral = '-'
        self.product = '-'
        self.notes = '-'
        self.alt_gene = '-'
        self.alt_org = '-'
        self.context_seq = '-'
        self.align_score = '-'
        self.date_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

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
