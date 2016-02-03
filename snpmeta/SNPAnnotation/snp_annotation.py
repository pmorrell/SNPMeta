#!/usr/bin/env python
"""Contains a class to store and print SNP annotation information."""

import datetime
from snpmeta.LookupTables.iupac import IUPAC
from Bio.Seq import Seq
from Bio import AlignIO


#   The verbose table header
table_header = [
    'SNPName',
    'Organism',
    'GenBankID',
    'ProteinID',
    'GeneShortName',
    'Position',
    'Downstream',
    'Upstream',
    'Silent',
    'AA1',
    'AA2',
    'GranthamScore',
    'CDSPosition',
    'Codon1',
    'Codon2',
    'AmbiguityCode',
    'ProductName',
    'Notes',
    'RelatedGene',
    'RelatedOrganism',
    'ContextSequence',
    'AlignScore',
    'DateTime'
]


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
        self.aligned_seq = None
        self.aligned_pos = None
        return


    def calculate_clen(self, cl, gbs, illumina):
        """Calculate the contextual sequence length, given the full SNP
        sequence, the estimated length, whether it came from GBS tags, and/or
        whether it came from an Illumina array."""
        if cl:
            #   The first check we perform is if the contextual sequence length
            #   was passed. If it was, we check the query sequence to see if it
            #   makes sense, and set it.
            if len(self.query_seq) > cl:
                #   check that the SNP sequence length is greater than the
                #   supplied context length. If it's not, there is an error.
                if self.query_seq[cl] in IUPAC:
                    #   First check from the left to see if the contextual
                    #   sequence length gives us an ambiguity
                    self.context = cl
                elif self.query_seq[-(cl + 1)] in IUPAC:
                    #   Next, check from the right. If this does give us an
                    #   ambiguity, convert the contextual length to a left-
                    #   based offset.
                    self.context = len(self.query_seq) - (cl + 1)
                else:
                    #   A last resort is to iterate through until we find the
                    #   first ambiguity. This is almost guaranteed to be wrong
                    for i, b in enumerate(self.query_seq):
                        if b in IUPAC:
                            self.context = i
                            break
        elif illumina:
            #   The next flag we check is whether or not the SNPs are given with
            #   Illumina formatting. If so, we search for the [ character that
            #   denotes the query states, and find the offset. We also remove
            #   the non-nucleotide characters from the sequence, as those will
            #   cause errors later.
            offset = re.search(r'\[', str(self.query_seq)).start()
            #   Build the fixed sequence
            seq_parts = re.split(r'[\[|\/|\]]', str(self.query_seq))
            for a, nona in IUPAC.iteritems():
                #   Search for the ambiguity that contains the two query SNP
                #   states, and drop it in
                if (seq_parts[1] in nona) and (seq_parts[2] in nona):
                    snp = a
                    break
            self.query_seq = Seq(seq_parts[0] + snp + seq_parts[-1])
            self.context = offset
        elif gbs:
            #   The final check is whether the SNPs are from GBS tags or not. If
            #   they are, we iterate until we find the ambiguity, and return the
            #   offset of the first one we find.
            for i, b in enumerate(self.query_seq):
                if b in IUPAC:
                    self.context = i
                    break
        return

    def get_aligned_pos(self, alignment):
        """Get the position of the SNP in the alignment."""
        #   First, extract the score of the alignment. This is stored on a line
        #   with the form
        #       #Score: xxx.x
        with open(alignment, 'r') as f:
            for line in f:
                if line.startswith('#Score:'):
                    self.align_score = line.split(':')[-1].strip()
                    break
        aln = AlignIO.read(alignment, 'emboss')
        #   First sequence is the GenBank sequence, we store it so that we can
        #   extract proper features from it later
        self.aligned_seq = aln[0].seq
        #   Don't bother storing the aligned SNP sequence long-term. We are just
        #   going to REs on it to calculate positions.
        snp_seq = str(temp_align[1].seq)
        snp_start = re.search('[^-]', snp_seq).start()
        #   Match the last instance of a non-gap character
        #   Two situations are:
        #       NNNNNNNNNN------$
        #   or  NNNNNNNNNNNNNNNN$
        #   +1 because of 0-based indexing
        search_pattern = '([ATCGRYMWKSN]-+$)|([ATGCRYMKWSN]$)'
        snp_end = re.search(search_pattern, snp_seq, re.IGNORECASE).start() + 1
        #   The complete range the aligned SNP sequence convers
        snp_coverage = range(snp_start, snp_end)
        #   Extract the relevant portion of the SNP sequence
        aligned_snp = snp_seq[snp_start:snp_end]
        #   These are the cases that are possible in the alignment:
        #   AGGATG-(*k)TTAGAGAT... [SNP] AGGAGTCGATG... | SNP with gap on left
        #   GATGGATGCCGTG... [SNP] AGCCTCAT(-*k)GATG... | SNP with gap on right
        #   GAAC(-*c)GACC... [SNP] AGTAGC(-*k)GGATGC... | SNP with gaps on both
        #   We remove the gaps from the sequence, and also remove those from
        #   the mapping positions
        coverage_pos = zip(snp_coverage, aligned_snp)
        actual_positions = [base[0] for base in coverage_pos if base[1] != '-']
        #   Count the number of expected bases in from the end to get the
        #   position of the SNP in the alignment
        self.aligned_pos = actual_positions[self.context]
        return

    def functional_class(self, gbseq, features):
        """Calculates where the SNP is coding or noncoding, and if applicable,
        the amino acid states."""
        pass

    def print_dbsnp(self, out):
        """Print the SNP annotation information in dbSNP report format."""
        pass

    def print_tabular(self, out):
        """Print the SNP annotation information in tabular format."""
        print(table_header)
        return
