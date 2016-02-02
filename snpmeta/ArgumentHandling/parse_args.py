#!/usr/bin/env python
"""Parses arguments using the argparse module."""

try:
    import argparse
except ImportError:
    print('Error! You need Python >=3.2 or the \'argparse\' library installed.')
    exit(1)

#   Set a description of the program
DESCR = """Annotates SNPs using sequences from a specified organism in NCBI's
nonredundant sequence database. Outputs metadata in a dbSNP submission format
or a tab-delimited verbose format. Requires Biopython and Python >=3.2."""


def parse_args():
    """Parse the arguments that are passed. See the manual for a more detailed
    description of all the arguments that are available."""
    parser = argparse.ArgumentParser(description=DESCR, add_help=True)

    #   Add a group for required arguments
    required = parser.add_argument_group(title='Required arguments')
    #   First required argument is the user's email address, for sending fetch
    #   requests to NCBI.
    required.add_argument(
        '-a',
        '--email',
        metavar='EMAIL',
        type=str,
        help='E-mail address to send with fetch queries to NCBI.')
    #   Within the required arguments parser, add a mutually exclusive group.
    #   The user will have to specify either a file with SNPs to annotate, or
    #   a directory with single record FASTA files to use for annotation.
    annotation_target = required.add_mutually_exclusive_group()
    annotation_target.add_argument(
        '-f',
        '--fasta-file',
        metavar='FASTA',
        default=None,
        help='A multi-record FASTA file containing SNPs to annotate. Mutually '
             'exclusive with -d.')
    annotation_target.add_argument(
        '-d',
        '--dir',
        metavar='DIR',
        default=None,
        help='A directory containing single-record FASTA files, each with a '
             'SNP to annotate. Mutually exclusive with -f.')

    #   Next, add another mutually exclusive group for specifying the contextual
    #   sequence length around the query SNPs. The user must either specify a
    #   numerical argument for a set contextual sequence length, or a flag that
    #   tells the program to search through Illumina design sequences to find
    #   the contextual sequence length, or that the SNPs are in GBS tags and
    #   there is no set contextual sequence length.
    context_length = parser.add_mutually_exclusive_group()
    context_length.add_argument(
        '-l',
        '--contextual-length',
        metavar='LENGTH',
        type=int,
        help='The length of the contextual sequence surrounding the SNP, if '
             'the SNP is stored as an IUPAC ambiguity in the FASTA file.')
    context_length.add_argument(
        '-g',
        '--gbs',
        action='store_true',
        default=False,
        help='SNP sequences are GBS tags, and the SNP can occur anywhere in '
             'the sequence. Default: False.')
    context_length.add_argument(
        '-i',
        '--illumina',
        action='store_true',
        default=False,
        help='SNP sequences are formatted like FASTA, but the query SNP is '
             'specified as the SNP states enclosed by brackets []. See the '
             'user manual for details on this format. Default: False.')

    #   Define an argument group for BLAST parameters, if we are running BLAST
    #   through SNPMeta. Users can specify either BLASTN or TBLASTX, and E-value
    #   threshold, the maximum number of hits to return, a path to a custom
    #   database, and/or an Entrez query to limit BLAST searches to a sepcific
    #   taxon.
    blast_args = parser.add_argument_group(
        title='BLAST arguments',
        description='Ignored if --no-blast is specified.')
    blast_args.add_argument(
        '-p',
        '--program',
        default='blastn',
        type=str,
        choices=['blastn', 'tblastx'],
        help='The BLAST program to use. Default: blastn')
    blast_args.add_argument(
        '-e',
        '--evalue',
        default='5e-10',
        help='Maximum E-value to accept as a BLAST hit. Default: 5e-10')
    blast_args.add_argument(
        '-m',
        '--max-hits',
        default=5,
        type=int,
        help='Maximum number of hits to return for a search. Default: 5')
    blast_args.add_argument(
        '-b',
        '--database',
        type=str,
        help='Path to a local BLAST database.')
    blast_args.add_argument(
        '-q',
        '--entrez-query',
        help='Entez query to limit BLAST searches to specific taxon.')

    #   Finally, we specify general SNPMeta options. The user may specify an
    #   output file, target organism(s), a flag for whether or not annotation
    #   should proceed against XML BLAST reports, and whether ot not the output
    #   should be in dbSNP or tabular format.
    parser.add_argument(
        '-o',
        '--output',
        metavar='OUTPUT_FILE',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='File to write the annotations. Default: STDOUT')
    parser.add_argument(
        '-t',
        '--target-organism',
        metavar='ORGANISM',
        action='append',
        help='The oganism(s) against which to annotate. May be specified '
             'multiple times.')
    parser.add_argument(
        '--no-blast',
        action='store_true',
        default=False,
        help='BLAST has already been run; annotate from pre-existing XML '
             'reports. Default: False')
    parser.add_argument(
        '--outfmt',
        type=str,
        choices=['dbsnp', 'tabular'],
        default='dbsnp',
        help='The output format in which to write annotations. Default: dbSNP')

    #   Then parse the arguments
    parsed_args = parser.parse_args()
    return(parsed_args)
