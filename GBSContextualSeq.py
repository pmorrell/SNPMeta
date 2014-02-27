#!/usr/bin/env python
"""
A script to take VCF files and build sample-specific contextual sequences
for SNPs, in the format required by SNPmeta. Takes a VCF and a reference FASTA
as input, and writes multiple FASTA files, depending on number of samples. 
"""
###############################################################################
#   July 3, 2012
#   Thomas Kono
#   Saint Paul, MN
#   Dependencies
#       Biopython:  http://biopython.org/
#       Argparse:   http://code.google.com/p/argparse/
#           Standard in python 2.7 and 3.2
#
###############################################################################

#   Import modules
#       For command line arguments
import argparse
#       For regular expression matching, for detanlging the FORMAT field
import re
#       Biopython modules
#       This lets us easily read FASTA files.
from Bio import SeqIO

#   Lookup Tables
#       IUPAC two-base ambiguities
IUPAC = {   'M': ('A', 'C'),
            'R': ('A', 'G'),
            'W': ('A', 'T'),
            'S': ('C', 'G'),
            'Y': ('C', 'T'),
            'K': ('G', 'T'),
        }

#   Argument parsing
#       A description of the program
DESCR = """A script to take a VCF file and a reference genome, and produce 
sample-specific FASTA files for input to SNPmeta. Depending on the size of the
genome and the VCF file, this can take up a LOT of memory."""
#       Start a new parser object, using a generated usage message and descr
Arguments = argparse.ArgumentParser(description=DESCR, add_help=True)
#       Add arguments - we only take three
Arguments.add_argument('-r',
    '--ref',
    metavar='REF_SEQ',
    help='Reference sequence, as FASTA.')
Arguments.add_argument('-f',
    '--vcf',
    metavar='VCF',
    help='VCF file from which to build contextual sequence')
Arguments.add_argument('-l',
    '--len',
    metavar='CONTEXT_LENGTH',
    type=int,
    default=100,
    help='The length of the contextual sequence around the SNP (Default: 100)')
#       Parse the arguments
ParsedArgs = Arguments.parse_args()

#   Functions
def check_arguments(args):
    """A function to check the arguments passed to the script."""
    #   Convert the arguments to a dictionary
    argdict = args.__dict__
    #   If our required arguments are 'None', then we should print help
    if not argdict['ref']\
    and not argdict['vcf']:
        Arguments.print_help()
        exit(1)
    #   Now go through the arguments one at a time, and check for validity
    if not argdict['ref']:
        print('Please supply a reference sequence!')
        exit(1)
    if not argdict['vcf']:
        print('Please supply a VCF file!')
        exit(1)
    #   context length should be an integer greater than 0
    #   We will print a warning if it is too short, less than 25.
    if argdict['len'] <= 0:
        print('Contextual length should be a positive integer!')
        exit(1)
    elif argdict['len'] < 25:
        print('Warning: Contextual lenghth is less than 25! BLAST Results might come from wrong organism!')
    #   Make sure we can open the files
    try:
        handle = open(argdict['ref'], 'r')
    except IOError:
        print('Error: The reference file does not exist, or is not readable!')
        exit(1)
    try:
        handle = open(argdict['vcf'], 'r')
    except IOError:
        print('Error: The VCF file does not exist, or is not readable!')
        exit(1)

def read_reference(ref):
    """A function to read the reference sequence as FASTA and return the
    sequence and the record names as a dictionary."""
    #   Write a little message that we are parsing the reference
    print('Reading ' + ref + ' ...')
    #   Use the BioPython parser to turn the sequence into a list
    parsedref = list(SeqIO.parse(ref, 'fasta'))
    #   We will turn it into a dictionary of the form
    #   {   name: sequence, name: sequence .... }
    refdict = {}
    for record in parsedref:
        #   'name' holds the names - for a reference, the chromosomes
        #   'seq' holds the sequence
        refdict[record.name] = record.seq
    #   Return the dictionary
    print('Done!')
    return(refdict)

def read_vcf(vcf, chromosomes):
    """A function to read a VCF file and strip out the SNP states and
    their chromosomal positions. It also handles separate samples.
    Also performs a check to make sure that the """
    #   A message that we are going through the VCF
    print('Reading ' + vcf + ' ...')
    #   Start an empty list to hold the data
    parsed_vcf_data = []
    #   Open the VCF, and read through it
    with open(vcf, 'r') as f:
        for line in f:
            #   We want to skip the lines that start with ##
            if line.startswith('##'):
                continue
            #   VCF has eight fixed fields, the first of which is #CHROM
            #   These are _always_ in the same order
            #   This is where the data starts
            elif line.startswith('#CHROM'):
                #   We want to save the samples, which are listed after FORMAT
                #   First, break it up into a list
                data_header = line.split()
                #   Get the index of the 'FORMAT' field
                #   Sometimes, VCF files won't have this - this means there is
                #   no genotype data in the file, and we cannot use it
                if 'FORMAT' in data_header:
                    #   list.index(value) returns the index of value, but it
                    #   has to exist, else it raises an error, which is why we
                    #   put it into an if block
                    format_field = data_header.index('FORMAT')
                else:
                    #   If not, we want to die with a message
                    print('This VCF file does not have any genotype data!')
                    exit(1)
                #   Take everything after FORMAT to the end - the samples
                samples = data_header[format_field+1:]
            #   We also don't want to mess with indels, since those are tricky
            elif 'INDEL' in line:
                continue
            #   Information is now mined out of the header
            #   Time to start accumulating genotype data
            else:
                #   Break the line up by entries
                variant_data = line.split()
                #   We want to make sure that the chromosomes listed in the VCF
                #   are actually the same as those in the reference sequence.
                #   If not, then we don't a proper combination of files.
                if variant_data[0] not in chromosomes:
                    print('The chromosomes in the VCF do not match those in the reference!')
                    exit(1)
                #   Save the relevant parts of the file as a tuple
                #   We need the chromosome (#CHROM), the position (POS), the 
                #   reference allele (REF), the alternate allele (ALT), and any
                #   sample information. We don't have to save the FORMAT field
                #   since we are only interested in the genotype calls, and VCF
                #   specifications are such that the first thing that must be 
                #   listed in the sample fields is the genotype information.
                #   We want to manipulate the REF and ALT alleles a little
                #   since they are listed as 0, 1, 2, and so on.
                #   The REF allele is always listed as 0 in genotype calls,
                #   and the ALT allels are 1, 2, ... comma-delimited
                #   If they are part of a list, then we can just use the index
                #   to get them, and it is much easier. split() gives a list,
                #   and + is a list concatenation operator.
                alleles = variant_data[3].split() + variant_data[4].split(',')
                relevant_data = (variant_data[0],
                                variant_data[1],
                                alleles,
                                variant_data[format_field+1:]
                                )
                #   Append the chromosomes, if they are not already in the list
                #   Append it to the beginning list
                parsed_vcf_data.append(relevant_data)
    print('Done!')
    #   And return the sample names and the huge list
    return(samples, parsed_vcf_data)

def build_sequences(samples, vcf_data, context_length, reference):
    """A function to build the query sequences based on the calls and positions
    given in the VCF file."""
    #   Print how many files we will create
    print('Will write ' + str(len(samples)) + ' files, one for each sample.')
    #   First, we want to iterate over samples, since those are the main unit
    #   We need to enumerate because the genotype information is in a separate
    #   list - we need to fetch by index.
    for sample_number, sample_name in enumerate(samples):
        #   Generate a filename and open a write handle to it
        snp_filename = sample_name + '_SNPs.fasta'
        handle = open(snp_filename, 'w')
        #   Print a little message
        print('Opening ' + snp_filename + ' for writing...')
        #   Then iterate over vcf_data - these are SNP calls
        #   again, enumerate() so we can count
        for index, call in enumerate(vcf_data):
            #   This is the genotype call for this SNP
            #   the fourth item of each entry in vcf_data is the genotype data
            #   stored as a list. The genotype data is stored as a colon-
            #   delimited list, with the actual calls as the first entry,
            #   as per the VCF standard.
            genotype = call[3][sample_number].split(':')[0]
            #   We break it into a list with regular expressions
            #   findall() returns a list of the regex matches
            #   r'\d+' searches for digits
            #   We do it this way because the VCF spec says that the calls can
            #   be separated by '/' (unphased) or '|' (phased)
            #   We go from list to set, as we can use it like a unique()
            #   function. Cast back to list so we can iterate over it
            sample_alleles = list(set(re.findall(r'\d+', genotype)))
            #   If they contain the reference allele, then we skip them, 
            #   since they are not SNPs in this sample. We also want to 
            #   skip hets, since it is unclear how we should handle those.
            #   If they are hets, then there will be multiple genotype calls,
            #   and the length of the list will be different from 1.
            if '0' in sample_alleles or len(sample_alleles) != 1:
                continue
            else:
                #   The actual two states of the SNP
                ref_state = call[2][0]
                #   sample_alleles is a list, so we just take the first item
                #   if it made it here, then it should just be one item long
                alt_state = call[2][int(sample_alleles[0])]
                #   Get the IUPAC ambiguity, given these two states
                #   Use a list comprehension and then take the value out
                #   Ugly, but it works.....
                ambi = [base for base in IUPAC if ref_state in IUPAC[base] and alt_state in IUPAC[base]][0]
                #   The name of the SNP, formatted for FASTA
                #   >sample_xxxx
                snp_name = '>' + sample_name + '_' + str(index+1)
                #   The lower bound of the SNP contextual sequence
                #   call[1] is the chromosomal position of the SNP, starting 
                #   from 1. Subtract 1 because of the way python counts
                lower_bound = int(call[1]) - context_length - 1
                #   And the upper bound
                upper_bound = int(call[1]) + context_length
                #   Get the right slice out of the reference sequence
                #   call[0] is the chromosome name. Cast to list so that we can
                #   edit the values
                snp_sequence_lst = list(reference[call[0]][lower_bound:upper_bound])
                #   Place the ambiguity
                snp_sequence_lst[context_length] = ambi
                #   Put it back together as a string
                snp_sequence = ''.join(snp_sequence_lst)
                #   Now write it all to the file
                handle.write(snp_name + '\n')
                handle.write(snp_sequence + '\n')
        print('Done!')
        handle.close()

#   Now we do all the work here
#       Check the arguments
check_arguments(ParsedArgs)
#       Parse the reference
ref = read_reference(ParsedArgs.ref)
#       Get the data out of the VCF
samples, vcf_data = read_vcf(ParsedArgs.vcf, ref.keys())
#       Build the sequences and write them to files
build_sequences(samples, vcf_data, ParsedArgs.len, ref)
