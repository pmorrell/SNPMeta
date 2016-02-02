#!/usr/bin/env python
"""
Annotates SNPs using sequences from a specified organism in NCBI's nonredundant
nucleotide sequence database. Produces a dbSNP submission report, or a tab-
delimited report. Requires Biopython.
"""

#   Import standard library modules
import sys
import os
import datetime
import re
import math
import time
import subprocess
import socket

#   Import the Biopython libraries
try:
    #   To run BLAST on NCBI's servers
    from Bio.Blast import NCBIWWW
    #   To run BLAST locally
    from Bio.Blast.Applications import NcbiblastnCommandline
    #   To work with NCBI's XML reports
    from Bio.Blast import NCBIXML
    #   To use NCBI's Efetch interface for data retrieval
    from Bio import Entrez
    #   To parse genbank flat files
    from Bio import SeqIO
    #   To perform operations on sequence features/annotations
    from Bio import SeqFeature
    #   This requires EMBOSS to be installed, and performs pairwise alignment
    from Bio.Emboss.Applications import NeedleCommandline
    #   To parse the alignment that results from EMBOSS tools
    from Bio import AlignIO
    #   To do translations, and handle translation tables
    from Bio.Seq import Seq
    #   To handle sequence types (DNA sequence, as opposed to RNA or protein)
    from Bio.Alphabet.IUPAC import unambiguous_dna
    from Bio.Alphabet.IUPAC import ambiguous_dna
except ImportError:
    print('This script requires the BioPython library to be installed.')
    exit(1)

#   Import package components
from snpmeta.LookupTables.iupac import IUPAC
from snpmeta.LookupTables.grantham import GSCORES
from snpmeta.ArgumentHandling import parse_args
from snpmeta.ArgumentHandling import validate_args

###############################################################################
#       CONSTANTS AND STRINGS
###############################################################################
#   The email address and tool name to send to NCBI
Entrez.tool = 'SNPMeta, using BioPython'
Entrez.email = ParsedArgs.email

#   The header of the table output
TABLE_HEADER = [    'SNPName',
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
#   The name of the temporary file to hold alignments
NEEDLE_OUTPUT = 'needle_output.txt'
#   The name of the temporary file to hold the genbank sequence
GENBANK_SEQ = 'gbseq.fa'
#   This is taken directly from the argument list
TARGET_ORGANISMS = ParsedArgs.target_organism
#   A class to store the SNP annotation data while the script runs
class SNPAnnotation:
    """A class to store the SNP annotation data."""
    def __init__(self):
        self.SNPName = '-'
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
        self.DateTime = '-'

###############################################################################
#       FUNCTIONS
#   These are the functions that do the bulk of the work
###############################################################################
def get_chunks(sequence, size):
    """Returns a list of lists, effectively splits a big list into a bunch of
    smaller lists of a certain size. We have to define this fuction because
    sometimes the URLs get too long to send to NCBI"""
    return(sequence[p:p+size] for p in xrange(0, len(sequence), size))


def get_dir_contents(no_blast):
    """Returns the list of files to work on. Depends on if --no-blast was given
    on the command line."""
    #   If --no-blast was supplied...
    if no_blast:
        #   If --no-blast was specified, operate on files ending in 'xml'
        suffix = 'xml'
    else:
        #   If not, then we want the FASTA to BLAST them
        suffix = 'fasta'
    #   The full path of the directory supplied
    #   We need this since the script looks relative to where it is located
    #   on the disk.
    full_path = os.path.abspath('.')
    #   This is the full directory listing of the current directory, as a list
    file_list = os.listdir(full_path)
    #   Trim the file list by those ending in the proper suffix
    file_list = [x for x in file_list if x.endswith(suffix)]
    return(file_list)

def calculate_context(seq, directory, context_len, gbs, illumina):
    """Calculates the contextual sequence length of the sequence passed to it, 
    given that the sequence came from either GBS or Illumina design."""
    if directory:
        #   Read the correct FASTA off disk
        query_seq = SeqIO.read(seq, 'fasta')
    #   Else, we have been served a record already
    else:
        query_seq = seq
    #   If we are not in the illumina format, then we can just return the 
    #   length that was provided
    if not illumina:
        if context_len:
            if len(seq) > context_len:
                if seq[context_len] in IUPAC:
                    return(query_seq, context_len)
                elif seq[-(context_len+1)] in IUPAC:
                    return(query_seq, len(query_seq)-(context_len+1))
                else:
                    for index, base in enumerate(query_seq):
                        if base in IUPAC:
                            return(query_seq, index)
            else:
                for index, base in enumerate(query_seq):
                    if base in IUPAC:
                        return(query_seq, index)
        #   If we are working with GBS, then we take the first SNP found
        if gbs:
            for index, base in enumerate(query_seq):
                if base in IUPAC:
                    return(query_seq, index)
    #   We are in illumina format
    else:
        #   Just search for the start, since this is the same position where
        #   the new IUPAC ambiguity will reside
        offset = re.search('\[', str(query_seq.seq)).start()
        #   Then build the new sequence.
        #   We split on [ / ], and the SNP states will always be the middle
        #   two elements
        seq_parts = re.split('[\[|\/|\]]', str(query_seq.seq))
        #   We replace it with the ambiguity
        for amb, stdnuc in IUPAC.items():
            if (seq_parts[1] in stdnuc) and (seq_parts[2] in stdnuc):
                snp = amb
                break
        query_seq.seq = seq_parts[0] + snp + seq_parts[-1]
        #   and now, we return it
        return(query_seq, offset)

def setup_env(arguments):
    """Sets up the 'environment' of the program - gets a list of SNPs to annotate,
    either from files, or from records stored in one file. Also sets up printing
    a table or a dbSNP-compatible report."""
    #   Are we annotating a directory or a single file?
    if arguments.directory:
        #   Change directories into the target
        os.chdir(arguments.directory)
        #   Get the list of files we are going to annotate
        record_list = get_dir_contents(arguments.no_blast)
        #   Drop a message if the directory is empty
        #   An empty list evaluates to Boolean False
        if not record_list:
            print('No SNPs in this directory!')
            exit(1)
    elif arguments.fasta_file:
        #   Get the generator object for the records in the FASTA file
        #   We use BioPython's SeqIO.parse() function, which returns
        #   a generator that spits out one FASTA record at a time
        record_list = SeqIO.parse(arguments.fasta_file, 'fasta')
    #   This section only runs if the user supplied the -v flag
    if arguments.verbose:
        #   Write the table header
        arguments.output.write('\t'.join(TABLE_HEADER))
        #   Add a newline
        arguments.output.write('\n')
    #   Send back the list of records to be annotated
    return(record_list)

def run_blast(sequence, program, evalue, max_return, database, equery):
    """Runs BLAST with the parameters given on the command line. Returns a
    handle to the BLAST results, in XML format."""
    #   Check the arguments - are we reading a file, or are we reading a SeqRecord?
    #   Since -f is mutually exclusive with -d, then we only need to check one
    #   Variable 'directory' will be True if we are reading a directory of SNPs
#    if directory:
        #   Read the sequence up as a string
        #       First, a handle to the FASTA with the sequence
#        handle = open(sequence, 'r')
        #       Read and parse the sequence
#        query = SeqIO.read(handle, 'fasta')
#        handle.close()
    #   Else, we are reading a large FASTA, and are being served one record at a time
#    else:
#        query = sequence
    #   Print a message to stderr, to keep track of what's going on
    sys.stderr.write('Running '+program+' on '+sequence.id+'...')
    #   Convert Illumina format to regular sequence format
    #   We split on [ / ] to make a list with four elements
    #   The middle two are always the query SNP
#    if illumina:
#        seq_parts = re.split('[\[|\/|\]]', str(query.seq))
#        #   We replace it with the ambiguity
#        for amb, stdnuc in IUPAC.items():
#            if (seq_parts[1] in stdnuc) and (seq_parts[2] in stdnuc):
#                snp = amb
#                break
#        query.seq = seq_parts[0] + snp + seq_parts[-1]
    #   If '-b' was supplied, we are running BLAST locally...
    if database:
        #   We have to write a temporary file, as local BLAST needs one
        SeqIO.write(sequence, 'tmp_blast.fasta', 'fasta')
        blast_command = NcbiblastnCommandline(query='tmp_blast.fasta',
                                            out='tmp_blast.xml',
                                            db=database,
                                            evalue=evalue,
                                            outfmt=5,
                                            num_descriptions=max_return,
                                            num_alignments=max_return)
        #   Run BLAST
        subprocess.call(str(blast_command).split())
        blast_handle = open('tmp_blast.xml', 'r')
    else:
        #   Also have to wrap the BLAST command into a try-except
        success = False
        while not success:
            try:
                #   If there is an Entrez query, we have to build it into the query
                if equery:
                    blast_handle = NCBIWWW.qblast(program,
                                'nr',
                                sequence.seq,
                                format_type='XML',
                                descriptions=max_return,
                                alignments=max_return,
                                hitlist_size=max_return,
                                expect=evalue,
                                entrez_query=equery)
                #   If there isn't an Extrez query, then we run BLAST just like this
                else:
                #   And now run the search
                    blast_handle = NCBIWWW.qblast(program,
                                'nr',
                                sequence.seq,
                                format_type='XML',
                                descriptions=max_return,
                                alignments=max_return,
                                hitlist_size=max_return,
                                expect=evalue)
                success = True
            except urllib2.HTTPError as e:
                sys.stderr.write('Caught' + str(e.code) + ' - retrying in 5 seconds...\n')
                time.sleep(5)
    #   Tell the user it is finished
    sys.stderr.write('Done!\n')
    #   Return the handle
    return(blast_handle)

def get_gbnumbers(no_blast, xmlfile):
    """Reads the XML report, extracts query name, Genbank ID, directions.
    Returns a tuple with this information for each hit"""
    #   If --no-blast was specified, then we have to read the XML report
    #   that is already on the disk
    if no_blast:
        blast_output = open(xmlfile, 'r')
    #   If not, then we already have a handle
    else:
        blast_output = xmlfile
    #   Start a parser that steps through each record
    blast_records = NCBIXML.parse(blast_output)
    #   List to hold information about our hits
    hits = []
    #   Step through the BLAST records
    for record in blast_records:
        #   Step through each alignment in each record
        for alignment in record.alignments:
            #   Then the HSPs in each alignment
            for hsp in alignment.hsps:
                #   The start and end positions of each hit
                hit_coords = (hsp.sbjct_start, hsp.sbjct_end)
                #   Split on the '|' character, genbank ID is last in the list
                #   have to use -2 instead, because of the trailing '|' in the
                #   XML report
                hit_gbid = alignment.title.split('|')[-2]
                #   Relative directions of the sequences
                hit_directions = hsp.frame
                break
                #   Put it all into a tuple and tack it onto hits.
                #   They _should_ be sorted from low to high by e-value
            hits.append((hit_gbid, hit_coords, hit_directions))
    #   Finished with this file
    blast_output.close()
    #   Return a tuple of the query name and the hits
    return(hits)

def extract_gbinfo(hit_info):
    """Fetch the GenBank records corresponding to each GenBank ID from
    the BLAST results. Return a tuple containing the record for annotation
    and a record for a related gene in a related organism."""
    #   Lists to hold the unpacked information
    genbank_ids = []
    coordinates = []
    directions = []
    #   Unpack all the information into the lists
    for hit in hit_info:
        genbank_ids.append(hit[0])
        coordinates.append(hit[1])
        directions.append(hit[2])
    #   Two lists to hold start positions and end positions
    starts = []
    ends = []
    #   Use the list of coordinates to build request ranges for each sequence
    for each in coordinates:
        #   We will request ~5kb total, 2.5kb upstream and 2.5k downstream
        #   5kb upstream
        hitstart = min(each) - 2500
        #   5kb downstream
        hitend = max(each) + 2500
        #   Weird things happen if a start request is negative,
        #   so set it to 1 if it is negative
        if hitstart < 0:
            hitstart = 1
        #   Store as strings, since they have to be sent to NCBI as strings
        starts.append(str(hitstart))
        ends.append(str(hitend))
    #   This is an unideal workaround, but it should work for now
    #   If there are more than 20 records to parse, then we could 
    #   Build the query string for NCBI
    #       We get data from the 'nucleotide' database
    query_db_name = 'nucleotide'
    #       We want data returns as genbank flat file
    query_ret_type = 'gb'
    #       And we want it returned as plain text
    query_ret_mode = 'text'
    #   Get the chunks that we will query NCBI with
    #   We choose 20 at a time, that should be reasonable
    genbank_chunks = get_chunks(genbank_ids, 20)
    start_chunks = get_chunks(starts, 20)
    end_chunks = get_chunks(ends, 20)
    #   Actually query NCBI!
    #   Entrez.efetch() forms the proper query and sends it to NCBI's servers
    #   Returns handle to data, which is contatenated genbank records
    #   We wrap this in a try-except block to catch some HTTP errors that
    #   NCBI gives during peak hours. Hopefully this helps
    #   We also give a counter, just in case. We don't want to get caught in 
    #   an infinite loop
    #   Keep an empty string to append the records to and parse it later
    records = []
    for chunk in zip(genbank_chunks, start_chunks, end_chunks):
        attempts = 0
        success = False
        while (not success) or (attempts < 5):
        #   If we have more than five attempts, then we quit
            if attempts > 5:
            #   Return a tuple of None if we weren't able
            #   to get the data from the servers
                return(None, None, None)
            try:
                handle = Entrez.efetch(db=query_db_name, 
                            id=','.join(chunk[0]),
                            rettype=query_ret_type,
                            retmode=query_ret_mode,
                            seq_start=','.join(chunk[1]),
                            seq_stop=','.join(chunk[2]))
                temp_genbank = SeqIO.parse(handle, 'genbank')
                records = records + list(temp_genbank)
                #   wait for a second, to avoid abusing the NCBI servers
                time.sleep(1)
                handle.close()
                success = True
            except urllib2.HTTPError as e:
                message = 'Caught ' + str(e.code) + ' - retrying in 5 seconds...\n'
                sys.stderr.write(message)
                time.sleep(5)
            attempts += 1
    #   The parser provied handles this
    #   but we must specify that it is reading a genbank file
    #records = SeqIO.parse(handle, 'genbank')
    #   We zip it with the hit directions for ease
    #   Put it into a list because SeqIO.parse() returns a generator object
    #   which can only be stepped through once. By putting it into a list,
    #   we can traverse it multiple times.
    dirs_recs = zip(directions, [record for record in records])
    #   An empty tuple for annotated regions, for if there are no annotations
    annotated_regions = ()
    annotated_genes = ()
    #   An empty tuple for related organisms and loci, 
    #   for if there is not a related named gene
    related_locus = ()
    #   Step through our list for the first time...
    #   On this pass, we will only get hits from our target organisms...
    #   And out of those, only the hits with an annotated CDS
    for direction, record in dirs_recs:
        #   We want to extract any genes that may be on the record
        annotated_genes = [each for each in record.features if each.type == 'gene']
        #   Extract the organism from the info
        organism = record.annotations['organism']
        #   If the user has supplied some target organisms...
        if TARGET_ORGANISMS:
            #   Skip uninteresting organisms
            if organism in TARGET_ORGANISMS:
                #   A list of all features in the record
                types = [feature.type for feature in record.features]
                #   We only care about those that have a CDS annotation
                #   Again, these should be sorted from low to high, by e-value
                if 'CDS' in types:
                    #   Save the name, organism of origin, the features, 
                    #   and the direction it hit the genbank record
                    annotated_regions = (record.name, organism, record.features, direction)
                    #   And save the genbank sequence for later
                    #   Recall that 'GENBANK_SEQ' was defined in 'CONSTANTS'
                    SeqIO.write(record, GENBANK_SEQ, 'fasta')
                    #   We popped the first good one off, so we break
                    break
        #   If not, then we just step through it, taking the best annotated hit
        #   This is the same as the above block, without the check for membership
        else:
            types = [feature.type for feature in record.features]
            if 'CDS' in types:
                annotated_regions = (record.name, organism, record.features, direction)
                SeqIO.write(record, GENBANK_SEQ, 'fasta')
                break
    #   The second pass through the list of hits, 
    #   this time just finding another gene annotation, 
    #   in case it is useful for something
    for direction, record in dirs_recs:
        #   Extract the organism in the same way as above
        organism = record.annotations['organism']
        #   Now, we are interested in all the features
        for feature in record.features:
            #   If there is a gene, then we will save it
            if feature.type == 'gene':
                #   Pull the name of the gene from the feature
                genename = feature.qualifiers.get('gene', '-')[0]
                #   Build a tuple of the related organism information
                related_locus = (genename, organism)
                #   And again, break once we find one
                break
    return(annotated_regions, related_locus, annotated_genes)

def sequence_align(seq, query, frames):
    """Uses the 'needle' program in the EMBOSS package to align the query
    sequence with the GenBank sequence. Takes a tuple of the hit
    directions to determine the proper orientation for alignment."""
    #   We need to check: are we reading a directory of FASTA files, or are we
    #   Reading a single file with many records?
    #   -f and -d are mutually exclusive, so we need to check just one
    #   We check the directory argument, should be True if reading a directory
#    if directory:
        #   Read the correct FASTA off disk
#        query_seq = SeqIO.read(query, 'fasta')
    #   Else, we have been served a record already
#    else:
#        query_seq = query
    #   Check the frames. If they are the same, then we can align as-is
    #   If not, then we we need to RC the query sequence
    #   We RC the query, since it would be complicated and 
    #   potentially introduce error to RC the genbank sequence
    if ('-' in str(frames[0]) and '-' not in str(frames[1])) or ('-' in str(frames[1]) and '-' not in str(frames[0])):
        snp_seq = Seq(query.seq)
    else:
        snp_seq = Seq(query.seq).reverse_complement()
    #   These are the defaults from the 'Geneious' aligner
    gapopen = 12
    gapextend = 3
    #   Specifying the options for the 'needle' program
    needle_cmd = NeedleCommandline()
    #   The first sequence in the genbank sequence
    needle_cmd.asequence = seq
    #   The second is the raw sequence from the FASTA
    #   it is raw, since we might have had to RC it
    needle_cmd.bsequence = 'asis:'+str(snp_seq)
    #   Pass the other options to the aligner
    needle_cmd.gapopen = gapopen
    needle_cmd.gapextend = gapextend
    #   The output file. This was defined in the 'CONSTANTS' section
    needle_cmd.outfile = NEEDLE_OUTPUT
    #   Run it!
    needle_cmd()

def get_snp_position(alignment, sequence, clength):
    """Get the position of the SNP in the alignment, based on the sepcified
    contextual sequence length. Returns the offset, the code for the SNP, 
    and any warnings about the alignment. Also extract the upstream and 
    downstream sequences, for dbSNP report generation."""
    #   Define ambiguous nucleotides to search for
    amb_nuc = ['R', 'K', 'M', 'S', 'Y', 'W']
    #   We want to get the score of the alignment
    #   The format of the line is 
    #   # score: xxx.x
    with open(alignment, 'r') as f:
        for line in f:
            if line.startswith('# Score:'):
                score = float(line.split(':')[-1].strip())
                break
    #   Read in the alignment, in 'emboss' format
    temp_align = AlignIO.read(alignment, 'emboss')
    #   if we are reading sequence out of a directory, then we have to actually
    #   hit the disk again
#    if directory:
#        queryseq = SeqIO.read(sequence, 'fasta')
#    else:
#        queryseq = sequence
    #   We then pass this off to a helper function to calculate the
    #   contextual sequence length
    #   If context_len is undefined, then we know that either GBS or Illumina
    #   format has been provided. 
#    if not context_len:
#        clength = calculate_context(queryseq, gbs, illumina)
#    else:
#        clength = context_len
    #   First sequence is the genbank sequence, maintain as a 'Seq' object
    #   so that we can extract regions with a 'Seqfeature' for later
    genbank_seq = temp_align[0].seq
    #   Second sequence is the SNP sequence
    #   Cast it to a string, to use regular expressions on it
    snp_seq = str(temp_align[1].seq)
    #   Match the first instance of a non-gap character
    snp_start = re.search('[^-]', snp_seq).start()
    #   Match the last instance of a non-gap character
    #   Two situations are:
    #       NNNNNNNNNN------$
    #   or  NNNNNNNNNNNNNNNN$
    #   Add one for the way python indexes
    search_pattern = '([ATCGRYMWKSN]-+$)|([ATGCRYMKWSN]$)'
    #   We want to ignore the case of letters
    snp_end = re.search(search_pattern, snp_seq, re.IGNORECASE).start() + 1
    #   The complete range the aligned SNP sequence convers
    snp_coverage = range(snp_start, snp_end)
    #   Extract the relevant portion of the SNP sequence
    aligned_snp = snp_seq[snp_start:snp_end]
    #   These are the cases that are possible in the alignment:
    #   AGGATG-(*k)TTAGAGAT... [SNP] AGGAGTCGATG... | SNP with gap on left
    #   GATGGATGCCGTG... [SNP] AGCCTCAT(-*k)GATG... | SNP with gap on right
    #   GAAC(-*c)GACC... [SNP] AGTAGC(-*k)GGATGC... | SNP with gaps on both
    #   We remove the gaps from the sequence, and also remove those positions
    #   from the mapping positions
    coverage_pos = zip(snp_coverage, aligned_snp)
    actual_positions = [base[0] for base in coverage_pos if base[1] != '-']
    #    We set these variables to blank, because they are not relevant
    #    in the non-GBS context
    fiveprime_seq = '-'
    threeprime_seq = '-'
    observed = '-'
    #   Count the number of expected bases in from the end to get the position
    #   of the SNP in the alignment
    #   Counting from the left
    fiveprime_offset = actual_positions[clength]
    #   And counting from the right
    threeprime_offset = actual_positions[-(clength+1)]
    #   And get the bases at these positions
    fiveprime_base = snp_seq[fiveprime_offset]
    threeprime_base = snp_seq[threeprime_offset]
    #   If the contextual sequence is just as expected
    #   The base grabbed should be an ambiguous base, and should be the same
    #   no matter which side it was grabbed from
    if fiveprime_base in amb_nuc and threeprime_base == fiveprime_base:
        #   There is no message about if this sequecne is badly formed
        badlength = '-'
        #   Store the relevant information
        amb_base = fiveprime_base
        offset = fiveprime_offset
    #   If not, then we have to do some searching
    #   We check if the 5' side is the right length
    elif fiveprime_base in amb_nuc and fiveprime_base != threeprime_base:
        badlength = 'Short on 3\' side!'
        amb_base = fiveprime_base
        offset = fiveprime_offset
    #   And then we check the 3' side
    elif threeprime_base in amb_nuc and fiveprime_base != threeprime_base:
        badlength = 'Short on 5\' side!'
        amb_base = threeprime_base
        offset = threeprime_offset
    #   The last resort is to find any ambiguity, starting from the 5' side
    #   These are prone to errors in annotation, but
    #   I don't know of any other way to handle this
    else:
        badlength = 'Short on both sides! Possibly an error!'
        for index, base in enumerate(snp_seq):
            if base in amb_nuc:
                amb_base = base
                offset = index
                #   Stop once we find one
                break
    #   If we are annotating from GBS, then we just search through the sequence
    #   We repeat the block, with a twist, for if our sequence is 'malformed', but it is supposed
    #   to look this way, so we don't give a warning
    #   We just have to be careful of sequences with multiple ambiguities
#    else:
#        #   Use lists to store the ambiguous bases
#        ambs = []
#        offsets = []
#        for index, base in enumerate(snp_seq):
#            if base in amb_nuc:
#                ambs.append(base)
#                offsets.append(index)
#        #   If there is more than 1 ambiguous base, then we issue a warning
#        if len(ambs) > 1:
#            badlength = 'Contextual sequence has multiple ambiguities!'
#        else:
#            badlength = '-'
#        #   And annotate off the first one encountered.
#        amb_base = ambs[0]
#        offset = offsets[0]
    #   Zip the SNP sequence and the GenBank sequence
    #   Step though it, and save the first ambiguity we find in the SNP seq
    #   Take the corresponding base from the GB seq
    #   Build the 'observed' field
    fiveprime_seq = str(sequence.seq[:clength])
    threeprime_seq = str(sequence.seq[clength+1:])
    observed = '/'.join(IUPAC[sequence[clength]])
#    combined = zip(genbank_seq, snp_seq)
#    fiveprime_seq = []
#    threeprime_seq = []
#    Found = False
#    for gbbase, snpbase in combined:
#        if snpbase in amb_nuc:
#            #   This is ugly, but it is the way to get the alternate state
#            #   Pull out the bases that the ambiguity covers, and get the one
#            #   that isn't the same as the genbank base
#            alt_base = [x for x in IUPAC[snpbase] if x != gbbase][0]
#            observed = gbbase + '/' + alt_base
#            Found = True
#        if not Found and snpbase != '-':
#            fiveprime_seq.append(snpbase)
#        elif Found and snpbase != '-' and snpbase not in amb_nuc:
#            threeprime_seq.append(snpbase)
#    #   Remove all the gaps from the SNP sequence
#    snp_seq = ''.join([base for base in snp_seq if base != '-'])
#    #   Use the offset provided to get the 5' and 3' flanking sequences
#    fiveprime_seq = ''.join(fiveprime_seq)
#    threeprime_seq = ''.join(threeprime_seq)
#    #   Return all the information
    return(amb_base, offset, genbank_seq, badlength, fiveprime_seq, threeprime_seq, observed, score)

def codon_position_trans(ambiguity, position, ref_seq, features):
    """Check if the SNP is coding or not, and return amino acids accordingly.
    """
    #   Raw reference sequence that we wrote to disk right after XML parsing
    handle = open(GENBANK_SEQ, 'r')
    nogap_ref = SeqIO.read(handle, 'fasta')
    handle.close()
    #   Step through the list of features in the genbank record
    for feature in features:
        #   We can only annnotate off CDS
        if feature.type == 'CDS':
            #   The positions of the reference sequence in the alignment, skipping gaps
            nogap_ref_positions = [index for index, base in enumerate(ref_seq) if base != '-']
            #   All positions in the alignment where the feature exists
            total_feature_range = [pos for index, pos in enumerate(nogap_ref_positions) if index in feature]
            feat_start = total_feature_range[0]
            feat_end = total_feature_range[-1]
            #   Check if our position is in the ungapped reference
            if position not in nogap_ref_positions:
                warning = 'SNP aligns to a gap in the GenBank sequence'
                coding = False
                in_feature = False
                snp_threeutr = '-'
                snp_fiveutr = '-'
                break
            #   Check if it is coding or not
            elif position in total_feature_range:
                coding_sequence = feature
                in_feature = True
                coding = True
                snp_threeutr = '-'
                snp_fiveutr = '-'
                break
            #   If the SNP position is not in the feature
            elif position not in total_feature_range:
                in_feature = False
                coding = False
                threeprime_dist = []
                fiveprime_dist = []
                #   Calculate a bunch of distances, to the end and start of the CDS
                for pos in total_feature_range:
                    threeprime_dist.append(math.fabs(position - pos))
                    fiveprime_dist.append(math.fabs(pos - position))
                #   Check which distance is closer, so we can get a number for the distance
                #   To the nearest CDS
                if min(threeprime_dist) <= min(fiveprime_dist):
                    snp_threeutr = int(min(threeprime_dist))
                    snp_fiveutr = '-'
                else:
                    snp_fiveutr = int(min(fiveprime_dist))
                    snp_threeutr = '-'
    #   If the position is in the CDS
    if in_feature:
        #   some CDS do not have a 'translation' qualifier, so we have to translate them beforehand
        fallback_translation = coding_sequence.extract(nogap_ref.seq).translate()
        #   The number of amino acids the CDS encodes (including STOP)
        trans_feat_len = len(coding_sequence.qualifiers.get('translation', fallback_translation)[0]) + 1
        #   The positions in the aligned reference sequence where there is CDS
        feat_positions = [base for index, base in enumerate(nogap_ref_positions) if index in coding_sequence]
        #   The start of the feature
        feat_start = feat_positions[0]
        #   If, for some reason, the position is not in the ungapped feature...
        if position not in feat_positions:
            #   Fill out blank or None for information that would have been collected later
            coding = False
            warning = 'SNP aligns to gap in CDS'
            snp_aa_pos = '-'
            snp_relpos = '-'
            snp_threeutr = '-'
            snp_fiveutr = '-'
            translations = None
            ordered_resolved = None
        #   If this character is in the location string, then the CDS starts
        #   before the position where the annotation starts, so we don't know the frame exactly...
        elif '<' in str(coding_sequence.location.start):
            warning = 'CDS annotation has a fuzzy start'
            snp_pos = feat_positions.index(position) - len(feat_positions)
        else:
            #   If the CDS starts on a normal position, then we are all-good
            snp_pos = feat_positions.index(position)
        #   Now, we only do the rest if we have a valid SNP alignment
        #   That is, if 'coding' is still True
        if coding:
            #   We want to update the translation table based on which we are supposed to use
            #   depending on the one annotated in the GenBank record
            #   We default to '1', the standard translation table
            table_id = coding_sequence.qualifiers.get('transl_table', '1')[0]
            #   If it is coding, then it will not be upstream nor downstream of the CDS
            snp_threeutr = '-'
            snp_fiveutr = '-'
            #   The position of the affected codon
            #   This can be negative because sometimes we count from the end
            #   That is okay, we can handle that
            #   Round up because any remainder is the next amino acid
            #   ie. quotient + 1
            #   If the position is positive, then we use math.ceil to round up
            if snp_pos >= 0:
                snp_aa_pos = int(math.ceil((snp_pos + 1)/3.0))
            else:
                snp_aa_pos = trans_feat_len + int(math.floor(snp_pos/3.0)) + 1
#            if snp_aa_pos < 0:
#                #   We take it from the end. But it is already negative, so we just add it normally
#                snp_aa_pos = trans_feat_len + snp_aa_pos + 1
                #   If it is positive, we use it as-is!
            #   Get the position in a codon
            #   If pos % 3 = 0, then we have first position
            #   If pos % 3 = 1, then we have second position
            #   If pos % 3 = 2, then we have third position
            snp_relpos = snp_pos % 3
            #   Get the codon based on the relative position
            #   A codon is given by
            #       [pos - (pos%3), pos + (2 - (pos%3))+1]
            #   If SNP is in 4th CDS position (index 3 in python)...
            #       [3-0, 3+(2-0)+1] -> [3, 6] -> [3, 4, 5] -> [4, 5, 6] (real counting)
            #   If SNP is in 5th CDS position (index 4 in python)...
            #       [4-1, 4+(2-1)+1] -> [3, 6] -> [3, 4, 5] -> [4, 5, 6] (real counting)
            #   If SNP is in 6th CDS position (index 5 in python)...
            #       [5-2, 5+(2-2)+1] -> [3, 6] -> [3, 4, 5] -> [4, 5, 6] (real counting)
            #   7th CDS position (6 in python)...
            #       [6-0, 6+(2-0)+1] -> [6, 9] -> [6, 7, 8] -> [7, 8, 9] (real counting)
            #   And so on...
            lower_bound = snp_pos - snp_relpos
            upper_bound = snp_pos + (2-snp_relpos) + 1
            #   The positions of the SNP-containing codon
            #   If the upper bound is 0, then we have to something different
            #   [-x:0] doesn't slice properly (gives empty list), so we just
            #   take until the end of the string.
            #   [-x:-y] does slice properly
            if upper_bound == 0:
                codon_positions = feat_positions[lower_bound:]
            else:
                codon_positions = feat_positions[lower_bound:upper_bound]
            #   Build the reference codon from these positions
            gb_codon = [ref_seq[pos] for pos in codon_positions]
            #   Check the strand of the feature
            #   If it is -1, then we have to RC the codon, and complement the ambiguity
            if coding_sequence.strand == -1:
                gb_codon = list(Seq(''.join(gb_codon), unambiguous_dna).reverse_complement())
                ambiguity = str(Seq(ambiguity, ambiguous_dna).complement())
                #   We also have to change the codon positions
                #   Third position becomes first, and first becomes third
                if snp_relpos == 2:
                    snp_relpos = 0
                elif snp_relpos == 0:
                    snp_relpos = 2
                #   The amino acid position has to be recalculated
                #   But it is just measuring from a different end
                snp_aa_pos = trans_feat_len - snp_aa_pos
            #   Test to see if the ambiguity represents a valid base in the genbank sequence
            #   That is, if the base in the genbank sequence is coded for in the ambiguity
            #   We have issues if it is not
            #   Recall that 'IUPAC' is defined in the 'LOOKUP' section
            if gb_codon[snp_relpos] not in IUPAC[ambiguity]:
                warning = 'SNP - Sequence Mismatch'
                #   We can't translate, nor get the codon states
                translations = None
                ordered_resolved = None
            #   Else, we are good to continue
            else:
                #   We made it this far - there should be no warnings from here on
                warning = '-'
                #   Build a copy for the SNP codon
                snp_codon = gb_codon[:]
                #   And drop the ambiguity into the right place
                snp_codon[snp_relpos] = ambiguity
                #   Translate the codon
                #   Can shortcut, since there should only be one ambiguity
                #   First, convert into a string
                gb_codon_str = ''.join(gb_codon)
                snp_codon_str = ''.join(snp_codon)
                #   Start a list to hold our unambiguous codons
                resolved_codons = []
                #   Next, resolve the ambiguities
                for stdnuc in IUPAC[ambiguity]:
                    #   Replace the ambiguity with a standard nucleotide
                    resolved_codons.append(snp_codon_str.replace(ambiguity, stdnuc))
                #   Which codon is the non-genbank SNP?
                non_gb_codon = [codon for codon in resolved_codons if codon != gb_codon_str][0]
                #   Now, build a tuple in the correct order
                #   The first item is the GenBank record codon
                #   The second is the mutated state
                #   We turn them into Seq objects for translation through BioPython's translation function
                #   A Seq object is a string with an attached 'alphabet', like DAN, RNA, Protein etc.
                #   We use the 'unambiguous_dna' alphabet, as ambiguities should have been removed
                ordered_resolved = (gb_codon_str, non_gb_codon)
                #   Then, we translate them
                translations = []
                for unambig_codon in ordered_resolved:
                    #   Use a built-in str() function to just get the string of the sequence
                    #   Pass the translation table ID to the translate() function
                    translations.append(str(Seq(unambig_codon, unambiguous_dna).translate(table=table_id)))
    #   If the SNP is NOT in the feature
    #   There is no information about the translations, the AA position, or the codons
    else:
        warning = '-'
        snp_aa_pos = '-'
        snp_relpos = '-'
        translations = None
        ordered_resolved = None
    #   Now we can return everything
    return(coding, warning, snp_aa_pos, snp_relpos, snp_threeutr, snp_fiveutr, translations, ordered_resolved)

###############################################################################
#       WORK
#   Now that we have defined all the functions, it is time to start annotating
###############################################################################
#   Do the simple argument check
check_arguments(ParsedArgs)
#   And set up the environment
file_list = setup_env(ParsedArgs)
#   Step through the list of files to operate on
for target in file_list:
    #   Start a new SNPAnnotation Class
    current_snp = SNPAnnotation()
    #   Current Date and Time
    current_snp.DateTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    #   Clear the contents of the temporary files, just in case
    #   Opening a file for writing clears its contents
    f = open(NEEDLE_OUTPUT, 'w')
    f.close()
    f = open(GENBANK_SEQ, 'w')
    f.close()
    #   If the directory was suppled, then we have to read the file to get the name
    if ParsedArgs.directory:
        #   If --no-blast was also specified, then we have to point to the FASTA file
        if ParsedArgs.no_blast:
            filename = target.replace('.xml', '.fasta')
        #   Now, we have to open the FASTA file
        #   This hits the disk to read the file, but there's not much else we can do
        else:
            filename = target
        current_snp.SNPName = SeqIO.read(filename, 'fasta').id
    #   If not, then we have already read the sequence.
    #   We just need to get its name
    else:
        current_snp.SNPName = target.id
        filename = target
    #   Were we provided with a length? Is it GBS? Illumina format?
    context = calculate_context(filename, ParsedArgs.directory, ParsedArgs.context_len, ParsedArgs.gbs, ParsedArgs.illumina)
    targetseq = context[0]
    context_len = context[1]
    #   Get the length of the sequence, for quality score purposes
    query_length = len(targetseq)
    #   Are we running BLAST?
    if ParsedArgs.no_blast:
        #   We are not, so we just read in the XML directly
        xml_report = target
    else:
        #   We are, so we need to BLAST it against nr, and return a handle to the results
        xml_report = run_blast(target, ParsedArgs.program, ParsedArgs.evalue, ParsedArgs.max_hits, ParsedArgs.database, ParsedArgs.entrez_query)
    #   Get the GenBank numbers out of the BLAST report
    blast_results = get_gbnumbers(ParsedArgs.no_blast, xml_report)
    #   If it is an empty tuple, then there are no BLAST hits
    if not blast_results:
        current_snp.Notes = 'No BLAST hits'
    #   If not, then we can continue
    else:
        #   Unpack the BLAST information
        annotations, related_locus, geneshortnames = extract_gbinfo(blast_results)
        #   If there is a related gene or organism, record that
        if related_locus:
            current_snp.AltGene, current_snp.AltOrg = related_locus
        #   If there are no annotations, take note. This is an empty tuple,
        #   which means we were able to download sequences but couldn't find
        #   any annotations
        if annotations == ():
            current_snp.Notes = 'No Annotations'
        #   If annotations is NoneType, then we weren't able to even fetch any
        #   sequences, so we print a message and move on
        elif annotations == None:
            current_snp.Notes = 'Could not fetch GenBank records'
        #   If there are annotations
        else:
            #   Store all the information from the annotations tuple
            current_snp.GenBankID, current_snp.Organism, gb_features, directions = annotations
            #   Use the 'needle' EMBOSS program to build an alignment
            sequence_align(GENBANK_SEQ, targetseq, directions)
            #   Get the code for the SNP, position in the alignment, and genbank seq
            snp_location = get_snp_position(NEEDLE_OUTPUT, targetseq, context_len)
            #   And unpack the information into the relevant fields
            current_snp.Ambiguity = snp_location[0]
            current_snp.ContextSeq = snp_location[3]
            current_snp.FiveFlank = snp_location[4]
            current_snp.ThreeFlank = snp_location[5]
            current_snp.Observed = snp_location[6]
            #   Calculate the score per seq. length, rounded to two decimal places
            current_snp.AlignScore = str(round(snp_location[7]/query_length, 2))
            #   Characterize the SNP according to coding sequence and amino acid states
            #   Send every piece of data from snp_location to the next function, up to the last four, in the same order
            snp_character = codon_position_trans(*snp_location[:-5], features=gb_features)
            #   And extract the information from the result
            #   There is a lot of data here
            coding, current_snp.Notes, current_snp.CDSPos, snppos, current_snp.ThreeUTR, current_snp.FiveUTR, aminoacids, codons = snp_character
            #   Search for some interesting features
            for each in gb_features:
                #   Genes are interesting
                if each.type == 'gene':
                    #   Get the short name of the gene.
                    #   For some reason, it is stored as a list
                    #   with only one element
                    current_snp.GeneShortName = each.qualifiers.get('gene', '-')[0]
                #   CDS with products are also interesting
                if each.type == 'CDS':
                    current_snp.Product = each.qualifiers.get('product', '-')[0]
                    #   Check to see if the CDS has a GenBank Protein ID
                    current_snp.ProteinID = each.qualifiers.get('protein_id', '-')[0]
            #   If the warnings are not any of the serious ones, then we may proceed
            #   This is an exact match, since the messages are generated by the program
            if current_snp.Notes != 'SNP - Sequence Mismatch' and current_snp.Notes != 'SNP aligns to a gap in CDS':
                #   Check if the SNP is coding
                if not coding:
                    current_snp.Position = 'non-coding'
                    current_snp.Silent = 'yes'
                else:
                    #   Add one because of the way python indexes lists
                    #   Cast it to str, since we want to print it out
                    current_snp.Position = str(snppos + 1)
                    #   And get some information about the amino acids and codons
                    current_snp.AA1, current_snp.AA2 = aminoacids
                    current_snp.Codon1, current_snp.Codon2 = codons
                    #   If they are the same, the change is synonymous
                    if aminoacids[0] == aminoacids[1]:
                        current_snp.Silent = 'yes'
                    else:
                        current_snp.Silent = 'no'
                        #   Get the Grantham Score if it is nonsynonymous
                        #   Recall that 'GSCORES' was defined in the 'LOOKUP' section
                        for pair in GSCORES:
                            if current_snp.AA1 in pair and current_snp.AA2 in pair:
                                current_snp.Grantham = str(GSCORES[pair])
                        #   Methionine mutations might be interesting
                        if 'M' in aminoacids:
                            current_snp.Notes = 'Methionine mutation'
                        #   If the first state is a stop codon, that could be interesting
                        if current_snp.AA1 == '*':
                            current_snp.Notes = 'Disrupted STOP codon'
                        #   If the second AA state is a stop codon, that is also interesting
                        if current_snp.AA2 == '*':
                            current_snp.Notes = 'Premature STOP codon'
    #   Back ALL the way to the base of the FOR loop
    #   The table of information to print, with the stuff cast to string where necessasry
    #   Check if the '-v' flag was set
    if ParsedArgs.verbose:
        current_annotation = [current_snp.SNPName,
                            current_snp.Organism,
                            current_snp.GenBankID,
                            current_snp.ProteinID,
                            current_snp.GeneShortName,
                            str(current_snp.Position),
                            str(current_snp.ThreeUTR),
                            str(current_snp.FiveUTR),
                            current_snp.Silent,
                            current_snp.AA1,
                            current_snp.AA2,
                            str(current_snp.Grantham),
                            str(current_snp.CDSPos),
                            current_snp.Codon1,
                            current_snp.Codon2,
                            current_snp.Ambiguity,
                            current_snp.Product,
                            current_snp.Notes,
                            current_snp.AltGene,
                            current_snp.AltOrg,
                            current_snp.ContextSeq,
                            str(current_snp.AlignScore),
                            current_snp.DateTime]
        #   Write it all to the output file
        ParsedArgs.output.write('\t'.join(current_annotation))
        #   Put a new line in there
        ParsedArgs.output.write('\n')
    else:
        #   We have to reassign some values here, for the dbSNP report
        #   namely, the synonymous and nonsynonymous data
        #   We will start with a list of three empty strings
        #   and build the comment field out of our information
        comment = ['', '', '']
        if str(current_snp.Position) == 'non-coding':
            comment[0] = 'non-coding'
        else:
            comment[1] = current_snp.AA1
            comment[2] = current_snp.AA2
            if current_snp.Silent == 'yes':
                comment[0] = 'synonymous'
            elif current_snp.Silent == 'no':
                comment[0] = 'nonsynonymous'
        current_annotation = [ 'SNP: ' + current_snp.SNPName,
                            'GENENAME: ' + current_snp.GeneShortName,
                            'ACCESSION:' + current_snp.GenBankID,
                            'COMMENT:' + ' '.join(comment),
                            'SAMPLESIZE: ',
                            'LENGTH: ?',
                            '5\'_FLANK: ' + current_snp.FiveFlank,
                            'OBSERVED: ' + current_snp.Observed,
                            '3\'_FLANK: ' + current_snp.ThreeFlank,
                            '||']
        #   Write out the block, with each entry on a newline
        ParsedArgs.output.write('\n'.join(current_annotation))
        #   And one more newline
        ParsedArgs.output.write('\n')
    #   Wait a second to keep from overloading NCBI's servers
    time.sleep(1)
