#!/usr/bin/env python
"""
Annotates SNPs using sequences from a specified organism in NCBI's nonredundant
nucleotide sequence database. Produces a dbSNP submission report, or a tab-
delimited report. Requires Biopython.
"""

#   Import standard library modules
import sys
import os

#   Import the Biopython libraries
try:
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
from snpmeta.Environment import setup_env
from snpmeta.ArgumentHandling import parse_args
from snpmeta.ArgumentHandling import validate_args
from snpmeta.SNPAnnotation.snp_annotation import SNPAnnotation
from snpmeta.SNPAnnotation.blast_search import BlastSearch
from snpmeta.SNPAnnotation.genbank import GenBankHandler


def main():
    """Main function."""
    #   First, parse arguments
    parser, args = parse_args.parse_args()
    #   And then validate
    error_code = validate_args.validate_args(args)
    if error_code:
        if error_code == 1:
            parser.print_help()
            sys.stderr.write('\nError! Some required inputs were not given.\n')
            exit(1)
        elif error_code == 2:
            sys.stderr.write('Email entered is invalid.\n')
            exit(1)
        elif error_code == 3:
            sys.stderr.write('Contextual sequence length entered is invalid.\n')
            exit(1)
        elif error_code == 4:
            sys.stderr.write(
                'Directory entered is not readable or does not exist.\n')
            exit(1)
        elif error_code == 5:
            sys.stderr.write(
                'FASTA file entered is an invalid FASTA file.\n')
            exit(1)

    #   Done checking arguments. Start annotating SNPs.
    for s in setup_env.build_targets(args.fasta_file, args.dir, args.no_blast):
        anno = SNPAnnotation(s)
        anno.calculate_clen(
            args.clen,
            args.gbs,
            args.illumina)
        blast = BlastSearch(
            args.program,
            args.evalue,
            args.max_hits,
            args.database,
            args.entrez_query)
        blast.build_commandline(s)
        genbank = GenBankHandler(None)
        #   Print some things to make sure our classes are doing what they are
        #   supposed to be doing
        print(blast.commandline)
        print(anno.snp_name, anno.query_seq)
        print(genbank.genbank_seq.name)
        #   Cleanup our temporary files
        os.remove(blast.blastin.name)
        os.remove(blast.blastout.name)
        os.remove(genbank.genbank_seq.name)
        os.remove(genbank.needle_out.name)
    return

main()
###############################################################################
#       CONSTANTS AND STRINGS
###############################################################################
#   The email address and tool name to send to NCBI
Entrez.tool = 'SNPMeta, using BioPython'
Entrez.email = ParsedArgs.email


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
