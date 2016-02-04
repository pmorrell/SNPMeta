#!/usr/bin/env python
"""
Annotates SNPs using sequences from a specified organism in NCBI's nonredundant
nucleotide sequence database. Produces a dbSNP submission report, or a tab-
delimited report. Requires Biopython.
"""

#   Import standard library modules
import sys
import os

#   Try to import the Biopython library
try:
    import Bio
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
        #   Start a SNPAnnotation class
        sys.stderr.write(
            'Annotating ' + s.id + '\n')
        anno = SNPAnnotation(s)
        #   Calculate the contextual sequence length
        anno.calculate_clen(
            args.clen,
            args.gbs,
            args.illumina)
        sys.stderr.write(
            '    Contextual sequence length: ' + str(anno.context) + '\n')

        #   Start a BlastSearch class
        blast = BlastSearch(
            args.program,
            args.evalue,
            args.max_hits,
            args.database,
            args.entrez_query)
        sys.stderr.write(
            '    Running BLAST with the following parameters: ... \n'
            '        Program:  ' + args.program + '\n'
            '        E-value:  ' + str(args.evalue) + '\n'
            '        Max hits: ' + str(args.max_hits) + '\n')
        #   Decide the command to run BLAST, web or local?
        blast.build_commandline(s)
        #   Run the search, return a handle to the results
        blastresults = blast.run_blast(s)
        #   Move on to next element if we could not get BLAST results
        if not blastresults:
            continue
        sys.stderr.write('    Done!\n')

        #   Start a GenBankHandler class to extract the GenBank record info
        genbank = GenBankHandler(args.email, blastresults)
        sys.stderr.write(
            '    Fetching GenBank records ... ')
        #   Fetch the relevant records
        genbank.fetch_gb_records()
        sys.stderr.write('Done!\n')
        #   And extract the regions and annotations
        sys.stderr.write(
            '    Extracting annotation information from GenBank records\n')
        gb_annotations = genbank.extract_annotations(args.target_organism)
        #   Set the organism name, gene name, and related locus info from the
        #   annotations.
        if gb_annotations['regions']:
            anno.organism = gb_annotations['regions'][1]
        if gb_annotations['genes']:
            anno.gene_short_name = gb_annotations['genes']
        if gb_annotations['related']:
            anno.alt_gene = gb_annotations['related'][0]
            anno.alt_org = gb_annotations['related'][1]
        #   And align the GenBank record to the SNP query sequence, using the
        #   hit frames from the extracted genbank information
        sys.stderr.write(
            '    Aligning GenBank record and query sequence ... ')
        genbank.align_genbank(s, gb_annotations['regions'][3])
        sys.stderr.write('Done!\n')

        #   Then use the alignment to calculate the aligned position of the
        #   query SNP. This method is back in the SNPAnnotation class.
        anno.get_aligned_pos(genbank.needle_out.name)
        #   Calculate whether or not the SNP is coding
        func_class = anno.functional_class(gb_annotations['regions'][2])
        #   If it is coding, get the amino acid states
        if func_class == 'Coding':
            anno.translate_codons(genbank.genbank_seq.name)
        elif func_class == 'Noncoding':
            pass
        elif func_class == 'Gap':
            pass

        print(vars(anno))

        #   Cleanup our temporary files
        os.remove(blast.blastin.name)
        os.remove(blast.blastout.name)
        os.remove(genbank.genbank_seq.name)
        os.remove(genbank.needle_out.name)
    return

main()


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
