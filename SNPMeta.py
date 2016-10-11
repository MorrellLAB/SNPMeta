#!/usr/bin/env python
"""
Annotates SNPs using sequences from a specified organism in NCBI's nonredundant
nucleotide sequence database. Produces a dbSNP submission report, or a tab-
delimited report. Requires Biopython.
"""

#   Import standard library modules
import sys
import os
import time

#   Make sure we are running Python 3
if sys.version_info.major is not 3:
    sys.exit('This script requires Python 3')

#   Try to import the Biopython library
try:
    import Bio
except ImportError:
    print('This script requires the BioPython library to be installed.')
    exit(1)

#   Import package components
from snpmeta.LookupTables.tabular_header import table_header
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

    #   Start writing the output file
    if args.outfmt == 'tabular':
        args.output.write('\t'.join(table_header) + '\n')

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
            sys.stderr.write(
                '    No BLAST hits! Try less stringent criteria.\n')
            continue
        sys.stderr.write('    Done!\n')

        #   Start a GenBankHandler class to extract the GenBank record info
        genbank = GenBankHandler(args.email, blastresults)
        sys.stderr.write(
            '    Fetching GenBank records ... ')
        #   Fetch the relevant records
        gb_success = genbank.fetch_gb_records()
        if not gb_success:
            sys.stderr.write('\n')
            anno.print_annotation(args.output, args.outfmt)
            continue
        sys.stderr.write('    Done!\n')
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
        if gb_annotations['regions']:
            genbank.align_genbank(s, gb_annotations['regions'][3])
        else:
            sys.stderr.write(
                '    No BLAST hits! Try less stringent criteria.\n')
            continue
        sys.stderr.write('Done!\n')

        #   Then use the alignment to calculate the aligned position of the
        #   query SNP. This method is back in the SNPAnnotation class.
        anno.get_aligned_pos(genbank.needle_out.name)
        #   Calculate whether or not the SNP is coding
        func_class = anno.functional_class(gb_annotations['regions'][2])
        #   If it is coding, get the amino acid states
        if func_class == 'Coding':
            anno.translate_codons(genbank.genbank_seq.name)
        elif func_class == 'Intronic':
            anno.notes = 'Intronic'
        elif func_class == 'Intergenic':
            #   Intergenic SNPs are done annotating
            pass
        elif func_class == 'Gap':
            #   As are SNPs that don't align nicely to the GenBank record
            pass
        #   Print the annotation information in the desired format
        anno.print_annotation(args.output, args.outfmt)

        #   Cleanup our temporary files
        os.remove(blast.blastin.name)
        os.remove(blast.blastout.name)
        os.remove(genbank.genbank_seq.name)
        os.remove(genbank.needle_out.name)
        #   Sleep for 1 second to avoid making requests too quickly
        time.sleep(1)
    #   Close the output handle
    args.output.close()
    return

main()
