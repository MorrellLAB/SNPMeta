#!/usr/bin/env python
"""Class to store SNP annotation information, calculate positions from a
pairwise alignment, get coding state using a CDS annotation, and translate
codons. Prints annotation output in two formats."""

import math
import re
import datetime
from snpmeta.LookupTables.iupac import IUPAC
from snpmeta.LookupTables.grantham import GSCORES
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Alphabet.IUPAC import ambiguous_dna


class SNPAnnotation(object):
    """A class to store and emit the SNP annotation data."""

    def __init__(self, sequence, arg_dict):
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
        self.context, self.amb = self.calculate_clen(
            arg_dict.clen,
            arg_dict.gbs, 
            arg_dict.illumina)
        self.aligned_seq = None
        self.aligned_pos = None
        self.cds_feat = None
        self.cds_snp_pos = None
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
                    return (cl, self.query_seq[cl])
                elif self.query_seq[-(cl + 1)] in IUPAC:
                    #   Next, check from the right. If this does give us an
                    #   ambiguity, convert the contextual length to a left-
                    #   based offset.
                    #   Drop a warning into the 'context_seq' field
                    self.context_seq = 'Short_on_3prime'
                    return (
                        len(self.query_seq) - (cl + 1),
                        self.query_seq[-(cl + 1)]
                        )
                else:
                    #   A last resort is to iterate through until we find the
                    #   first ambiguity. This is almost guaranteed to be wrong
                    self.context_seq = 'Shor_on_both'
                    for i, b in enumerate(self.query_seq):
                        if b in IUPAC:
                            return (i, b)
        elif illumina:
            #   The next flag we check is whether or not the SNPs are given with
            #   Illumina formatting. If so, we search for the [ character that
            #   denotes the query states, and find the offset. We also remove
            #   the non-nucleotide characters from the sequence, as those will
            #   cause errors later.
            offset = re.search(r'\[', str(self.query_seq)).start()
            #   Build the fixed sequence
            seq_parts = re.split(r'[\[|\/|\]]', str(self.query_seq))
            for a, nona in IUPAC.items():
                #   Search for the ambiguity that contains the two query SNP
                #   states, and drop it in
                if (seq_parts[1] in nona) and (seq_parts[2] in nona):
                    snp = a
                    break
            self.query_seq = Seq(seq_parts[0] + snp + seq_parts[-1])
            return (offset, snp)
        elif gbs:
            #   The final check is whether the SNPs are from GBS tags or not. If
            #   they are, we iterate until we find the ambiguity, and return the
            #   offset of the first one we find.
            self.context_seq = 'GBS'
            for i, b in enumerate(self.query_seq):
                if b in IUPAC:
                    return (i, b)

    def get_aligned_pos(self, alignment):
        """Get the position of the SNP in the alignment."""
        #   First, extract the score of the alignment. This is stored on a line
        #   with the form
        #       #Score: xxx.x
        with open(alignment, 'r') as f:
            for line in f:
                if line.startswith('# Score:'):
                    self.align_score = line.split(':')[-1].strip()
                    break
        aln = AlignIO.read(alignment, 'emboss')
        #   First sequence is the GenBank sequence, we store it so that we can
        #   extract proper features from it later
        self.aligned_seq = aln[0].seq
        #   And the ID should be the same as from the GenBank record
        self.genbank_id = aln[0].id
        #   Don't bother storing the aligned SNP sequence long-term. We are just
        #   going to REs on it to calculate positions.
        snp_seq = str(aln[1].seq)
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

    def functional_class(self, features):
        """Calculates where the SNP is coding or noncoding."""
        #   Iterate through the features extracted from the GenBank record
        for f in features:
            if f.type == 'gene':
                #   Get the short name of the gene.
                #   For some reason, it is stored as a list
                #   with only one element
                self.gene_short_name = f.qualifiers.get('gene', '-')[0]
            if f.type == 'CDS':
                #   Save the postions of the alignment that are not gaps in the
                #   GenBank sequence
                nogap_gb = [i
                            for
                            i, b
                            in
                            enumerate(self.aligned_seq)
                            if b != '-']
                #   Then, save the corresponding feature positions
                cds_pos = [pos
                           for
                           i, pos
                           in
                           enumerate(nogap_gb)
                           if i in f]
                #   Save the organism, product, and protien ID information from
                #   the CDS annnotation.
                self.product = f.qualifiers.get('product', '-')[0]
                #   Check to see if the CDS has a GenBank Protein ID
                self.protein_id = f.qualifiers.get('protein_id', '-')[0]
                #   Check that the aligned position does not occur in a gap
                if self.aligned_pos not in nogap_gb:
                    #   Return text that lets us know in the main program what
                    #   sort of errors there are, if any.
                    self.cds_feat = None
                    self.cds_snp_pos = None
                    return 'Gap'
                elif self.aligned_pos in cds_pos:
                    #   If it is in the feature, we have a coding SNP. Set the
                    #   feature and return text that says we coding
                    self.cds_feat = f
                    #   Calculate where in the aligned feature the SNP is. We
                    #   have to check the feature start, since it can be fuzzy.
                    #   With fuzzy starts, we cannot be sure of the frame of the
                    #   CDS. We do know, however, that the CDS must end on a
                    #   codon triplet.
                    if '>' in str(f.location.start):
                        self.cds_snp_pos = (f.location.end -
                                            cds_pos.index(self.aligned_pos))
                    else:
                        self.cds_snp_pos = cds_pos.index(self.aligned_pos)
                    return 'Coding'
                elif self.aligned_pos not in cds_pos:
                    #   If not in the feature, it is noncoding, but we can
                    #   estimate how far upstream or downstream the nearest
                    #   piece of coding sequence is.
                    self.cds_feat = None
                    self.cds_snp_pos = None
                    self.silent = 'Yes'
                    #   The distance from the start of the CDS to the aligned
                    #   position
                    up_dist = int(min(cds_pos) - self.aligned_pos)
                    #   The distance of the aligned position, past the end of
                    #   the CDS
                    down_dist = int(self.aligned_pos - max(cds_pos))
                    #   If the SNP is intronic, both of these will be negative
                    if up_dist < 0 and down_dist < 0:
                        return 'Intronic'
                    elif up_dist < 0 and down_dist > 0:
                        #   If the upstream distance is negative and the
                        #   downstream distance is positive, then the SNP is
                        #   3' of the gene
                        self.three_flank = down_dist
                        return 'Intergenic'
                    elif up_dist > 0 and down_dist < 0:
                        #   If upstream distance is positive and downstream
                        #   distance is negative, it is 5' of the gene.
                        self.five_flank = up_dist
                        return 'Intergenic'

    def translate_codons(self, gbseq):
        """Get the codon that contains the SNP, and translate the two states."""
        #   Read the GenBank sequence back off the disk. We need to do this
        #   to get an accurate translation of the coding sequence.
        ref = SeqIO.read(gbseq, 'fasta')
        #   Get the translation table that the feature uses, so we can be sure
        #   that the amino acid sequence is correct. Default to 1.
        table_id = self.cds_feat.qualifiers.get('transl_table', '1')[0]
        #   Extract the coding sequence from the GenBank sequence, and translate
        #   We use the extract() method of SeqFeature objects to get the CDS.
        cds_tr = self.cds_feat.extract(ref.seq).translate(table=table_id)
        #   Then, calculate the amino acid residue that the SNP occurs in
        self.cds_pos = int(math.floor(self.cds_snp_pos/3))
        #   And get the position within the codon
        #   If pos % 3 = 0, then we have first position
        #   If pos % 3 = 1, then we have second position
        #   If pos % 3 = 2, then we have third position
        self.position = self.cds_snp_pos % 3
        #   Get the codon based on the relative position
        #   A codon is given by
        #       [pos - (pos%3), pos + (2 - (pos%3))+1]
        #   If SNP is in 4th CDS position (index 3 in Python)
        #       [3-0, 3+(2-0)+1] -> [3, 6] -> [3, 4, 5]
        #   If SNP is in 5th CDS position (index 4 in Python)
        #       [4-1, 4+(2-1)+1] -> [3, 6] -> [3, 4, 5]
        #   If SNP is in 6th CDS position (index 5 in Python)
        #       [5-2, 5+(2-2)+1] -> [3, 6] -> [3, 4, 5]
        #   7th CDS position (6 in Python)
        #       [6-0, 6+(2-0)+1] -> [6, 9] -> [6, 7, 8]
        #   And so on...
        lower_bound = self.cds_snp_pos - self.position
        upper_bound = self.cds_snp_pos + (2 - self.position) + 1
        codon_positions = list(self.cds_feat)[lower_bound:upper_bound]
        #   This is a list because we need to 'mutate' it later.
        gb_codon = [ref[pos] for pos in codon_positions]
        #   Also determine the two nucleotide staes
        #   Check the strand of the feature. If it is -1, then we have to RC the
        #   codon, and complement the ambiguity.
        if self.cds_feat.strand == -1:
            gb_codon = list(Seq(''.join(gb_codon)).reverse_complement())
            self.amb = str(Seq(self.amb, ambiguous_dna).complement())
            #   We also have to change the codon positions
            #   Third position becomes first, and first becomes third
            if self.position == 2:
                self.position = 0
            elif self.position == 0:
                self.position = 2
            #   The amino acid position has to be recalculated
            #   But it is just measuring from a different end
            self.cds_pos = len(cds_tr) - self.cds_pos
        else:
            #   If the strand is +1, then we don't have to re-count things
            pass
        #   Build a copy of the GenBank sequence codon, so we have a 'ref' and
        #   and 'alt' state
        snp_codon = gb_codon[:]
        #   And drop the ambiguity into the right place
        snp_codon[self.position] = self.amb
        #   Translate the codon
        #   Can shortcut, since there should only be one ambiguity
        #   First, convert into a string
        gb_codon_str = ''.join(gb_codon)
        snp_codon_str = ''.join(snp_codon)
        #   Start a list to hold our unambiguous codons
        resolved_codons = []
        #   Next, resolve the ambiguities
        for stdnuc in IUPAC[self.amb]:
            #   Replace the ambiguity with a standard nucleotide
            resolved_codons.append(snp_codon_str.replace(self.amb, stdnuc))
        #   Which codon has the non-GenBank state?
        non_gb_codon = [c for c in resolved_codons if c != gb_codon_str][0]
        self.codon_1 = gb_codon_str
        self.codon_2 = non_gb_codon
        #   Then, we translate them
        self.aa_1 = str(Seq(gb_codon_str).translate(table=table_id))
        self.aa_2 = str(Seq(non_gb_codon).translate(table=table_id))
        #   Last, if the amino acids are the same, then the SNP is silent
        if self.aa_1 == self.aa_2:
            self.silent = 'Yes'
        else:
            self.silent = 'No'
            #   Get Grantham scores for nonsnyonmyous SNPs
            for m, s in GSCORES.items():
                if self.aa_1 in m and self.aa_2 in m:
                    self.grantham = s
                    break
        return

    def print_annotation(self, outhandle, outformat):
        """Print the SNP annotation information."""
        if outformat == 'tabular':
            towrite = '\t'.join(
                [
                    self.snp_name,
                    self.organism,
                    self.genbank_id,
                    self.protein_id,
                    self.gene_short_name,
                    str(self.position),
                    str(self.three_flank),
                    str(self.five_flank),
                    self.silent,
                    self.aa_1,
                    self.aa_2,
                    str(self.grantham),
                    str(self.cds_pos),
                    self.codon_1,
                    self.codon_2,
                    self.amb,
                    self.product,
                    self.notes,
                    self.alt_gene,
                    self.alt_org,
                    self.context_seq,
                    self.align_score,
                    self.date_time
                ]) + '\n'
        elif outformat == 'dbsnp':
            #   We have to reassign some values here, for the dbSNP report
            #   namely, the synonymous and nonsynonymous data
            #   We will start with a list of three empty strings
            #   and build the comment field out of our information
            comment = ['', '', '']
            if not self.cds_feat:
                comment[0] = 'non-coding'
            else:
                comment[1] = self.aa_1
                comment[2] = self.aa_2
            if self.silent == 'Yes':
                comment[0] = 'synonymous'
            elif self.silent == 'No':
                comment[0] = 'nonsynonymous'
            fiveflank = len(self.query_seq) - (self.context + 1)
            threeflank = self.context
            towrite = '\n'.join(
                [
                    'SNP: ' + self.snp_name,
                    'GENENAME: ' + self.gene_short_name,
                    'ACCESSION:' + self.genbank_id,
                    'COMMENT:' + ' '.join(comment),
                    'SAMPLESIZE:',
                    'LENGTH: ?',
                    '5\'_FLANK: ' + str(fiveflank),
                    'OBSERVED: ' + '/'.join(IUPAC[self.amb]),
                    '3\'_FLANK: ' + str(threeflank),
                    '||']) + '\n'
        outhandle.write(towrite)
        outhandle.flush()
        return
