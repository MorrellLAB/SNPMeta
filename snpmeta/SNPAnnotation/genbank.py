#!/usr/bin/env python
"""A class and functions to handle fetching, parsing, and aligning sequences
fetched from GenBank."""

import sys
import time
import tempfile
from urllib.error import HTTPError
from urllib.error import URLError
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Blast import NCBIXML
from Bio import Entrez


def get_chunks(s, k):
    """Break up sequence s into chunks of k size, return a nested list. This is
    just a helper function so we don't send too many GenBank IDs in one request
    to NCBI."""
    return [s[i:i+k] for i in range(0, len(s), k)]


class GenBankHandler(object):
    """Handles the fetching and parsing of GenBank records, as well as aligning
    them to the query SNP sequence."""

    #   These are defined as class variables rather than instance variables
    #   since they will not change from run to run. We will always fetch
    #   GenBank flat files in text mode from the 'nucleotide' database.
    rettype = 'gb'
    retmode = 'text'
    database = 'nucleotide'
    #   Alignment parameters, will not change from instance to instance
    gapopen = 12
    gapextend = 3

    #   The name of the program that is being used, for NCBI's purposes.
    Entrez.tool = 'SNPMeta'

    def __init__(self, email, blast):
        """Initialize the class with the BLAST report."""
        self.gb_ids = []
        self.hit_coords = []
        self.hit_directions = []
        self.genbank_records = []
        self.get_gb_info(blast)
        Entrez.email = email
        self.genbank_seq = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='SNPMeta_GenBankSeq_',
            suffix='.fasta',
            delete=False)
        self.needle_out = tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix='SNPMeta_NeedleAlign_',
            suffix='.txt',
            delete=False)
        return

    def get_gb_info(self, resultshandle):
        """Extracts the GenBank record IDs, the hit positions, and the sequence
        orientations from the BLAST report."""
        #   Start a parser that steps through each record
        blast_records = NCBIXML.parse(resultshandle)
        #   List to hold information about our hits
        #   Step through the BLAST records
        for record in blast_records:
            #   Step through each alignment in each record
            for alignment in record.alignments:
                #   Then the HSPs in each alignment
                for hsp in alignment.hsps:
                    #   The start and end positions of each hit
                    hit_coords = (hsp.sbjct_start, hsp.sbjct_end)
                    #   Split on the '|' character, genbank ID is last in the
                    #   list have to use -2 instead, because of the trailing '|'
                    #   in the XML report
                    hit_gbid = alignment.title.split('|')[-2]
                    #   Relative directions of the sequences
                    hit_directions = hsp.frame
                    break
                #   Tack the IDs, coordinates, and directions onto our lists
                self.gb_ids.append(hit_gbid)
                self.hit_coords.append(hit_coords)
                self.hit_directions.append(hit_directions)
        #   Finished with this file
        resultshandle.close()
        return

    def fetch_gb_records(self):
        """Fetches the GenBank records from NCBI."""
        #   If we get here without any BLAST hits, then we return a False
        if len(self.gb_ids) == 0:
            return False
        #   Build lists to hold the request ranges for each record
        starts = []
        ends = []
        for c in self.hit_coords:
            #   We will request ~5kb total, 2.5kb upstream and 2.5k downstream
            #   5kb upstream
            hitstart = min(c) - 2500
            #   5kb downstream
            hitend = max(c) + 2500
            #   Weird things happen if a start request is negative,
            #   so set it to 1 if it is negative
            if hitstart < 0:
                hitstart = 1
        #   Store as strings, since they have to be sent to NCBI as strings
        starts.append(str(hitstart))
        ends.append(str(hitend))
        #   This is an unideal workaround, but it should work for now
        #   If there are more than 20 records to parse, then we could
        #   Get the chunks that we will query NCBI with
        #   We choose 20 at a time, that should be reasonable
        genbank_chunks = get_chunks(self.gb_ids, 20)
        start_chunks = get_chunks(starts, 20)
        end_chunks = get_chunks(ends, 20)
        #   Actually query NCBI!
        #   Entrez.efetch() forms the proper query and sends it to NCBI. Returns
        #   a handle to concatenated GenBank records. Just like with BLAST
        #   searching, we wrap in try/except/finally with a counter.
        #   Keep an empty string to append the records to and parse it later
        batch = 1
        for gid, st, en in zip(genbank_chunks, start_chunks, end_chunks):
            sys.stderr.write(
                '\n        Requesting batch ' + str(batch) + ' of ' +
                str(len(genbank_chunks)) + ' (' + str(len(self.gb_ids)) +
                ' records total)' + '\n')
            tries = 0
            success = False
            while not (success or tries >= 3):
                tries += 1
                try:
                    #   Fetch the records
                    handle = Entrez.efetch(
                        db=self.database,
                        id=','.join(gid),
                        rettype=self.rettype,
                        retmode=self.retmode,
                        seq_start=','.join(st),
                        seq_stop=','.join(en))
                    #   Parse them
                    temp_genbank = SeqIO.parse(handle, 'genbank')
                    #   And tack them onto our list
                    self.genbank_records += list(temp_genbank)
                    #   wait for a second, to avoid abusing the NCBI servers
                    time.sleep(1)
                    handle.close()
                    success = True
                except HTTPError as e:
                    sys.stderr.write(
                        'Caught ' +
                        str(e.code) +
                        ' with reason ' +
                        str(e.reason) +
                        '. Retrying in 5 seconds ...\n')
                    time.sleep(5)
                except (URLError, ConnectionRefusedError) as e:
                    #   Sometimes our connection is refused. We handle this by
                    #   retrying after 5 seconds. The reason for this error is
                    #   usually that the service that listens on the requested
                    #   port can't handle the request.
                    sys.stderr.write(
                        'Caught URL error, connection refused. Retrying in 5 '
                        'seconds ...\n')
                    time.sleep(5)
                finally:
                    #   Print out a message if we could not fetch the list of
                    #   GenBank records
                    if not success and tries >= 3:
                        sys.stderr.write(
                            'Failed to fetch GenBank records after three '
                            'attempts. Moving on ...\n')
        return success

    def extract_annotations(self, target):
        """Extracts the sequences and annotated regions from the GenBank
        records. Also sets the secondary hit information. Filters based on
        whether records come from target organisms or not."""
        #   Empty dict to hold the extracted annotation information
        annotations = {
            'records': None,
            'regions': None,
            'genes': None,
            'related': None
            }
        #   Step through our list for the first time
        #   On this pass, we will only get hits from our target organisms, and
        #   out of those, we will only take the hits with an annotated CDS.
        for d, r in zip(self.hit_directions, self.genbank_records):
            #   Save the IDs, since we will want to print them out later
            annotations['records'] = r.id
            #   We want to extract any genes that may be on the record
            annotations['genes'] = [f.qualifiers.get('gene', [None])[0] for f in r.features if f.type == 'gene']
            #   Extract the organism from the info
            organism = r.annotations['organism']
            #   If the user has supplied some target organisms...
            if target:
                #   Skip uninteresting organisms
                if organism in target:
                    #   A list of all features in the record
                    types = [f.type for f in r.features]
                    #   We only care about those that have a CDS annotation
                    if 'CDS' in types:
                        #   Save the name, organism of origin, the features,
                        #   and the direction it hit the genbank record
                        annotations['regions'] = (
                            r.name,
                            organism,
                            r.features,
                            d)
                        #   And save the genbank sequence for later
                        SeqIO.write(r, self.genbank_seq.name, 'fasta')
                        #   We popped the first good one off, so we break
                        break
            else:
                #   If no target organism specifed, we just take the best hit
                types = [f.type for f in r.features]
                if 'CDS' in types:
                    annotations['regions'] = (r.name, organism, r.features, d)
                    SeqIO.write(r, self.genbank_seq.name, 'fasta')
                    break

        #   The second pass through the list of hits, where we will find genes
        #   in other similar sequences, in case this information is useful for
        #   interpretation.
        for d, r in zip(self.hit_directions, self.genbank_records):
            #   Extract the organism in the same way as above
            organism = r.annotations['organism']
            #   Now, we are interested in all the features
            for f in r.features:
                #   If there is a gene, then we will save it
                if f.type == 'gene':
                    #   Pull the name of the gene from the feature. The get
                    #   method will use the second argument as a default if the
                    #   first does not exist as a key in the features dict.
                    genename = f.qualifiers.get('gene', '-')[0]
                    #   Build a tuple of the related organism information
                    annotations['related'] = (genename, organism)
                    #   And again, break once we find one
                    break
        #   Return our information. This will be fed into the SNPAnnotation
        #   class to get estimated functional class.
        return annotations

    def align_genbank(self, query, frames):
        """Aligns the chosen GenBank record to the query sequence with the
        'needle' program in EMBOSS."""
        #   First, we check the frames of the hits from the BLAST record. If
        #   they are different, then we have to RC one of the sequences. We
        #   choose to RC the query, since it is much simpler, and will not
        #   introduce errors. The frames are given as integers with the values
        #   {-3, -2, -1, 1, 2, 3}. If they are the same sign, don't RC, if they
        #   are not, then RC.
        if (frames[0] > 0 and frames[1] > 0):
            snp_seq = query.seq
        else:
            snp_seq = query.seq.reverse_complement()
        #   Build the 'needle' command
        needle_cmd = NeedleCommandline(
            asequence=self.genbank_seq.name,
            bsequence='asis:' + str(snp_seq),
            gapopen=self.gapopen,
            gapextend=self.gapextend,
            outfile=self.needle_out.name)
        #   Run it!
        needle_cmd()
        return
