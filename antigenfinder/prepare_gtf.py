import os
import re

import gffutils
import pyfaidx
from Bio.Seq import Seq


class TranscriptCache:
    end_pattern = r"\.[0-9]{1,2}$"

    def __init__(self, gff_file: str, fasta_file: str, cache_dir: str = None, skip_if_exists: bool = True):
        if cache_dir is not None:
            os.makedirs(cache_dir, exist_ok=True)
            gff_out_dir = os.path.join(cache_dir, 'gffdb')
        else:
            gff_out_dir = ':memory:'

        self.gffdb = parse_gff(gff_file, gff_out_dir, skip_if_exists = skip_if_exists)
        self.fasta = parse_genome_fasta(fasta_file)

        self.name_to_seq = {}


    def find_transcript_sequence(self, tid):
        transcript_version = None
        if not tid in self.name_to_seq.keys() and re.search(self.end_pattern, tid):
            transcript_version = tid.split('.')[-1]
            tid = re.sub(self.end_pattern, '', tid)

        if not tid in self.name_to_seq.keys():
            cds_features = list(self.gffdb.children(tid, featuretype='CDS', order_by='start'))
            if transcript_version:
                cds_features = list(filter(lambda x: x.attributes.get('transcript_version', [None])[0] == transcript_version, cds_features))
                if not cds_features:
                    raise Exception(f'No CDS features remained after filtering on transcript_id: {tid} / {transcript_version}')

            if not cds_features:
                raise Exception(f'Missing CDS for: {tid} / {transcript_version}')

            # Order by exon_number,
            cds_features = sorted(cds_features, key=lambda x: int(x.attributes.get('exon_number', [None])[0]))

            cds_seq = []
            for cds in cds_features:
                cds_seq.append(cds.sequence(self.fasta, use_strand=True))

            nt_seq = Seq("".join(cds_seq))

            #NOTE: the score column denotes a transcript that starts out of frame. this matters for the first codon:
            offset = cds_features[0].frame
            if offset and offset != '.':
                offset = int(offset)
                nt_seq = nt_seq[offset:]

            # NOTE: biopython complains if the input length isnt divisible by three
            if len(nt_seq) % 3 != 0:
                to_add = 3 - (len(nt_seq) % 3)
                nt_seq = nt_seq + ('N' * to_add)

            aa_seq = nt_seq.translate()

            self.name_to_seq[tid] = str(aa_seq)

        if not tid in self.name_to_seq.keys():
            raise Exception(f'Missing sequence for: {tid}')

        return self.name_to_seq[tid]


def parse_gff(gff_file, database_filename = ':memory:', skip_if_exists=True) :
    if skip_if_exists and database_filename != ':memory:' and os.path.exists(database_filename):
        print('Cached gffutils DB exists, reusing: {}'.format(database_filename))
    else:
        print('Creating gffutils DB: {}'.format(database_filename))
        gffutils.create_db(gff_file,
                           database_filename,
                           merge_strategy='create_unique',
                           disable_infer_transcripts=True,
                           disable_infer_genes=True
       )

    return gffutils.FeatureDB(database_filename)

def parse_genome_fasta(genome_fasta_file) -> pyfaidx.Fasta:
    print('Loading reference genome FASTA')
    return pyfaidx.Fasta(genome_fasta_file)
