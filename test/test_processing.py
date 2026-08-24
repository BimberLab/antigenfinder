import os
import unittest

import pysam

from antigenfinder.prepare_gtf import TranscriptCache
from parse_vcf import process_vcf


class TestDataProcessing(unittest.TestCase):

    def test_gtf_parsing(self):
        gtf = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/Homo_sapiens.GRCh38.108.gtf'
        fasta = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/129_Human_GRCh38.p13_Ensembl.fasta'
        cache_dir = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/antigenfinder.cache'
        #cache_dir = tempfile.mkdtemp(prefix = 'antigenfinder.')

        tc = TranscriptCache(gff_file=gtf, fasta_file=fasta, cache_dir=cache_dir)

        seq = tc.find_transcript_sequence('ENST00000374222.6')

        self.assertTrue(seq is not None)


    def test_vcf_parsing(self):
        gtf = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/Homo_sapiens.GRCh38.108.gtf'
        fasta = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/129_Human_GRCh38.p13_Ensembl.fasta'
        cache_dir = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/antigenfinder.cache'

        tc = TranscriptCache(gff_file=gtf, fasta_file=fasta, cache_dir=cache_dir)

        vcf = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/SHIELD.802761.vcf.gz'
        if not os.path.exists(vcf + '.tbi'):
            print('Making VCF index')
            pysam.tabix_index(vcf, preset="vcf", force=True)

        out_file = '/mnt/c/Users/Bimber/Downloads/SHIELD_VCF/SHIELD.802761.output.txt'
        results = process_vcf(vcf_file=vcf, source_sample='SHIELD_002_Recipient', transcript_cache=tc, max_records_to_process=250000, output_file=out_file)
        self.assertEqual(len(results), 5486, 'Incorrect number of results')


if __name__ == "__main__":
    unittest.main()
