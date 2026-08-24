import os
import unittest
import tempfile
import urllib.request
import pathlib
import pysam

from antigenfinder.prepare_gtf import TranscriptCache
from parse_vcf import process_vcf


def download_inputs():
    data_dir = pathlib.Path(__file__).parent / "data"
    local_cache = data_dir / "local_cache"
    local_cache.mkdir(parents=True, exist_ok=True)

    genome_url = 'https://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
    local_genome = local_cache / 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
    print('Downloading Genome FASTA to: {}'.format(local_genome))
    if not os.path.exists(local_genome):
        urllib.request.urlretrieve(genome_url, local_genome)
    else:
        print('using existing file')

    gtf_url = 'https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz'
    local_gtf = local_cache / 'Homo_sapiens.GRCh38.108.gtf.gz'
    print('Downloading GTF to: {}'.format(local_gtf))
    if not os.path.exists(local_gtf):
        urllib.request.urlretrieve(gtf_url, local_gtf)
    else:
        print('using existing file')

    return [str(local_genome), str(local_gtf)]


class TestDataProcessing(unittest.TestCase):

    def test_gtf_parsing(self):
        cache_dir = tempfile.mkdtemp(prefix = 'antigenfinder.')
        downloaded_genomes = download_inputs()
        tc = TranscriptCache(gff_file=downloaded_genomes[1], fasta_file=downloaded_genomes[0], cache_dir=cache_dir)

        self.assertEqual(len(tc.find_transcript_sequence('ENST00000374222.6')), 520)


    def test_vcf_parsing(self):
        downloaded_genomes = download_inputs()

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
