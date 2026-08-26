import gzip
import os
import pathlib
import shutil
import tempfile
import unittest
import urllib.request

import pysam

from antigenfinder.prepare_gtf import TranscriptCache
from antigenfinder.parse_vcf import process_vcf


def get_local_cache():
    data_dir = pathlib.Path(__file__).parent / "data"
    local_cache = data_dir / "local_cache"
    local_cache.mkdir(parents=True, exist_ok=True)
    return local_cache

def download_inputs():
    local_cache = get_local_cache()
    genome_url = 'https://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
    local_genome_gzip = local_cache / 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
    local_genome = local_cache / 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fasta'
    print('Downloading Genome FASTA to: {}'.format(local_genome))
    if not os.path.exists(local_genome):
        urllib.request.urlretrieve(genome_url, local_genome_gzip)
        print('unzipping genome FASTA')
        with gzip.open(local_genome_gzip, 'rb') as f_in:
            with open(local_genome, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(local_genome_gzip)
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

def ensure_vcf_index(vcf):
    if not os.path.exists(vcf + '.tbi'):
        print('Making VCF index')
        pysam.tabix_index(vcf, preset="vcf", force=True)

class TestDataProcessing(unittest.TestCase):

    def test_gtf_parsing(self):
        cache_dir = tempfile.mkdtemp(prefix = 'antigenfinder.')
        downloaded_genomes = download_inputs()
        tc = TranscriptCache(gff_file=downloaded_genomes[1], fasta_file=downloaded_genomes[0], cache_dir=cache_dir)

        self.assertEqual(len(tc.find_transcript_sequence('ENST00000374222.6')), 520)

    def test_vcf_parsing(self):
        data_dir = pathlib.Path(__file__).parent / "data"
        downloaded_genomes = download_inputs()

        cache_dir = str(get_local_cache() / 'antigenfinder.cache')
        tc = TranscriptCache(gff_file=downloaded_genomes[1], fasta_file=downloaded_genomes[0], cache_dir=cache_dir)

        vcf = str(data_dir / 'test.vcf.gz')
        ensure_vcf_index(vcf)

        out_file = 'antigenfinder.output.txt'
        stats_collector = process_vcf(vcf_file=vcf, source_sample='Sample1', transcript_cache=tc, output_file=out_file)
        self.assertEqual(stats_collector.total_sites_inspected, 248873, 'Incorrect total_sites_inspected')
        self.assertEqual(stats_collector.total_discordant_sites, 716, 'Incorrect total_discordant_sites')
        self.assertEqual(stats_collector.total_protein_altering, 271, 'Incorrect total_protein_altering')
        self.assertEqual(len(stats_collector.unique_transcripts), 656, 'Incorrect unique_transcripts')
        self.assertEqual(len(stats_collector.unique_gene_names), 191, 'Incorrect unique_gene_names')

if __name__ == "__main__":
    unittest.main()
