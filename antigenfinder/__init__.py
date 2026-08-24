import argparse
from importlib.metadata import version as get_version

from antigenfinder.parse_vcf import iterate_variants
from antigenfinder.prepare_gtf import TranscriptCache

if __name__ == "__main__":
    try:
        pkg_version = get_version('antigenfinder')
    except Exception:
        pkg_version = 'unknown'

    parser = argparse.ArgumentParser(description='antigenfinder')
    parser.add_argument('-v', '--version', action='version', version=f'antigenfinder {pkg_version}')

    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    prepare_gtf_parser = subparsers.add_parser('prepare_gtf')
    prepare_gtf_parser.add_argument('--gtf_file', help='The path to the GTF file.', type=str, default=[])
    prepare_gtf_parser.add_argument('--fasta_file', help='The path to the genome FASTA file.', type=str, default=[])
    prepare_gtf_parser.add_argument('--cache_dir', help='The path where the output will be created.', type=str, default=[])

    process_vcf_parser = subparsers.add_parser('process_vcf')
    process_vcf_parser.add_argument('--vcf_file', help='The path to the VCF file.', type=str, default=[])
    process_vcf_parser.add_argument('--sample_name', help='The name of the sample to test for neoantigens.', type=str, default=[])

    process_vcf_parser.add_argument('--gtf_file', help='The path to the GTF file.', type=str, default=[])
    process_vcf_parser.add_argument('--fasta_file', help='The path to the genome FASTA file.', type=str, default=[])
    process_vcf_parser.add_argument('--cache_dir', help='The path where the parsed GTF/FASTA output read or created.', type=str, default=[])
    process_vcf_parser.add_argument('--flank_window', help='The number of amino acids to report on either side of the SNV', type=int, default=10)
    process_vcf_parser.add_argument('--output_file', help='The path where the output TSV will be written.', type=str, default=[])

    args = parser.parse_args()

    if args.subcommand == 'prepare_gtf':
        TranscriptCache(gff_file = args.gtf_file, fasta_file = args.fasta_file, cache_dir = args.cache_dir, skip_if_exists = False)
    elif args.subcommand == 'process_vcf':
        tc = TranscriptCache(gff_file = args.gtf_file, fasta_file = args.fasta_file, cache_dir = args.cache_dir, skip_if_exists = False)
        results = iterate_variants(vcf_file = args.vcf_file, sample_name = args.sample_name, transcript_cache = tc, flank_window = args.flank_window)