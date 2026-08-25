import argparse
from importlib.metadata import version as get_version

import argcomplete

from antigenfinder.parse_vcf import process_vcf
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
    prepare_gtf_parser.add_argument('--gtf_file', help='The path to the GTF file.', type=str, required=True)
    prepare_gtf_parser.add_argument('--fasta_file', help='The path to the genome FASTA file.', type=str, required=True)
    prepare_gtf_parser.add_argument('--cache_dir', help='The path where the output will be created.', type=str, required=False)
    prepare_gtf_parser.add_argument('--skip_processing_if_cached', help='An optional flag to skip re-processing. If the GFF DB exists in --cache_dir, it will be re-used, instead of re-processing the GTF, ', action="store_true", default=False, required=False)
    prepare_gtf_parser.add_argument('--debug_output', help='If provided, all inferred AA sequences will be written to this FASTA file. This can help troubleshoot processing of a GTF.', type=str, required=False)

    process_vcf_parser = subparsers.add_parser('process_vcf')
    process_vcf_parser.add_argument('--vcf_file', help='The path to the VCF file.', type=str, required=True)
    process_vcf_parser.add_argument('--source_sample', help='The name of the sample to test for neoantigens.', type=str, required=True)
    process_vcf_parser.add_argument('--aa_flank_window', help='The number of amino acids to report on either side of the variant', type=int, default=10, required=False)
    process_vcf_parser.add_argument('--output_file', help='The path where the output TSV will be written.', type=str, required=True)

    process_vcf_parser.add_argument('--gtf_file', help='The path to the GTF file.', type=str, required=True)
    process_vcf_parser.add_argument('--fasta_file', help='The path to the genome FASTA file.', type=str, required=True)
    process_vcf_parser.add_argument('--cache_dir', help='The path where the parsed GTF/FASTA output read or created.', type=str, required=False)
    process_vcf_parser.add_argument('--skip_processing_if_cached', help='An optional flag to skip re-processing. If the GFF DB exists in --cache_dir, it will be re-used, instead of re-processing the GTF, ', action="store_true", default=False, required=False)

    args = parser.parse_args()

    if args.subcommand == 'prepare_gtf':
        tc = TranscriptCache(gff_file = args.gtf_file, fasta_file = args.fasta_file, cache_dir = args.cache_dir, skip_if_exists = args.skip_processing_if_cached)

        if args.debug_output:
            print('Writing FASTA with all inferred AA sequences to {}'.format(args.debug_output))
            with open(args.output_file, mode='w') as out:
                total = 0
                for transcript in tc.gffdb.features_of_type('transcript'):
                    out.write(transcript.id + '\n')
                    aa = tc.find_transcript_sequence(transcript.id)
                    out.write(aa + '\n')
                    total += 1

                print('Total written: {}'.format(total))

    elif args.subcommand == 'process_vcf':
        tc = TranscriptCache(gff_file = args.gtf_file, fasta_file = args.fasta_file, cache_dir = args.cache_dir, skip_if_exists = args.skip_processing_if_cached)
        results = process_vcf(vcf_file = args.vcf_file, source_sample = args.source_sample, transcript_cache = tc, aa_flank_window = args.aa_flank_window, output_file=args.output_file)