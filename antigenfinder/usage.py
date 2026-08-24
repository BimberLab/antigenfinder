import sys

USAGE_STRING = """
Documentation is at https://github.com/BimberLab/antigenfinder/blob/main/README.md

usage:
  python -m antigenfinder prepare_gtf <gtf_file>

  python -m antigenfinder parse_vcf <vcf_file> <

  python -m antigenfinder help"""


def print_usage_and_exit():
    print(USAGE_STRING)
    sys.exit()
