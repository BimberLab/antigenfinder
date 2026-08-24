# antigenfinder

[![Run Tests](https://github.com/BimberLab/antigenfinder/actions/workflows/ci.yml/badge.svg)](https://github.com/BimberLab/antigenfinder/actions/workflows/ci.yml)

## Overview

antigenfinder is a specialized tool designed to identify potentially immunogenetic variants between samples (e.g., the donor/recipient pair in transplantation). In addition to identifying variants unique to the source sample, it inspects each variant, identifies linked flanking variants, and reports the AA sequence within the local region. The result is a table of putative T cell epitopes, which can be further processed using HLA prediction tools or other analyses.

## Usage

antigenfinder expects a VCF with at least two samples. That VCF must be processed with both [SnpEff](https://pcingola.github.io/SnpEff/) to calculate AA consequences, and [whatshap](https://whatshap.readthedocs.io/) to provide read-backed phasing.

Basic usage requires two steps:

First, you can optionally prepare the FASTA/GTF. This runs the initial steps to parse the CDS from the GTF and can store the resulting file on disk. Doing this will save time if you run the tool multiple times. You can optionally provide '--debug_output', in which case the tool will write a FASTA with the inferred CDS AA sequences for each transcript.
```
python -m antigenfinder prepare_gtf \
    --gtf_file myGtf.gtf \
    --fasta_file myFasta.fasta \
    --cache_dir /path/to/cached_files \
    --skip_processing_if_cached \
    --debug_output /optional/output_with_cds.fasta
```

Now run the tool against the VCF. Note: if you already cached the GTF/FASTA info to --cache_dir using 'prepare_gtf', then it will be re-used. Note, by default it will report AAs with a 10 AA window on either side; however, this is configurable with '--aa_flank_window'. The '--source_sample' argument refers to the sample in the VCF to test for potential neoantigens. See below for more information on variant processing.   
```
python -m antigenfinder process_vcf \
    --vcf_file /path/to/myVcf.vcf.gz \
    --source_sample TheSample \
    --output_file /path/to/myOutput.txt \
    --gtf_file myGtf.gtf \
    --fasta_file myFasta.fasta \
    --cache_dir /path/to/cached_files \
    --skip_processing_if_cached \
    --debug_output /optional/output_with_cds.fasta

```

## Processing of Variants:
The tool follows the following logic:
- For each position, it will compare the variants in the --source_sample all other samples in the VCF. If you want to perform a pairwise contrast, only include two samples in the VCF.
- Filtered variants are ignored. Likewise, if the source_sample's genotype is no-call, that position will be skipped. Also, at least one sample besides the source sample must have a called genotype for that site to be included. 
- Any NT present in the source_sample and absent in all other samples is considered. This could include the wild-type allele.
- Per variant, the tool inspects the SnpEff annotation. It iterates all annotated protein consequences. In practice, this can mean a variant is reported more than once, when multiple isoforms exist for a given protein.
- If the variant alters protein coding, the tool builds the local AA sequence around this variant. The size of that region is defined by --aa_flank_window (if the window is 10AA, it will inspect +/-10 AA from the variant).
- When creating that AA sequence, it uses the whatshap annotations to find linked variants. If there is a linked variant, it will include this change in the resulting AA.
- It produces a tab-delimited file with one line for each protein change, and lots of information about that site. In this file, it reports the AA_Sequence, where any AAs altered by the variant are shown in lowercase  
- Because there are many reasons for redundant lines in the output (including overlapping isoforms), in many cases it will make sense to subset to the core columns (e.g., just GeneName and AA_Sequence) and remove duplicates.
- Note: there is a final column called 'Messages', which contains information about sites that it could not parse. This primarily occurs for complex INDELs. The tool always reports the SnpEff annotations as a column, in case there is question about the predicted AA.