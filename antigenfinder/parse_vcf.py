import math

import pysam
from pysam import VariantRecord

from antigenfinder.prepare_gtf import TranscriptCache
from snpeff import SnpEffAnn


class Hit:
    def __init__(self, record: VariantRecord, sample_name: str, allele_idx: int):
        self.sample_name = sample_name
        self.chrom = record.chrom
        self.pos = record.pos
        self.ref = record.ref
        self.allele_idx = allele_idx
        self.nt_change = record.alleles[allele_idx]

def process_variant(record: VariantRecord, source_sample: str, record_buffer: dict[int, list[VariantRecord]], transcript_cache: TranscriptCache, window_size_nt: int) -> list[Hit]:
    annotation = record.info.get('ANN')
    if not annotation:
        raise Exception('No ANN value for variant: {}'.format(record))

    ret = []
    for allele_idx in get_unique_alleles(record, source_sample):
        if allele_idx is None:
            print('Allele IDX is empty: {}'.format(allele_idx))
            print(record.samples[source_sample]['GT'])
            print(get_unique_alleles(record, source_sample))

        ann = SnpEffAnn(annotation, record, allele_idx)
        tids = ann.get_transcript_ids()
        if not tids:
            continue

        flanking_variants = []
        ps = record.samples[source_sample].get('PS')
        if ps:
            min_nt_pos = record.pos - window_size_nt
            max_nt_pos = record.pos + window_size_nt
            for nt_pos in record_buffer.keys():
                if min_nt_pos <= nt_pos <= max_nt_pos:
                    for fv in record_buffer[nt_pos]:
                        if fv.samples[source_sample].get('PS') == ps:
                            flanking_variants.append(fv)

        for tid in tids:
            aa_seq = transcript_cache.find_transcript_sequence(tid)
            aa_positions1 = ann.get_aa_positions(tid)
            if not aa_positions1:
                # TODO: capture this better
                continue

            aa_changes = ann.get_aa_changes(tid)
            if len(aa_changes) != len(aa_positions1):
                raise Exception('AA changes not equal to positions: {}, {}, {}'.format(tid, aa_positions1, aa_changes))

            ref_aas = ann.get_ref_aas(tid)
            if len(ref_aas) != len(aa_positions1):
                raise Exception('ref_aas not equal to positions: {}, {}, {}'.format(tid, aa_positions1, ref_aas))

            idx = 0
            for aa_pos1 in aa_positions1:
                aa_pos0 = aa_pos1 - 1
                edited_seq = list(aa_seq.upper())
                aa_change = aa_changes[idx]
                expected_ref = ref_aas[idx]
                idx += 1

                # Example: p.Ter603Glnext*?
                if aa_pos0 >= len(aa_seq):
                    continue
                    #raise Exception("AA position is longer than AA sequence, aa_pos0: {}, len(aa_seq): {}, TID: {}, aa_cons: {}", aa_pos0, len(aa_seq), tid, ann.get_aa_cons(tid))

                if aa_pos0 >= len(edited_seq):
                    raise Exception("AA index more than seq length: {}".format(aa_pos0))

                if edited_seq[aa_pos0] != expected_ref:
                    ss = aa_seq[max(1, aa_pos0-5):min(len(aa_seq), aa_pos0+5)]
                    print('Error: expected_ref: {}, found: {}, aa_pos1: {}, cons: {}, aa_len: {}, {}, tid: {}, ref: {}, alt: {}'.format(expected_ref, edited_seq[aa_pos0], aa_pos1, ann.get_aa_cons(tid), len(edited_seq), ss, tid, record.ref, record.alts))
                    #print(edited_seq)
                    continue
                    #raise Exception(f'AA position doesnt match expected ref, found: {edited_seq[aa_pos0]} / expected: {expected_ref}')

                # This should automatically select the correct REF/ALT, based on how the SnpEffRecord was created:
                if not aa_change:
                    print('Unable to calculate AA_change: expected_ref: {}, found: {}, aa_pos1: {}, cons: {}, aa_len: {}, tid: {}, cons: {}, ref: {}, alt: {}'.format(expected_ref, edited_seq[aa_pos0], aa_pos1, ann.get_aa_cons(tid), len(edited_seq), tid, ann.get_aa_cons(tid), record.ref, record.alts))
                    continue

                if expected_ref == aa_change:
                    #print('Synonymous variant, skipping')
                    continue

                edited_seq[aa_pos0] = aa_change.lower()

            for fv in flanking_variants:
                # TODO
                allele_idx = ann.allele_idx




            # Use one less than the window size, to ensure this position is contained in the window
            # window_size_aa = int(window_size_nt / 3)
            # region_start1 = max(1, aa_pos1 - window_size_aa)
            # region_end1 = min(len(aa_seq), aa_pos0 + window_size_aa)
            # orig_seq = aa_seq[region_start1-1:region_end1]

            ret.append(Hit(record, sample_name=source_sample, allele_idx=allele_idx))

    return ret

def process_variant_queue(variants_to_process: list[VariantRecord], source_sample: str, collected_hits: list[Hit], record_map: dict[int, list[VariantRecord]], max_position_to_process: float, transcript_cache: TranscriptCache, window_size_nt: int):
    return_variants = []
    for variant in variants_to_process:
        if variant.pos < max_position_to_process:
            collected_hits.extend(process_variant(variant, source_sample, record_map, transcript_cache, window_size_nt=window_size_nt))
        else:
            return_variants.append(variant)

    #print(f'After processing: {len(return_variants)}, collected_hits: {len(collected_hits)}')
    return return_variants

def is_passing(record):
    return not record.filter or 'PASS' in record.filter

def is_discordant_genotype(record, source_sample):
    return get_unique_alleles(record, source_sample)

def get_unique_alleles(record, source_sample):
    if not is_passing(record):
        return []

    gts = list(record.samples[source_sample]['GT'])
    gts = list(filter(lambda x: x is not None, gts))
    if not gts:
        return []

    for sample in record.samples:
        if sample == source_sample:
            continue

        if not record.samples[sample]['GT']:
            continue

        for gt in list(record.samples[sample]['GT']):
            if gt in gts:
                gts.remove(gt)

    return set(gts)


def process_vcf(vcf_file: str, source_sample: str, transcript_cache: TranscriptCache, aa_flank_window: int = 10):
    print('Iterating variants')
    with pysam.VariantFile(vcf_file) as vcf_in:
        samples = list(vcf_in.header.samples)
        if not source_sample in samples:
            raise Exception('Sample name not found: {}'.format(source_sample))

        if not 'ANN' in vcf_in.header.info:
            raise Exception('The VCF must contain the ANN field, created by SnpEff')

        if not 'PS' in vcf_in.header.formats:
            raise Exception('The VCF must contain the PS format field, created by whatshap')

        # Convert AA -> NT, with extra:
        buffer_window = aa_flank_window * 4
        nt_flank_window = aa_flank_window * 3
        n_processed = 0
        record_map = {}
        variants_to_process = []
        collected_hits = []
        current_chrom = None
        for record in vcf_in:
        #for record in vcf_in.fetch("1", 914100, 10355805):
            n_processed+=1

            if n_processed % 20000 == 0:
                print('Processed {} variants, # hits: {}'.format(n_processed, len(collected_hits)))

            # Always process queue when switching contigs:
            if current_chrom is not None and record.chrom != current_chrom:
                print(f'Chromosome change: {current_chrom} -> {record.chrom}')
                if variants_to_process:
                    map(lambda x: collected_hits.extend(process_variant(x, source_sample, record_map, transcript_cache, window_size_nt=nt_flank_window)), variants_to_process)
                    variants_to_process.clear()
                record_map.clear()

            current_chrom = record.chrom

            if not record.pos in record_map.keys():
                record_map[record.pos] = []

            record_map[record.pos].append(record)

            if is_discordant_genotype(record, source_sample=source_sample):
                variants_to_process.append(record)

            max_position_to_process = record.pos - buffer_window
            variants_to_process = process_variant_queue(variants_to_process, source_sample, collected_hits, record_map, max_position_to_process, transcript_cache, window_size_nt=nt_flank_window)

            # Prune the buffer:
            min_position = min(map(lambda x: x.pos, variants_to_process)) - buffer_window if variants_to_process else record.pos - buffer_window
            record_map = {k: v for k, v in record_map.items() if k >= min_position}

            if n_processed == 500000:
                break

        # Process anything left in the queue:
        process_variant_queue(variants_to_process, source_sample, collected_hits, record_map, math.inf, transcript_cache, window_size_nt=nt_flank_window)

        print(f'Total hits: {len(collected_hits)}')
        return collected_hits

