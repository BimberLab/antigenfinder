import math
import csv
import pysam
from pysam import VariantRecord

from antigenfinder import utils
from antigenfinder.prepare_gtf import TranscriptCache
from snpeff import SnpEffAnn

def fix_gt(gt: tuple) -> tuple[str, ...]:
    return tuple("." if item is None else str(item) for item in gt)

class Hit:
    def __init__(self, record: VariantRecord, sample_name: str, allele_idx: int, haplotype_idx: int, tid: str, aa_length: int, gene_name: str, aa_positions: list[int]):
        self.sample_name = sample_name
        self.chrom = record.chrom
        self.pos = record.pos
        self.ref = record.ref
        self.alt = ';'.join([str(item) for item in record.alts])
        self.allele_idx = allele_idx
        self.haplotype_idx = haplotype_idx
        self.tid = tid
        self.gene_name = gene_name
        self.aa_length = aa_length
        self.immunogenic_allele = record.alleles[allele_idx]
        self.aa_positions = aa_positions
        self.sample_gt = '/'.join(fix_gt(record.samples[self.sample_name].alleles))

        other_sample_names = [val for val in record.samples.keys() if val != sample_name]
        self.other_gts = ';'.join(map(lambda x: '/'.join(fix_gt(record.samples[x].alleles)), other_sample_names))

        self.wt_region: str|None = None
        self.seq_region: str|None = None
        self.region_start1: int|None = None
        self.region_end1: int|None = None

        self.messages: list[str] = []
        self.consequences_applied: list[str] = []

def update_seq(tid: str, record: VariantRecord, ann: SnpEffAnn, aa_positions1, aa_seq: list[str], tracker: Hit):
    aa_changes = ann.get_aa_changes(tid)
    if len(aa_changes) != len(aa_positions1):
        raise Exception('AA changes not equal to positions: {}, {}, {}'.format(tid, aa_positions1, aa_changes))

    ref_aas = ann.get_ref_aas(tid)
    if len(ref_aas) != len(aa_positions1):
        raise Exception('ref_aas not equal to positions: {}, {}, {}'.format(tid, aa_positions1, ref_aas))

    idx = 0
    orig_messages = len(tracker.messages)
    orig_consequences = len(tracker.consequences_applied)
    for aa_pos1 in aa_positions1:
        aa_pos0 = aa_pos1 - 1
        aa_change = aa_changes[idx]
        expected_ref = ref_aas[idx]
        idx += 1

        # Example: p.Ter603Glnext*?
        if aa_pos0 >= len(aa_seq):
            continue
            #raise Exception("AA position is longer than AA sequence, aa_pos0: {}, len(aa_seq): {}, TID: {}, aa_cons: {}", aa_pos0, len(aa_seq), tid, ann.get_aa_cons(tid))

        if aa_pos0 >= len(aa_seq):
            raise Exception("AA index more than seq length: {}".format(aa_pos0))

        # This will occur if there is an indel. Rather than look up the true REF, rely on '.' to indicate the WT allele:
        if expected_ref == '.':
            expected_ref = aa_seq[aa_pos0]

        # Same idea as above
        if ann.allele_idx == 0 and aa_change == '.':
            aa_change = aa_seq[aa_pos0]

        if aa_seq[aa_pos0] != expected_ref:
            if aa_seq[aa_pos0].islower():
                tracker.messages.append('Position already edited: AA Pos: {}; Variant POS: {}; Expected REF: {}; Found: {}'.format(aa_pos1, record.pos, expected_ref, aa_seq[aa_pos0]))
            else:
                ss = aa_seq[max(1, aa_pos0-5):min(len(aa_seq), aa_pos0+5)]
                raise Exception('Error: incorrect reference AA! expected_ref: {}, found: {}, aa_pos1: {}, cons: {}, aa_len: {}, {}, tid: {}, ref: {}, alt: {}'.format(expected_ref, aa_seq[aa_pos0], aa_pos1, ann.get_aa_cons(tid), len(aa_seq), ss, tid, record.ref, record.alts))

        # This should automatically select the correct REF/ALT, based on how the SnpEffRecord was created:
        if not aa_change:
            tracker.messages.append('Unable to calculate AA_change: ref: {}; aa_pos1: {}; cons: {}; ref: {}; alt: {}'.format(aa_seq[aa_pos0], aa_pos1, ann.get_aa_cons(tid), record.ref, record.alts))
            continue

        # Skip synonymous changes:
        if expected_ref == aa_change:
            continue

        tracker.consequences_applied.append('NT-{}: {}'.format(record.pos, ann.get_aa_cons(tid)))
        aa_seq[aa_pos0] = aa_change.lower()

    return orig_consequences != len(tracker.consequences_applied) or orig_messages != len(tracker.messages)

def process_variant(record: VariantRecord, source_sample: str, record_buffer: dict[int, list[VariantRecord]], transcript_cache: TranscriptCache, window_size_nt: int) -> list[Hit]:
    if not record.info.get('ANN'):
        raise Exception('No ANN value for variant: {}'.format(record))

    ret = []
    for allele_idx in get_unique_alleles(record, source_sample):
        if allele_idx is None:
            print('Allele IDX is empty: {}'.format(allele_idx))
            print(record.samples[source_sample]['GT'])
            print(get_unique_alleles(record, source_sample))

        ann = SnpEffAnn(record, allele_idx)
        tids = ann.get_transcript_ids()
        if not tids:
            continue

        flanking_variants: list[VariantRecord] = []
        ps = record.samples[source_sample].get('PS')
        if ps:
            min_nt_pos = record.pos - window_size_nt
            max_nt_pos = record.pos + window_size_nt
            for nt_pos in record_buffer.keys():
                if min_nt_pos <= nt_pos <= max_nt_pos:
                    for fv in record_buffer[nt_pos]:
                        if is_passing(fv) and fv.pos != record.pos and fv.samples[source_sample].get('PS') == ps:
                            flanking_variants.append(fv)

        for tid in tids:
            aa_seq = transcript_cache.find_transcript_sequence(tid)
            aa_positions1 = ann.get_aa_positions(tid)
            if not aa_positions1:
                raise Exception('No passing AA positions: tid: {}, pos: {}'.format(tid, record.pos))

            gene_name = ann.get_gene_name(tid)

            # Determine which haplotypes to process:
            gts = record.samples[source_sample].get('GT')
            indices = [i for i, val in enumerate(gts) if val == allele_idx]
            for haplotype_idx in indices:
                hit = Hit(
                    tid = tid,
                    gene_name = gene_name,
                    aa_length = len(aa_seq),
                    record=record,
                    sample_name=source_sample,
                    allele_idx=allele_idx,
                    aa_positions=aa_positions1,
                    haplotype_idx = haplotype_idx + 1
                )

                edited_seq = list(aa_seq.upper())
                changes_made = update_seq(tid, record, ann, aa_positions1, edited_seq, hit)
                if not changes_made:
                    continue

                for fv in flanking_variants:
                    # haplotype_index determines whether we are inspecting H1 or H2.
                    # we need to inspect the GTs at this position, and translate haplotype_index to allele
                    fgt = fv.samples[source_sample].get('GT')
                    if fgt is None or None in fgt:
                        continue

                    allele_idx = fgt[haplotype_idx]

                    # WT, nothing to do:
                    if allele_idx == 0:
                        continue

                    fa = SnpEffAnn(fv, allele_idx)
                    faa = fa.get_aa_positions(tid)
                    if faa:
                        try:
                            update_seq(tid, fv, fa, faa, edited_seq, hit)
                        except Exception as e:
                            print('ERROR: pos: {}, parent_aa: {}, faa: {}, FV NT: {}'.format(record.pos, aa_positions1, faa, fv.pos))
                            raise e

                # Use one less than the window size, to ensure this position is contained in the window
                window_size_aa = int(window_size_nt / 3)
                region_start1 = max(1, min(aa_positions1) - window_size_aa)
                region_end1 = min(len(aa_seq), min(aa_positions1)-1 + window_size_aa)

                hit.region_start1 = region_start1
                hit.region_end1 = region_end1
                hit.seq_region = ''.join(edited_seq[region_start1-1:region_end1])
                hit.wt_region = aa_seq[region_start1-1:region_end1]

                ret.append(hit)
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

    total_passing = 0
    for sample in record.samples:
        if sample == source_sample:
            continue

        if record.samples[sample]['GT'] is None or any(g is None for g in record.samples[sample]['GT']):
            continue

        total_passing += 1
        for gt in list(record.samples[sample]['GT']):
            if gt in gts:
                gts = [x for x in gts if x != gt]

    if not total_passing:
        return []

    return set(gts)


def process_vcf(vcf_file: str, source_sample: str, output_file: str, transcript_cache: TranscriptCache, aa_flank_window: int = 10, max_records_to_process = -1) -> list[Hit]:
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

            if n_processed == max_records_to_process:
                break

        # Process anything left in the queue:
        process_variant_queue(variants_to_process, source_sample, collected_hits, record_map, math.inf, transcript_cache, window_size_nt=nt_flank_window)

        print(f'Total hits: {len(collected_hits)}')

        if output_file:
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(['SampleName', 'TranscriptId', 'GeneName', 'AA_Length', 'Chrom', 'NT_Pos', 'NT_Change', 'Sample_GT', 'Other_GTs', 'ImmunogenicAllele', 'Category', 'HaplotypeNumber', 'AA_Pos', 'RegionStart', 'RegionEnd', 'AA_Seq', 'WT-Region', 'AA_Changes', 'NT_Translation', 'Warnings', 'Messages'])
                for hit in collected_hits:
                    writer.writerow([
                        hit.sample_name,
                        hit.tid,
                        hit.gene_name,
                        hit.aa_length,
                        hit.chrom,
                        hit.pos,
                        hit.ref + '>' + hit.alt,
                        hit.sample_gt,
                        hit.other_gts,
                        hit.immunogenic_allele,
                        'Wild-Type' if hit.allele_idx == 0 else 'Variant',
                        hit.haplotype_idx,
                        ';'.join(str(num) for num in hit.aa_positions),
                        hit.region_start1,
                        hit.region_end1,
                        hit.seq_region,
                        hit.wt_region,
                        ';'.join(set(hit.consequences_applied)),
                        utils.aa_to_nt(hit.seq_region),
                        'REVIEW' if hit.allele_idx > 0 and hit.seq_region == hit.wt_region else '',
                        ';'.join(set(hit.messages))
                    ])
        return collected_hits

