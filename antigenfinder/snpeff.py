import re

from Bio.SeqUtils import seq1
from pysam import VariantRecord


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


# Converts SnpEff's 3-digit code to one letter:
def convert_codon(aa):
    if aa == '*':
        return aa
    elif aa == 'Ter':
        return '*'
    elif aa == 'ext':
        return '?'
    elif aa == 'del':
        return '-'
    elif aa == 'Ins':
        return ' '
    elif aa == 'fs':
        return '~'

    return seq1(aa)


def convert_aa(aa):
    output = []
    while aa:
        if aa.startswith('fs'):
            output.append('~')
            aa = aa[2:]
        elif aa.startswith('ins'):
            # I think we can just ignore this...
            aa = aa[3:]
        elif aa.startswith('*'):
            # I think we can just ignore this...
            aa = aa[1:]
        elif aa.startswith('?'):
            # This indicates unknown ending, like a frameshift
            aa = aa[1:]
        else:
            codon = aa[:3]
            output.append(convert_codon(codon))
            aa = aa[3:]

    return ''.join(output)

class AaConsequence:
    def __init__(self):
        self.ref = None
        self.pos = None
        self.alt = None

def parse_compound_consequence(match: re.Match, raw_cons: str) -> list[AaConsequence]:
    positions = []

    ref1 = match.group(1)
    pos_start = int(match.group(2))
    ref2 = match.group(3)
    if len(ref2) != 3:
        print('Unexpected ref2: {}, {}'.format(ref1, raw_cons))
    pos_end = int(match.group(4))
    alt = match.group(5)

    cons = AaConsequence()
    cons.ref = convert_aa(ref1)
    cons.pos = pos_start
    cons.alt = convert_aa(alt)
    positions.append(cons)

    # NOTE: leave ALT blank, since we will attribute the entire change to the first position
    current_pos = pos_start + 1
    while current_pos <= pos_end:
        cons_new = AaConsequence()
        cons_new.ref = convert_aa(ref2) if current_pos == pos_end else '.'
        cons_new.pos = current_pos
        cons_new.alt = '-'
        positions.append(cons_new)
        current_pos += 1

    expected_length = pos_end - pos_start + 1
    if len(positions) != expected_length:
        print('Incorrect complex consequence output: {}, {}'.format(expected_length, len(positions)))

    return positions

def parse_simple_consequence(match: re.Match) -> list[AaConsequence]:
    ref = match.group(1)
    pos = match.group(2)
    alt = match.group(3)

    cons = AaConsequence()
    cons.ref = convert_aa(ref)
    cons.pos = int(pos)
    cons.alt = convert_aa(alt)

    return [cons]

def parse_consequence(val: str) -> list[AaConsequence]:
    val = re.sub(r'^p\.', '', val)

    match_simple = re.search(r'^([a-zA-Z?\\*]+)([0-9]+)([a-zA-Z?\\*]+)$', val)
    if match_simple:
        return parse_simple_consequence(match_simple)

    match_complex = re.search(r'^([a-zA-Z?\\*]+)([0-9]+)_([a-zA-Z?\\*]+)([0-9]+)([a-zA-Z?\\*]+)([\\?]{0,1})$', val)
    if match_complex:
        return parse_compound_consequence(match_complex, val)

    print('No match: {}'.format(val))

    return []


class SnpEffRecord:
    def __init__(self, raw_ann: str, variant: VariantRecord):
        self.variant = variant
        tokens = raw_ann.split('|')
        self.allele = tokens[0]
        self.impact = tokens[1]
        self.effect = tokens[2]
        self.gene = tokens[3]
        self.gene_id = tokens[4]
        self.transcript_id = tokens[6]
        self.nt_consequence = tokens[9]
        self.aa_consequence = tokens[10]

        self.aa_positions: list[int] = []
        self.aa_ref: list[str] = []
        self.aa_alt: list[str] = []

        self.nt_pos: int|None = None
        self.nt_ref: str|None = None
        self.nt_alt: str|None = None


    def parse_aa_consequence(self):
        if not self.aa_consequence:
            return

        consequences = parse_consequence(self.aa_consequence)

        self.aa_ref = []
        self.aa_alt = []
        self.aa_positions = []
        for consequence in consequences:
            self.aa_ref.append(consequence.ref)
            self.aa_positions.append(consequence.pos)
            self.aa_alt.append(consequence.alt)

    def get_aa_positions(self) -> list[int]:
        if not self.aa_consequence:
            return []

        if not self.aa_positions:
            self.parse_aa_consequence()

        return self.aa_positions

    def get_alt_aa(self) -> list[str]:
        if not self.aa_consequence:
            return []

        if not self.aa_alt:
            self.parse_aa_consequence()

        return self.aa_alt

    def get_ref_aa(self) -> list[str]:
        if not self.aa_consequence:
            return []

        if not self.aa_ref:
            self.parse_aa_consequence()

        return self.aa_ref


def make_unique(list_of_lists):
    o = []
    for e in list_of_lists:
        if e not in o:
            o.append(e)

    return o


class SnpEffAnn:
    def __init__(self, record: VariantRecord, allele_idx: int):
        self.allele_idx = allele_idx
        self.allele = record.alleles[allele_idx]

        self.annotations = {}
        self.alleles = record.alleles
        for raw_ann in record.info.get('ANN'):
            x = SnpEffRecord(raw_ann, record)
            if not x.allele in self.annotations.keys():
                self.annotations[x.allele] = []

            self.annotations[x.allele].append(x)

    def get_annotations(self, tid = None) -> list[SnpEffRecord]:
        annotations = []
        if self.allele_idx == 0:
            for x in self.annotations.keys():
                annotations.extend(self.annotations[x])
        else:
            annotations.extend(self.annotations[self.allele])

        if tid:
            annotations = [x for x in annotations if x.transcript_id == tid]
            if len(annotations) == 0:
                raise Exception('No annotation found for tid for filter {}'.format(tid))

        return annotations

    def get_transcript_ids(self):
        tids = []
        for ann in self.get_annotations():
            if ann.aa_consequence:
                tids.append(ann.transcript_id)

        tids = list(set(tids))
        if len(tids) == 0:
            return []

        return tids

    def get_aa_cons(self, tid: str|None):
        ret = []
        for ann in self.get_annotations(tid):
            ret.append(ann.aa_consequence)

        return ','.join(ret)

    def get_gene_name(self, tid: str):
        ret = []

        for ann in self.get_annotations(tid):
            ret.append(ann.gene)

        ret = list(set(ret))

        return ','.join(ret)

    def get_aa_positions(self, tid: str) -> list[int]:
        ret = []
        anns = self.get_annotations(tid)
        for ann in anns:
            aa = ann.get_aa_positions()
            if aa:
                ret.append(aa)

        if len(ret) == 1:
            return ret[0]

        ret = make_unique(ret)
        if len(ret) == 0:
            # Only report this if there are annotations for this transcript:
            return []
        elif len(ret) > 1:
            print('Multiple AA positions found: {}, {}, {}, {}'.format(tid, ret, self.get_aa_cons(tid), self.get_aa_cons(None)))

        return ret[0]

    def get_ref_aas(self, tid: str) -> list[str]:
        ret = []
        for ann in self.get_annotations(tid):
            aa = ann.get_ref_aa()
            if aa:
                ret.append(aa)

        ret = make_unique(ret)
        if len(ret) == 0:
            print('No REF AA positions found: {}, {}'.format(tid, ret))
            return []
        elif len(ret) > 1:
            print('Multiple REF AA positions found: {}, {}'.format(tid, ret))

        return ret[0]

    def is_synonymous(self, tid: str) -> bool:
        idx = 0
        for ann in self.get_annotations(tid):
            aa1 = ann.get_ref_aa()
            aa2 = ann.get_alt_aa()

            if aa1 != aa2:
                return False
            else:
                idx += 1

        return idx > 0

    def get_aa_changes(self, tid: str) -> list[str]:
        ret = []
        for ann in self.get_annotations(tid):
            aa = ann.get_ref_aa() if self.allele_idx == 0 else ann.get_alt_aa()
            if aa:
                ret.append(aa)

        ret = make_unique(ret)
        if len(ret) == 0:
            print('No AA changes found: {}, {}'.format(tid, ret))
            return []
        elif len(ret) > 1:
            print('Multiple AA changes found: {}, {}'.format(tid, ret))

        return ret[0]

