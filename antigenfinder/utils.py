import Bio.Data.CodonTable

def aa_to_nt(aa_seq: str) -> str:
    bt = Bio.Data.CodonTable.standard_dna_table.back_table

    nt = []
    aa_seq = list(aa_seq.upper())
    for aa in aa_seq:
        if aa == '-':
            nt.append('---')
        elif aa == '~':
            nt.append('~~~')
        elif aa in bt.keys():
            nt.append(bt[aa])
        else:
            raise Exception('Unknown AA: {}'.format(aa))

    # TODO: consider including this:
    #cai_obj = CodonUsage.CodonAdaptationIndex()

    return ''.join(nt)