import pandas as pd
import re
import iedb

hla_pattern = r'^HLA-([ABC]+[0-9]*)\*([0-9]{2}):([0-9]{2})$'

def predict_epitopes(df: pd.DataFrame, hla_type: list[str], source_field_name = 'Peptide', peptide_lengths_i = [8,9,10,11], peptide_lengths_ii = [10,11,12,13,14,15], max_percentile_rank = 2):
    if not source_field_name in df.columns:
        raise ValueError('Source field name not found in dataframe: {}'.format(source_field_name))

    for hla in hla_type:
        if not re.search(pattern=hla_pattern, string=hla):
            raise ValueError('{} is not a valid HLA allele name'.format(hla))

    unique_peptides = pd.DataFrame({
        'Peptide': list(set(df[source_field_name]))
    })
    unique_peptides['seq_num'] = range(1, len(df) + 1)

    mhci_res = iedb.query_mhci_binding(
        method = "recommended",
        sequence = unique_peptides['Peptide'],
        allele = hla_type,
        length = peptide_lengths_i,
    )

    mhci_res = pd.DataFrame(mhci_res)
    mhci_res['seq_num'] = mhci_res['seq_num'].astype(int)
    mhci_res['percentile_rank'] = mhci_res['percentile_rank'].astype(float)

    mhci_res = unique_peptides.merge(mhci_res, on='seq_num', how='left')

    print(mhci_res.info())
    mhci_res = (
        mhci_res.filter(lambda row: row['percentile_rank'] < max_percentile_rank)
    )
    print(mhci_res.info())

    return mhci_res

    # mhcii_res = iedb.query_mhcii_binding(
    #     method = "recommended",
    #     sequence = unique_peptides,
    #     allele = hla_type,
    #     length = peptide_lengths_ii,
    # )
    #
    # mhcii_res = pd.DataFrame(mhcii_res)
    # print(mhcii_res.columns)


