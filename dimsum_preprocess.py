import pandas as pd
from sequence_info import DMS_REGION

variant_counts = pd.read_table('variant_counts.txt')

POSITION_OFFSET = 45
WT_SEQ = DMS_REGION[27:]


def generate_full_read(row: pd.Series) -> str:
    '''Reconstruct full-length AF6 reads from variant count data.'''
    position = (row['position'] - POSITION_OFFSET) * 3
    variant_codon = row['variant_codon']
    return f'{WT_SEQ[:position]}{variant_codon}{WT_SEQ[position + 3:]}'


def generate_variant_read_counts() -> pd.DataFrame:
    variant_count_table = (
        variant_counts.assign(nt_seq=lambda x: x.apply(generate_full_read, axis=1))[
            ['nt_seq', 'replicate', 'p0', 'p3', 'p6']
        ]
        .melt(id_vars=['nt_seq', 'replicate'], var_name='sample_id', value_name='n')
        .assign(
            sample_id=lambda x: x.sample_id.map(
                {'p0': 'input1', 'p3': 'output1', 'p6': 'output2'}
            ),
            replicate=lambda x: x.replicate.str.lower(),
        )
        .assign(
            sample_name=lambda x: x.sample_id.str[:-1]
            + x.replicate
            + x.sample_id.str[-1:]
        )
        .pivot(index='nt_seq', columns='sample_name', values='n')
        .reset_index()
        .rename_axis(columns=None)[
            [
                'nt_seq',
                'input1a1',
                'output1a1',
                'output1a2',
                'input1b1',
                'output1b1',
                'output1b2',
                'input2a1',
                'output2a1',
                'output2a2',
                'input2b1',
                'output2b1',
                'output2b2',
                'input2c1',
                'output2c1',
                'output2c2',
            ]
        ]
    )
    return variant_count_table


def generate_wt_read_counts() -> pd.DataFrame:
    wt_counts = pd.read_table('indexed_wt_counts.txt')
    wt_count_table = (
        wt_counts[['sample_number', 'n']]
        .assign(
            nt_seq=WT_SEQ,
            sample_name=lambda x: x.sample_number.map(
                {
                    1: 'input1a1',
                    2: 'output1a1',
                    3: 'output1a2',
                    4: 'input1b1',
                    5: 'output1b1',
                    6: 'output1b2',
                    7: 'input2a1',
                    8: 'output2a1',
                    9: 'output2a2',
                    10: 'input2b1',
                    11: 'output2b1',
                    12: 'output2b2',
                    13: 'input2c1',
                    14: 'output2c1',
                    15: 'output2c2',
                }
            ),
        )
        .dropna()
        .drop(columns='sample_number')
        .pivot(index='nt_seq', columns='sample_name', values='n')
        .rename_axis(columns=None)
        .reset_index()
    )
    return wt_count_table


def main():
    variant_table = generate_variant_read_counts()
    wt_table = generate_wt_read_counts()

    count_table = (
        pd.concat([variant_table, wt_table])
        .set_index('nt_seq')
        .fillna(1)
        .astype(int)
        .reset_index()
    )
    count_table.to_csv('dimsum/variant_count_table.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
