import pandas as pd
from sample_info import SAMPLE_NAMES
from paths import VARIANT_COUNTS_DIR


def get_indexed_wt_dict(flank: int) -> dict[int, int]:
    df = pd.read_table('indexed_wt_counts.txt').groupby('flank').get_group(flank)
    return dict(zip(df.sample_number, df.n))


def combine_sample_timepoints(start_index: int, replicate_name: str) -> pd.DataFrame:
    TIMEPOINTS = ['p0', 'p3', 'p6']
    samples = [sample for sample in SAMPLE_NAMES[start_index : start_index + 3]]
    variant_count_dfs = [
        pd.read_table(VARIANT_COUNTS_DIR / f'{sample}.counts.txt') for sample in samples
    ]

    return (
        # concatenate variant count dataframes for each replicate and assign replicate
        # name / timepoint
        pd.concat(
            [
                table.assign(sample=timepoint, replicate=replicate_name)
                for timepoint, table in zip(TIMEPOINTS, variant_count_dfs)
            ]
        )
        .pivot(
            index=['replicate', 'position', 'wt_aa', 'variant_codon', 'variant_aa'],
            columns='sample',
            values='n',
        )
        .reset_index()
        .loc[lambda x: ~x.position.isin([108, 328])]  # positions not in input library
    )


def main():
    sample_numbers = [
        (0, '1A'),
        (3, '1B'),
        (6, '2A'),
        (9, '2B'),
        (12, '2C'),
    ]

    # combine p0, p3, p6 data from each replicate
    combined_timepoints = [
        combine_sample_timepoints(sample_number, replicate)
        for sample_number, replicate in sample_numbers
    ]
    all_counts = pd.concat(combined_timepoints).rename_axis(columns=None)
    all_counts.to_csv('variant_counts.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
