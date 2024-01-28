import pandas as pd


def main():
    POSITION_OFFSET = 44
    
    combined_fitness_data = (
        pd.read_table('dimsum/ndufaf6_dms/fitness_singles.txt', sep=' ')[
            ['Pos', 'WT_AA', 'Mut', 'fitness', 'sigma']
        ]
        .rename(columns={'Pos': 'position', 'WT_AA': 'wt_aa', 'Mut': 'variant_aa'})
        .assign(position=lambda x: x.position + POSITION_OFFSET)
        .sort_values('position')
        .reset_index(drop=True)
        .loc[lambda x: x.position != 116]
    )

    combined_fitness_data.to_csv('data/dms_fitness.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
