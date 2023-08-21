import gzip

from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd

from sequence_info import DMS_REGION, VECTOR_BACKBONE
from sample_info import SAMPLE_NAMES
from paths import ORIENTED_DATA_DIR


FLANKING_BASES = 15
FORWARD_INDEX = DMS_REGION[-FLANKING_BASES:] + 'CT'
REVERSE_INDEX = 'CT' + VECTOR_BACKBONE[:FLANKING_BASES]


def is_indexed_wt(sequence: Seq) -> bool:
    return (FORWARD_INDEX in sequence) or (REVERSE_INDEX in sequence)


def count_wt_reads(sample: str) -> int:
    print(f'counting wild-type reads from {sample}')
    with gzip.open(ORIENTED_DATA_DIR / f'{sample}.oriented.fasta.gz', 'rt') as handle:
        sequences = (record.seq for record in SeqIO.parse(handle, 'fasta'))
        return sum(is_indexed_wt(sequence) for sequence in sequences)


def main():
    wt_counts = {sample: count_wt_reads(sample) for sample in SAMPLE_NAMES}
    wt_counts_table = pd.DataFrame(
        ((key.partition('_')[0], value) for key, value in wt_counts.items()),
        columns=['sample_number', 'n'],
    )
    wt_counts_table.to_csv('indexed_wt_counts.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
