from pathlib import Path
import gzip
import sys
from typing import Generator

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from filter_variants import is_wt
from sequence_info import F_PRIMER, R_PRIMER, WT_SEQUENCE, WT_SEQUENCE_RC
from utils import set_output_dir

_start = F_PRIMER[:8]
_start_rc = str(Seq(_start).reverse_complement())
_end = R_PRIMER[-8:]
_end_rc = str(Seq(_end).reverse_complement())

ff_junction = f'{_end}{_start}'
fr_junction = f'{_end}{_end_rc}'
rr_junction = f'{_start_rc}{_end_rc}'
rf_junction = f'{_start_rc}{_start}'


def split_ff_junction(sequence: Seq) -> list[Seq]:
    return sequence.split(ff_junction)


def split_fr_junction(sequence: Seq) -> list[Seq]:
    fragment_1, fragment_2 = sequence.split(fr_junction)
    return fragment_1, fragment_2.reverse_complement()


def split_rr_junction(sequence: Seq) -> list[Seq]:
    fragment_1, fragment_2 = sequence.split(rr_junction)
    return fragment_1.reverse_complement(), fragment_2.reverse_complement()


def split_rf_junction(sequence: Seq) -> list[Seq]:
    fragment_1, fragment_2 = sequence.split(rf_junction)
    return fragment_1.reverse_complement(), fragment_2


def keep_split_sequence(sequence: Seq) -> bool:
    return (not is_wt(sequence)) and (len(sequence) > 50)


def forward_alignment_score(sequence: Seq) -> int:
    return pairwise2.align.localxs(
        WT_SEQUENCE, sequence, -10, -1, one_alignment_only=True, score_only=True
    )


def reverse_alignment_score(sequence: Seq) -> int:
    return pairwise2.align.localxs(
        WT_SEQUENCE_RC, sequence, -10, -1, one_alignment_only=True, score_only=True
    )


def orient(sequence: Seq) -> Seq:
    forward_score = forward_alignment_score(sequence)
    reverse_score = reverse_alignment_score(sequence)

    return sequence if forward_score >= reverse_score else sequence.reverse_complement()


def orient_sequences(
    records: Generator[SeqRecord, None, None]
) -> Generator[SeqRecord, None, None]:
    def get_junction_type(sequence: Seq) -> str:
        junctions = {
            'ff': ff_junction,
            'fr': fr_junction,
            'rr': rr_junction,
            'rf': rf_junction,
        }
        for junction_type, junction_sequence in junctions:
            if sequence.count(junction_sequence) == 1:
                return junction_type
        else:
            return 'no_junction'

    junction_splitter = {
        'ff': split_ff_junction,
        'fr': split_fr_junction,
        'rr': split_rr_junction,
        'rf': split_rf_junction,
    }

    for record in records:
        sequence = record.seq
        junction_type = get_junction_type(sequence)

        if junction_type == 'no_junction':
            yield SeqRecord(orient(sequence), id=record.id)
        else:
            fragments = [
                fragment
                for fragment in junction_splitter[junction_type](sequence)
                if keep_split_sequence(fragment)
            ]
            for i, fragment in enumerate(fragments, 1):
                yield SeqRecord(fragment, id=f'{record.id}_{i}')


def main(fastq_path: Path) -> None:
    OUTPUT_DIR = set_output_dir('oriented-data/')

    with gzip.open(fastq_path, 'rt') as handle:
        generator = tqdm(orient_sequences(SeqIO.parse(handle, format='fastq')))
        file_prefix = (
            fastq_path.name.rpartition(".")[0].rpartition(".")[0].rpartition(".")[0]
        )

        SeqIO.write(
            generator,
            handle=OUTPUT_DIR / f'{file_prefix}.oriented.fasta',
            format='fasta',
        )


if __name__ == "__main__":
    fastq_path = Path(sys.argv[1])
    main(fastq_path)
