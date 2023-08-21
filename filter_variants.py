import gzip
from pathlib import Path
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from sequence_info import REFERENCE


def main(fastq_path: Path) -> None:
    OUTPUT_DIR = Path('variant-data/')

    with gzip.open(fastq_path, 'rt') as handle:
        generator = tqdm(
            record
            for record in SeqIO.parse(handle, format='fastq')
            if not is_wt(record.seq)
        )
        file_prefix = (
            fastq_path.name.rpartition(".")[0].rpartition(".")[0].rpartition(".")[0]
        )
        SeqIO.write(
            generator,
            handle=OUTPUT_DIR / f'{file_prefix}.variants.fastq',
            format='fastq',
        )


def is_wt(sequence: Seq) -> bool:
    return str(sequence.strip('N')) in REFERENCE


def get_non_wt_records(records: list[SeqRecord]) -> SeqRecord:
    return [record for record in records if not is_wt(record.seq)]


if __name__ == "__main__":
    fastq_path = Path(sys.argv[1])
    main(fastq_path)
