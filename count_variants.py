from collections import Counter
import gzip
import itertools
from pathlib import Path
import re
import sys
from typing import Generator, Pattern

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from tqdm import tqdm

from sequence_info import DMS_REGION, VECTOR_BACKBONE
from utils import set_output_dir

variant_match = tuple[int, str]


def get_match_sequences(flank: int) -> tuple[str, str, str]:
    def _get_sequences(position: int, flank: int):
        bp_position = (position - 36) * 3
        upstream_seq = TEMPLATE[bp_position - flank : bp_position]
        wt_codon = TEMPLATE[bp_position : bp_position + 3]
        downstream_seq = TEMPLATE[bp_position + 3 : bp_position + flank + 3]
        return upstream_seq, wt_codon, downstream_seq

    TEMPLATE = DMS_REGION + VECTOR_BACKBONE
    return {position: _get_sequences(position, flank) for position in range(45, 334)}


def get_regex_patterns(match_sequences) -> list[tuple[Pattern, Pattern]]:
    return {
        position: (
            re.compile(rf'{sequence[0]}(\w{{3}})'),
            re.compile(rf'(\w{{3}}){sequence[2]}'),
        )
        for position, sequence in match_sequences.items()
    }


def get_variants(records: list[SeqRecord], flank: int) -> list[variant_match]:
    def _find_variant(seq: str, position: int) -> variant_match:
        upstream_pattern, downstream_pattern = REGEX_PATTERNS[position]
        match = re.search(upstream_pattern, seq) or re.search(downstream_pattern, seq)
        if match:
            return position, match.group(1)

    def _scan(seq: str) -> list[variant_match]:
        variants = []
        for position, sequences in MATCH_SEQUENCES.items():
            upstream_seq, wt_codon, downstream_seq = sequences

            # if query sequence is not in read, skip to next position
            if (upstream_seq not in seq) and (downstream_seq not in seq):
                continue

            # if codon at position matches wt codon, skip to next position
            if (f'{upstream_seq}{wt_codon}' in seq) or (
                f'{wt_codon}{downstream_seq}' in seq
            ):
                continue

            variant = _find_variant(seq, position)
            if variant:
                variants.append(variant)

        return variants

    def _iterate_sequence(records: Generator[SeqRecord, None, None]) -> str:
        for record in records:
            sequence = str(record.seq)
            yield list(_scan(sequence))

    MATCH_SEQUENCES = get_match_sequences(flank)
    REGEX_PATTERNS = get_regex_patterns(MATCH_SEQUENCES)

    return [variant_matches for variant_matches in tqdm(_iterate_sequence(records))]


def get_match_stats(variants: list[variant_match]) -> pd.DataFrame:
    match_counter = Counter(len(matches) for matches in variants)
    return pd.DataFrame(match_counter.items(), columns=['matches', 'n']).sort_values(
        'matches'
    )


def get_allowed_variant_codons() -> set[str]:
    codon_table = pd.read_table('variant_library_codons.txt')
    return {
        position: set(table.variant_codon)
        for position, table in codon_table.groupby('aa_position')
    }


def filter_allowed_variants(
    variant_matches: list[variant_match],
) -> list[variant_match]:
    allowed_codons = get_allowed_variant_codons()
    return [
        [
            (position, match)
            for position, match in matches
            if match in allowed_codons[position]
        ]
        for matches in variant_matches
    ]


def filter_single_variants(variant_matches: list[variant_match]) -> list[variant_match]:
    return [matches for matches in variant_matches if len(matches) == 1]


def filter_null_variants(variant_matches: list[variant_match]) -> list[variant_match]:
    return [matches for matches in variant_matches if matches]


def translate_codon(codon: str) -> str:
    return str(Seq(codon).translate())


def get_wt_aa(position: int) -> str:
    WT_SEQ = str(Seq(DMS_REGION).translate())
    return WT_SEQ[position - 36]


def make_count_table(variant_matches: list[variant_match]) -> pd.DataFrame:
    variant_counts = Counter(itertools.chain(*variant_matches))
    return pd.DataFrame(
        [(*key, value) for key, value in variant_counts.items()],
        columns=['position', 'variant_codon', 'n'],
    ).assign(
        variant_aa=lambda df: df.variant_codon.map(translate_codon),
        wt_aa=lambda df: df.position.map(get_wt_aa),
    )[
        ['position', 'wt_aa', 'variant_codon', 'variant_aa', 'n']
    ]


def main(fastq_path: Path) -> None:
    def _save_table(df: pd.DataFrame, file_suffix: str) -> None:
        df.to_csv(OUTPUT_DIR / f'{FILE_PREFIX}.{file_suffix}', sep='\t', index=False)

    OUTPUT_DIR = set_output_dir('variant-counts/')
    FILE_PREFIX = '.'.join(fastq_path.name.split('.')[:-3])

    with gzip.open(fastq_path, 'rt') as handle:
        records = SeqIO.parse(handle, 'fasta')
        variant_matches = get_variants(records, flank=15)

    single_matches = filter_single_variants(variant_matches)
    allowed_matches = filter_allowed_variants(single_matches)
    filtered_variant_matches = filter_null_variants(allowed_matches)

    unfiltered_match_stats = get_match_stats(variant_matches)
    filtered_match_stats = get_match_stats(filtered_variant_matches)

    variant_table = make_count_table(filtered_variant_matches)

    unfiltered_match_stats.pipe(_save_table, 'unfiltered_stats.txt')
    filtered_match_stats.pipe(_save_table, 'filtered_stats.txt')
    variant_table.pipe(_save_table, 'counts.txt')


if __name__ == "__main__":
    fastq_path = Path(sys.argv[1])
    main(fastq_path)
