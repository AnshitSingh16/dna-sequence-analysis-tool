from collections import Counter

from modules.sequence_io import validate_dna_sequence


STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


def extract_codons(sequence: str) -> list:
    """Extract complete, unambiguous codons from a DNA sequence."""

    sequence = validate_dna_sequence(sequence)

    codons = []

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]

        if all(base in "ATGC" for base in codon):
            codons.append(codon)

    return codons


def calculate_codon_frequency(sequence: str) -> dict:
    """Calculate the frequency of each unambiguous codon."""

    codons = extract_codons(sequence)

    counts = Counter(codons)

    return dict(sorted(counts.items()))


def calculate_codon_statistics(sequence: str) -> dict:
    """
    Calculate comprehensive codon usage statistics.

    Ambiguous codons and incomplete trailing nucleotides
    are excluded from codon frequency calculations.
    """

    sequence = validate_dna_sequence(sequence)

    total_possible_codons = len(sequence) // 3

    codons = extract_codons(sequence)

    counts = Counter(codons)

    total_valid_codons = len(codons)

    frequencies = {
        codon: round(
            (count / total_valid_codons) * 100, 2
        )
        for codon, count in sorted(counts.items())
    } if total_valid_codons else {}

    most_frequent = (
        max(counts.items(), key=lambda item: item[1])[0]
        if counts else None
    )

    least_frequent = (
        min(counts.items(), key=lambda item: item[1])[0]
        if counts else None
    )

    return {
        "total_possible_codons": total_possible_codons,
        "total_valid_codons": total_valid_codons,
        "excluded_codons": total_possible_codons - total_valid_codons,
        "codon_counts": dict(sorted(counts.items())),
        "codon_percentages": frequencies,
        "start_codon_count": counts.get(START_CODON, 0),
        "stop_codon_counts": {
            codon: counts.get(codon, 0)
            for codon in sorted(STOP_CODONS)
        },
        "most_frequent_codon": most_frequent,
        "least_frequent_codon": least_frequent
    }


def calculate_codon_usage(sequence: str) -> dict:
    """Return codon usage as counts and percentages."""

    return calculate_codon_statistics(sequence)