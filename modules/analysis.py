from Bio.Seq import Seq
from collections import Counter


# Standard IUPAC DNA nucleotide codes
VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")

# Canonical DNA nucleotides
CANONICAL_BASES = set("ATGC")


def clean_sequence(sequence: str) -> str:
    """Remove whitespace and convert a DNA sequence to uppercase."""
    return "".join(sequence.split()).upper()


def validate_sequence(sequence: str) -> str:
    """Validate a DNA sequence using standard IUPAC nucleotide codes."""

    sequence = clean_sequence(sequence)

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    invalid_bases = set(sequence) - VALID_DNA_BASES

    if invalid_bases:
        raise ValueError(
            f"Invalid nucleotide characters: {', '.join(sorted(invalid_bases))}"
        )

    return sequence


def calculate_sequence_length(sequence: str) -> int:
    """Calculate the total length of a DNA sequence."""

    sequence = validate_sequence(sequence)

    return len(sequence)


def calculate_nucleotide_composition(sequence: str) -> dict:
    """Calculate the frequency of canonical DNA nucleotides."""

    sequence = validate_sequence(sequence)

    counts = Counter(sequence)

    return {
        "A": counts.get("A", 0),
        "T": counts.get("T", 0),
        "G": counts.get("G", 0),
        "C": counts.get("C", 0)
    }


def calculate_ambiguous_nucleotide_count(sequence: str) -> int:
    """Count ambiguous nucleotide symbols in a DNA sequence."""

    sequence = validate_sequence(sequence)

    return sum(
        1 for nucleotide in sequence
        if nucleotide not in CANONICAL_BASES
    )


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content as a percentage.

    Ambiguous nucleotides are excluded from the denominator.
    """

    sequence = validate_sequence(sequence)

    counts = calculate_nucleotide_composition(sequence)

    canonical_count = sum(counts.values())

    if canonical_count == 0:
        return 0.0

    gc_count = counts["G"] + counts["C"]

    return round((gc_count / canonical_count) * 100, 2)


def calculate_at_content(sequence: str) -> float:
    """
    Calculate AT content as a percentage.

    Ambiguous nucleotides are excluded from the denominator.
    """

    sequence = validate_sequence(sequence)

    counts = calculate_nucleotide_composition(sequence)

    canonical_count = sum(counts.values())

    if canonical_count == 0:
        return 0.0

    at_count = counts["A"] + counts["T"]

    return round((at_count / canonical_count) * 100, 2)


def generate_reverse_complement(sequence: str) -> str:
    """Generate the reverse complement of a DNA sequence."""

    sequence = validate_sequence(sequence)

    return str(Seq(sequence).reverse_complement())


def transcribe_dna(sequence: str) -> str:
    """Transcribe a DNA sequence into RNA."""

    sequence = validate_sequence(sequence)

    return str(Seq(sequence).transcribe())


def translate_dna(sequence: str, reading_frame: int = 1) -> str:
    """
    Translate DNA into a protein sequence.

    Reading frames:
    1: Starts at nucleotide 1.
    2: Starts at nucleotide 2.
    3: Starts at nucleotide 3.

    Uses the standard genetic code.
    """

    sequence = validate_sequence(sequence)

    if (
        not isinstance(reading_frame, int)
        or isinstance(reading_frame, bool)
        or reading_frame not in (1, 2, 3)
    ):
        raise ValueError("Reading frame must be 1, 2, or 3.")

    adjusted_sequence = sequence[reading_frame - 1:]

    # Remove incomplete trailing codons
    complete_length = len(adjusted_sequence) - (len(adjusted_sequence) % 3)

    adjusted_sequence = adjusted_sequence[:complete_length]

    if not adjusted_sequence:
        return ""

    return str(Seq(adjusted_sequence).translate())


def translate_six_frames(sequence: str) -> dict:
    """
    Translate a DNA sequence in all six reading frames.

    Forward frames: +1, +2, +3
    Reverse frames: -1, -2, -3

    Uses the standard genetic code.
    """

    sequence = validate_sequence(sequence)

    reverse_sequence = str(Seq(sequence).reverse_complement())

    return {
        "+1": translate_dna(sequence, 1),
        "+2": translate_dna(sequence, 2),
        "+3": translate_dna(sequence, 3),
        "-1": translate_dna(reverse_sequence, 1),
        "-2": translate_dna(reverse_sequence, 2),
        "-3": translate_dna(reverse_sequence, 3)
    }