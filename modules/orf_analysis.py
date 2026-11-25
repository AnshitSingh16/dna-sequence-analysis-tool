from Bio.Seq import Seq

from modules.sequence_io import validate_dna_sequence


# Standard genetic code start and stop codons
START_CODON = "ATG"

STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs_in_frame(
    sequence: str,
    reading_frame: int,
    strand: str = "+"
) -> list:
    """
    Find complete ORFs in one reading frame.

    Every ATG is considered independently.
    Overlapping ORFs are therefore supported.

    Positions are reported using 1-based coordinates
    relative to the original DNA sequence.
    """

    if reading_frame not in (1, 2, 3):
        raise ValueError("Reading frame must be 1, 2, or 3.")

    if strand not in ("+", "-"):
        raise ValueError("Strand must be '+' or '-'.")

    sequence = validate_dna_sequence(sequence)

    original_length = len(sequence)

    # Reverse complement for negative strand
    if strand == "-":
        working_sequence = str(
            Seq(sequence).reverse_complement()
        )
    else:
        working_sequence = sequence

    frame_offset = reading_frame - 1

    orfs = []

    # Examine every codon in the selected reading frame
    for start in range(
        frame_offset,
        len(working_sequence) - 2,
        3
    ):

        start_codon = working_sequence[start:start + 3]

        if start_codon != START_CODON:
            continue

        # Find the first in-frame stop codon
        for stop in range(
            start + 3,
            len(working_sequence) - 2,
            3
        ):

            codon = working_sequence[stop:stop + 3]

            if codon in STOP_CODONS:

                # Include the stop codon in the ORF DNA sequence
                orf_sequence = working_sequence[start:stop + 3]

                # Translate without including the terminal stop
                protein_sequence = str(
                    Seq(orf_sequence).translate(to_stop=True)
                )

                if strand == "+":

                    start_position = start + 1
                    end_position = stop + 3

                else:

                    start_position = original_length - (stop + 2)
                    end_position = original_length - start

                orfs.append({
                    "strand": strand,
                    "reading_frame": reading_frame,
                    "frame": f"{strand}{reading_frame}",
                    "start_position": start_position,
                    "end_position": end_position,
                    "start_codon": START_CODON,
                    "stop_codon": codon,
                    "length": len(orf_sequence),
                    "protein_length": len(protein_sequence),
                    "sequence": orf_sequence,
                    "protein_sequence": protein_sequence,
                    "is_complete": True
                })

                # Continue searching for other start codons.
                # This allows overlapping ORFs.
                break

    return orfs


def find_all_orfs(sequence: str) -> list:
    """
    Find complete ORFs across all six reading frames.

    Forward strand: +1, +2, +3
    Reverse strand: -1, -2, -3

    Returns a list of ORF dictionaries.
    """

    sequence = validate_dna_sequence(sequence)

    all_orfs = []

    for strand in ("+", "-"):

        for reading_frame in (1, 2, 3):

            frame_orfs = find_orfs_in_frame(
                sequence,
                reading_frame,
                strand
            )

            all_orfs.extend(frame_orfs)

    # Sort ORFs by their position on the original sequence
    all_orfs.sort(
        key=lambda orf: (
            orf["start_position"],
            orf["end_position"]
        )
    )

    return all_orfs


def calculate_orf_statistics(orfs: list) -> dict:
    """Calculate summary statistics for detected ORFs."""

    if not orfs:
        return {
            "total_orfs": 0,
            "longest_orf": 0,
            "shortest_orf": 0,
            "average_orf_length": 0.0,
            "total_protein_length": 0
        }

    lengths = [orf["length"] for orf in orfs]

    protein_lengths = [
        orf["protein_length"] for orf in orfs
    ]

    return {
        "total_orfs": len(orfs),
        "longest_orf": max(lengths),
        "shortest_orf": min(lengths),
        "average_orf_length": round(
            sum(lengths) / len(lengths), 2
        ),
        "total_protein_length": sum(protein_lengths)
    }