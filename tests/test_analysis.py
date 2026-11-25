from modules.analysis import (
    calculate_sequence_length,
    calculate_nucleotide_composition,
    calculate_gc_content,
    calculate_at_content,
    generate_reverse_complement,
    transcribe_dna,
    translate_dna
)


def test_sequence_length():
    assert calculate_sequence_length("ATGC") == 4


def test_nucleotide_composition():
    result = calculate_nucleotide_composition("AATTGGCC")

    assert result == {
        "A": 2,
        "T": 2,
        "G": 2,
        "C": 2
    }


def test_gc_content():
    assert calculate_gc_content("ATGCGC") == 66.67


def test_at_content():
    assert calculate_at_content("ATGCGC") == 33.33


def test_reverse_complement():
    assert generate_reverse_complement("ATGC") == "GCAT"


def test_transcription():
    assert transcribe_dna("ATGC") == "AUGC"


def test_translation():
    assert translate_dna("ATGGCC") == "MA"