import pytest

from modules.sequence_io import (
    clean_sequence,
    validate_dna_sequence,
    parse_fasta
)


def test_clean_sequence():
    assert clean_sequence("atgc atgc") == "ATGCATGC"


def test_valid_dna_sequence():
    assert validate_dna_sequence("ATGCGCTA") == "ATGCGCTA"


def test_lowercase_dna_sequence():
    assert validate_dna_sequence("atgc") == "ATGC"


def test_invalid_dna_sequence():
    with pytest.raises(ValueError):
        validate_dna_sequence("ATGCXYZ")


def test_empty_sequence():
    with pytest.raises(ValueError):
        validate_dna_sequence("")


def test_parse_fasta():
    records = parse_fasta("data/sample_sequences.fasta")

    assert len(records) == 3
    assert records[0]["id"] == "Sequence_1"
    assert records[0]["length"] == 37