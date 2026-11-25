import pytest

from modules.codon_analysis import (
    extract_codons,
    calculate_codon_frequency,
    calculate_codon_statistics,
    calculate_codon_usage
)


def test_extract_codons():
    assert extract_codons("ATGAAATAA") == [
        "ATG", "AAA", "TAA"
    ]


def test_codon_frequency():
    result = calculate_codon_frequency("ATGATGAAATAA")

    assert result == {
        "AAA": 1,
        "ATG": 2,
        "TAA": 1
    }


def test_codon_statistics():
    result = calculate_codon_statistics("ATGAAATAA")

    assert result["total_possible_codons"] == 3
    assert result["total_valid_codons"] == 3
    assert result["excluded_codons"] == 0
    assert result["start_codon_count"] == 1

    assert result["stop_codon_counts"] == {
        "TAA": 1,
        "TAG": 0,
        "TGA": 0
    }


def test_ambiguous_codon_exclusion():
    result = calculate_codon_statistics("ATGNNNTAA")

    assert result["total_possible_codons"] == 3
    assert result["total_valid_codons"] == 2
    assert result["excluded_codons"] == 1


def test_incomplete_trailing_bases():
    result = calculate_codon_statistics("ATGAAAT")

    assert result["total_possible_codons"] == 2
    assert result["total_valid_codons"] == 2


def test_empty_frequency_result():
    result = calculate_codon_frequency("NNNNNN")

    assert result == {}


def test_empty_statistics():
    result = calculate_codon_statistics("NNNNNN")

    assert result["total_valid_codons"] == 0
    assert result["most_frequent_codon"] is None
    assert result["least_frequent_codon"] is None


def test_codon_usage():
    result = calculate_codon_usage("ATGAAATAA")

    assert result["start_codon_count"] == 1


def test_invalid_sequence():
    with pytest.raises(ValueError):
        calculate_codon_statistics("ATGXYZ")