import pytest

from modules.orf_analysis import (
    find_orfs_in_frame,
    find_all_orfs,
    calculate_orf_statistics
)


def test_find_orf_in_forward_strand():
    sequence = "ATGAAATAA"

    orfs = find_orfs_in_frame(sequence, 1, "+")

    assert len(orfs) == 1
    assert orfs[0]["start_position"] == 1
    assert orfs[0]["end_position"] == 9
    assert orfs[0]["protein_sequence"] == "MK"
    assert orfs[0]["length"] == 9


def test_find_orf_in_reverse_strand():
    sequence = "TTATTTCAT"

    orfs = find_orfs_in_frame(sequence, 1, "-")

    assert len(orfs) == 1
    assert orfs[0]["strand"] == "-"
    assert orfs[0]["start_position"] == 1
    assert orfs[0]["end_position"] == 9
    assert orfs[0]["protein_sequence"] == "MK"


def test_find_overlapping_orfs():
    sequence = "ATGATGTAA"

    orfs = find_orfs_in_frame(sequence, 1, "+")

    assert len(orfs) == 2

    assert orfs[0]["start_position"] == 1
    assert orfs[1]["start_position"] == 4


def test_find_orf_in_reading_frame_2():
    sequence = "AATGAAATAA"

    orfs = find_orfs_in_frame(sequence, 2, "+")

    assert len(orfs) == 1
    assert orfs[0]["reading_frame"] == 2
    assert orfs[0]["start_position"] == 2


def test_no_orf_found():
    sequence = "ATGAAAGGG"

    orfs = find_all_orfs(sequence)

    assert orfs == []


def test_invalid_reading_frame():
    with pytest.raises(ValueError):
        find_orfs_in_frame("ATGAAATAA", 4, "+")


def test_invalid_strand():
    with pytest.raises(ValueError):
        find_orfs_in_frame("ATGAAATAA", 1, "?")


def test_orf_statistics():
    orfs = [
        {"length": 9, "protein_length": 2},
        {"length": 12, "protein_length": 3}
    ]

    result = calculate_orf_statistics(orfs)

    assert result["total_orfs"] == 2
    assert result["longest_orf"] == 12
    assert result["shortest_orf"] == 9
    assert result["average_orf_length"] == 10.5
    assert result["total_protein_length"] == 5


def test_empty_orf_statistics():
    result = calculate_orf_statistics([])

    assert result["total_orfs"] == 0
    assert result["longest_orf"] == 0
    assert result["shortest_orf"] == 0