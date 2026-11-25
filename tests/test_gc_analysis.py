import pytest

from modules.gc_analysis import (
    calculate_window_gc,
    calculate_gc_statistics
)


def test_basic_window_gc():
    result = calculate_window_gc(
        "ATGCGCAT",
        window_size=4,
        step_size=4
    )

    assert len(result) == 2
    assert result[0]["gc_content"] == 50.0
    assert result[1]["gc_content"] == 50.0


def test_window_positions():
    result = calculate_window_gc(
        "ATGCGCAT",
        window_size=4,
        step_size=2
    )

    assert result[0]["start_position"] == 1
    assert result[0]["end_position"] == 4
    assert result[1]["start_position"] == 3
    assert result[1]["end_position"] == 6


def test_sequence_shorter_than_window():
    result = calculate_window_gc(
        "ATGC",
        window_size=100
    )

    assert len(result) == 1
    assert result[0]["gc_content"] == 50.0


def test_ambiguous_nucleotide_handling():
    result = calculate_window_gc(
        "ATGCNN",
        window_size=6
    )

    assert result[0]["valid_bases"] == 4
    assert result[0]["ambiguous_bases"] == 2
    assert result[0]["gc_content"] == 50.0


def test_all_ambiguous_window():
    result = calculate_window_gc(
        "NNNN",
        window_size=4
    )

    assert result[0]["gc_content"] == 0.0
    assert result[0]["valid_bases"] == 0


def test_invalid_window_size():
    with pytest.raises(ValueError):
        calculate_window_gc("ATGC", window_size=0)


def test_invalid_step_size():
    with pytest.raises(ValueError):
        calculate_window_gc("ATGC", step_size=0)


def test_gc_statistics():
    results = [
        {"gc_content": 50.0},
        {"gc_content": 75.0},
        {"gc_content": 25.0}
    ]

    statistics = calculate_gc_statistics(results)

    assert statistics["total_windows"] == 3
    assert statistics["average_gc"] == 50.0
    assert statistics["maximum_gc"] == 75.0
    assert statistics["minimum_gc"] == 25.0


def test_empty_gc_statistics():
    result = calculate_gc_statistics([])

    assert result["total_windows"] == 0
    assert result["average_gc"] == 0.0