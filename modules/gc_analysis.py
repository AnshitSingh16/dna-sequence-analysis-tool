from modules.sequence_io import validate_dna_sequence


def calculate_window_gc(
    sequence: str,
    window_size: int = 100,
    step_size: int = 10
) -> list:
    """
    Calculate GC content across sliding windows.

    Parameters:
        sequence: DNA sequence.
        window_size: Number of nucleotides per window.
        step_size: Number of nucleotides to move between windows.

    Returns:
        A list of dictionaries containing:
        - start_position
        - end_position
        - gc_content
        - valid_bases
        - ambiguous_bases
    """

    sequence = validate_dna_sequence(sequence)

    if (
        not isinstance(window_size, int)
        or isinstance(window_size, bool)
        or window_size <= 0
    ):
        raise ValueError("Window size must be a positive integer.")

    if (
        not isinstance(step_size, int)
        or isinstance(step_size, bool)
        or step_size <= 0
    ):
        raise ValueError("Step size must be a positive integer.")

    results = []

    # If sequence is shorter than the window,
    # analyse the complete sequence as one window.
    if len(sequence) <= window_size:
        windows = [(0, sequence)]

    else:
        windows = []

        for start in range(
            0,
            len(sequence) - window_size + 1,
            step_size
        ):
            window = sequence[start:start + window_size]
            windows.append((start, window))

    for start, window in windows:

        gc_count = window.count("G") + window.count("C")

        valid_bases = sum(
            1 for base in window
            if base in "ATGC"
        )

        ambiguous_bases = len(window) - valid_bases

        if valid_bases == 0:
            gc_content = 0.0

        else:
            gc_content = round(
                (gc_count / valid_bases) * 100,
                2
            )

        results.append({
            "start_position": start + 1,
            "end_position": start + len(window),
            "gc_content": gc_content,
            "valid_bases": valid_bases,
            "ambiguous_bases": ambiguous_bases
        })

    return results


def calculate_gc_statistics(window_results: list) -> dict:
    """Calculate summary statistics from sliding-window GC results."""

    if not window_results:
        return {
            "total_windows": 0,
            "average_gc": 0.0,
            "maximum_gc": 0.0,
            "minimum_gc": 0.0
        }

    gc_values = [
        window["gc_content"]
        for window in window_results
    ]

    return {
        "total_windows": len(window_results),
        "average_gc": round(
            sum(gc_values) / len(gc_values), 2
        ),
        "maximum_gc": max(gc_values),
        "minimum_gc": min(gc_values)
    }