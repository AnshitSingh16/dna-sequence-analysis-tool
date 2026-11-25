import json
from pathlib import Path
from datetime import datetime, timezone
from copy import deepcopy
from uuid import uuid4
import tempfile

from modules.sequence_io import validate_dna_sequence

from modules.analysis import (
    calculate_sequence_length,
    calculate_nucleotide_composition,
    calculate_gc_content,
    calculate_at_content
)


# Project directories
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Use temporary writable storage on Streamlit Cloud.
# Use the project's data directory during local development.
if Path("/mount/src").exists():
    DATA_DIRECTORY = (
        Path(tempfile.gettempdir())
        / "dna-sequence-analysis-tool"
    )
else:
    DATA_DIRECTORY = PROJECT_ROOT / "data"

HISTORY_FILE = DATA_DIRECTORY / "analysis_history.json"
SAVED_FILE = DATA_DIRECTORY / "saved_sequences.json"
SETTINGS_FILE = DATA_DIRECTORY / "settings.json"


# Default application settings
DEFAULT_SETTINGS = {
    "theme": "Light",
    "reading_frame": 1,
    "export_format": "PDF"
}


def initialise_directory():
    """Create the data directory if it does not exist."""

    DATA_DIRECTORY.mkdir(parents=True, exist_ok=True)


def save_json(file_path, data):
    """Save data to JSON using a temporary file."""

    file_path = Path(file_path)

    file_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_file = file_path.with_suffix(
        file_path.suffix + ".tmp"
    )

    try:
        with open(temporary_file, "w", encoding="utf-8") as file:
            json.dump(
                data,
                file,
                indent=4,
                ensure_ascii=False
            )

        temporary_file.replace(file_path)

    except OSError:
        if temporary_file.exists():
            temporary_file.unlink()

        raise


def load_json(file_path, default_data):
    """Load JSON data with safe fallback handling."""

    initialise_directory()

    file_path = Path(file_path)

    if not file_path.exists():
        return deepcopy(default_data)

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)

        if not isinstance(data, type(default_data)):
            return deepcopy(default_data)

        return data

    except (json.JSONDecodeError, OSError):
        return deepcopy(default_data)


def initialise_storage():
    """Create storage files with default data when missing."""

    initialise_directory()

    storage_files = [
        (HISTORY_FILE, []),
        (SAVED_FILE, []),
        (SETTINGS_FILE, DEFAULT_SETTINGS)
    ]

    for file_path, default_data in storage_files:

        if not file_path.exists():
            save_json(file_path, default_data)


def get_history():
    """Retrieve analysis history."""

    initialise_storage()

    history = load_json(HISTORY_FILE, [])

    return [
        record for record in history
        if isinstance(record, dict)
    ]


def add_history(sequence, source="Manual Input"):
    """Add a DNA sequence analysis record to history."""

    initialise_storage()

    sequence = validate_dna_sequence(sequence)

    history = get_history()

    now = datetime.now(timezone.utc)

    record = {
        "id": str(uuid4()),
        "date": now.astimezone().isoformat(timespec="seconds"),
        "source": source,
        "sequence": sequence,
        "length": calculate_sequence_length(sequence),
        "gc_content": calculate_gc_content(sequence),
        "at_content": calculate_at_content(sequence),
        "composition": calculate_nucleotide_composition(sequence)
    }

    history.insert(0, record)

    save_json(HISTORY_FILE, history)

    return record


def delete_history(record_id):
    """Delete an individual history record."""

    history = get_history()

    history = [
        record for record in history
        if record.get("id") != record_id
    ]

    save_json(HISTORY_FILE, history)


def clear_history():
    """Delete all analysis history."""

    initialise_storage()

    save_json(HISTORY_FILE, [])


def get_saved_sequences():
    """Retrieve saved DNA sequences."""

    initialise_storage()

    saved = load_json(SAVED_FILE, [])

    return [
        record for record in saved
        if isinstance(record, dict)
    ]


def save_sequence(name, sequence):
    """Save a named DNA sequence."""

    initialise_storage()

    if not isinstance(name, str) or not name.strip():
        raise ValueError("Sequence name cannot be empty.")

    sequence = validate_dna_sequence(sequence)

    saved = get_saved_sequences()

    record = {
        "id": str(uuid4()),
        "name": name.strip(),
        "sequence": sequence,
        "length": len(sequence)
    }

    saved.insert(0, record)

    save_json(SAVED_FILE, saved)

    return record


def delete_saved_sequence(record_id):
    """Delete an individual saved sequence."""

    saved = get_saved_sequences()

    saved = [
        record for record in saved
        if record.get("id") != record_id
    ]

    save_json(SAVED_FILE, saved)


def get_settings():
    """Retrieve application settings with default values."""

    initialise_storage()

    settings = load_json(
        SETTINGS_FILE,
        DEFAULT_SETTINGS
    )

    for key, value in DEFAULT_SETTINGS.items():
        settings.setdefault(key, value)

    return settings


def update_settings(new_settings):
    """Update application settings."""

    if not isinstance(new_settings, dict):
        raise TypeError("Settings must be provided as a dictionary.")

    current_settings = get_settings()

    current_settings.update(new_settings)

    save_json(SETTINGS_FILE, current_settings)

    return current_settings