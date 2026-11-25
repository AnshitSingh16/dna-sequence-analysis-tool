from Bio import SeqIO
from pathlib import Path
from io import StringIO


# Supported IUPAC DNA nucleotide codes
VALID_DNA_BASES = set("ACGTRYSWKMBDHVN")


def clean_sequence(sequence: str) -> str:
    """Remove whitespace and convert a DNA sequence to uppercase."""

    if not isinstance(sequence, str):
        raise TypeError("DNA sequence must be a string.")

    return "".join(sequence.split()).upper()


def validate_dna_sequence(sequence: str) -> str:
    """Validate a DNA sequence and return its cleaned version."""

    sequence = clean_sequence(sequence)

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    invalid_bases = set(sequence) - VALID_DNA_BASES

    if invalid_bases:
        raise ValueError(
            f"Invalid nucleotide characters: {', '.join(sorted(invalid_bases))}"
        )

    return sequence


def parse_fasta(file_path) -> list:
    """
    Read and validate FASTA sequences.

    Supports:
    - Regular file paths
    - Streamlit uploaded files
    - Multiple FASTA records

    Returns a list of dictionaries containing:
    - id
    - description
    - sequence
    - length
    """

    # Handle Streamlit uploaded files
    if hasattr(file_path, "getvalue"):

        content = file_path.getvalue()

        if isinstance(content, bytes):
            try:
                content = content.decode("utf-8-sig")
            except UnicodeDecodeError as error:
                raise ValueError(
                    "The uploaded FASTA file must use UTF-8 encoding."
                ) from error

        elif not isinstance(content, str):
            raise TypeError("Unsupported uploaded file content.")

    # Handle regular file paths
    else:

        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(
                f"FASTA file not found: {file_path}"
            )

        if not file_path.is_file():
            raise ValueError(
                f"The provided path is not a file: {file_path}"
            )

        try:
            content = file_path.read_text(encoding="utf-8-sig")
        except UnicodeDecodeError as error:
            raise ValueError(
                "The FASTA file must use UTF-8 encoding."
            ) from error
        except OSError as error:
            raise OSError(
                f"Could not read FASTA file: {file_path}"
            ) from error

    # Check whether the file contains content
    if not content.strip():
        raise ValueError("The FASTA file is empty.")

    # Check for FASTA headers
    if not any(
        line.lstrip().startswith(">")
        for line in content.splitlines()
    ):
        raise ValueError(
            "Invalid FASTA format: no FASTA header (>) found."
        )

    # Parse FASTA records
    text_handle = StringIO(content)

    records = []

    try:
        for record in SeqIO.parse(text_handle, "fasta"):

            sequence = validate_dna_sequence(str(record.seq))

            records.append({
                "id": record.id,
                "description": record.description,
                "sequence": sequence,
                "length": len(sequence)
            })

    except ValueError:
        raise

    except Exception as error:
        raise ValueError(
            f"Could not parse FASTA file: {error}"
        ) from error

    if not records:
        raise ValueError("No valid FASTA sequences found.")

    return records