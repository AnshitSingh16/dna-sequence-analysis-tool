# DNA Sequence Analysis Tool

A Python-based bioinformatics application for analysing DNA sequences through nucleotide composition, GC content, codon usage, open reading frame identification, translation and interactive visualisations.

Built with Python and Streamlit, this project aims to make fundamental DNA sequence analysis accessible through an interactive web interface.

## Features

### 1. DNA Sequence Analysis
- DNA sequence validation and cleaning.
- Sequence length calculation.
- Nucleotide composition analysis.
- GC and AT content calculation.
- Ambiguous nucleotide identification.
- DNA reverse complement generation.
- DNA-to-RNA transcription.
- Protein translation.

### 2. Advanced Bioinformatics Analysis
- Six-frame translation.
- Open Reading Frame (ORF) identification.
- Start and stop codon detection.
- Codon frequency analysis.
- Codon usage statistics.
- Sliding-window GC content analysis.

### 3. Interactive Visualisation
- Nucleotide composition charts.
- Codon frequency charts.
- GC content distribution across sequence windows.
- ORF visualisation.

### 4. FASTA Support
- Upload FASTA files.
- Parse multiple FASTA records.
- Select individual sequences for analysis.
- Validate uploaded sequences.

### 5. Sequence Management
- Analysis history.
- Save frequently used DNA sequences.
- Reanalyse saved sequences.
- Delete individual records.
- Clear history and saved sequences.

### 6. Reports and Export
Export analysis results in:
- PDF
- TXT
- CSV

PDF reports can include additional scientific analysis sections, including ORFs, codon statistics, six-frame translations and sliding-window GC analysis.

### 7. Application Settings
- Light and Dark themes.
- Default reading-frame selection.
- Preferred export format.
- Persistent application settings.

## Technology Stack

| Technology | Purpose |
|---|---|
| Python | Core programming language |
| Streamlit | Interactive web application |
| Biopython | Bioinformatics library |
| Pandas | Data processing and CSV export |
| Plotly | Interactive data visualisation |
| ReportLab | PDF report generation |
| Pytest | Automated testing |

## Project Structure

```text
dna-sequence-analysis-tool/
│
├── app.py
├── requirements.txt
├── README.md
├── .gitignore
│
├── modules/
│   ├── __init__.py
│   ├── sequence_io.py
│   ├── analysis.py
│   ├── orf_analysis.py
│   ├── codon_analysis.py
│   ├── gc_analysis.py
│   ├── visualisation.py
│   ├── pdf_report.py
│   └── storage.py
│
├── tests/
│   ├── __init__.py
│   ├── test_analysis.py
│   ├── test_sequence_io.py
│   ├── test_orf_analysis.py
│   ├── test_codon_analysis.py
│   ├── test_gc_analysis.py
│
└── data/
    ├── sample_sequences.fasta
    ├── analysis_history.json
    ├── saved_sequences.json
    └── settings.json
```

## Installation and Setup

### Prerequisites

- Python 3.10 or later.
- pip.
- Git (optional).

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/dna-sequence-analysis-tool.git
```

Navigate into the project directory:

```bash
cd dna-sequence-analysis-tool
```

### 2. Create a virtual environment

```bash
python -m venv .venv
```

### 3. Activate the environment

Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
```

Windows Command Prompt:

```cmd
.venv\Scripts\activate.bat
```

Linux/macOS:

```bash
source .venv/bin/activate
```

### 4. Install dependencies

```bash
pip install -r requirements.txt
```

### 5. Run the application

```bash
streamlit run app.py
```

The application will open in your default browser.

## Example DNA Sequence

You can test the application using:

```text
ATGAAATTTGGGCCCTAA
```

This sequence contains:

- 18 nucleotides.
- 38.89% GC content.
- 61.11% AT content.
- An ATG start codon.
- A TAA stop codon.
- A six-codon-length reading frame, including the terminal stop codon.

## Running Tests

The project uses Pytest for automated testing.

Run the complete test suite:

```bash
pytest
```

To display detailed test results:

```bash
pytest -v
```

## Scientific Considerations

- DNA sequences are validated using standard IUPAC nucleotide symbols.
- Ambiguous nucleotides are excluded from canonical GC and AT content calculations.
- Translation is performed using the standard genetic code.
- Incomplete trailing codons are excluded from translation.
- ORF identification considers forward and reverse strands.
- Codon frequency calculations exclude ambiguous codons.
- Sliding-window GC analysis uses configurable window and step sizes.

Results should be interpreted in the context of the biological sequence and the assumptions of each analysis method.

## Future Development

Potential improvements include:

- DNA sequence alignment.
- BLAST integration.
- Genome annotation.
- Sequence similarity analysis.
- Mutation and variant analysis.
- Additional biological visualisations.
- Expanded testing with real genomic datasets.

## Disclaimer

This application is an educational and research-oriented bioinformatics project. Its outputs should not be treated as clinical, diagnostic or medical advice.

## Author

Developed as an independent bioinformatics and computational biology project.

**Domain:** Bioinformatics | Computational Biology | DNA Sequence Analysis