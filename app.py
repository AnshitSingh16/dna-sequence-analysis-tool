import csv
import io
from typing import Any

import pandas as pd
import plotly.express as px
import streamlit as st

from modules.visualisation import (
    create_nucleotide_chart,
    create_codon_frequency_chart,
    create_gc_window_chart,
    create_orf_visualisation,
)
from modules.pdf_report import generate_pdf_report
from modules.sequence_io import validate_dna_sequence, parse_fasta
from modules.storage import (
    add_history,
    get_history,
    delete_history,
    clear_history,
    get_saved_sequences,
    save_sequence,
    delete_saved_sequence,
    get_settings,
    update_settings,
)
from modules.analysis import (
    calculate_sequence_length,
    calculate_nucleotide_composition,
    calculate_ambiguous_nucleotide_count,
    calculate_gc_content,
    calculate_at_content,
    generate_reverse_complement,
    transcribe_dna,
    translate_dna,
    translate_six_frames,
)
from modules.orf_analysis import find_all_orfs, calculate_orf_statistics
from modules.codon_analysis import calculate_codon_statistics
from modules.gc_analysis import calculate_window_gc, calculate_gc_statistics


# ==================================================
# PAGE CONFIGURATION
# ==================================================

st.set_page_config(
    page_title="Genome Sequence Analysis Tool",
    page_icon="🧬",
    layout="wide",
)


# ==================================================
# SESSION STATE
# ==================================================

DEFAULT_SESSION_STATE = {
    "dna_sequence": "",
    "fasta_records": [],
    "active_page": "Dashboard",
    "confirm_clear_history": False,
    "confirm_clear_saved": False,
}

for state_key, default_value in DEFAULT_SESSION_STATE.items():
    if state_key not in st.session_state:
        st.session_state[state_key] = default_value


# ==================================================
# APPLICATION SETTINGS — VALIDATED SAFELY
# ==================================================

settings = get_settings()

current_theme = settings.get("theme", "Light")
if current_theme not in ["Light", "Dark"]:
    current_theme = "Light"

reading_frame = settings.get("reading_frame", 1)
if reading_frame not in [1, 2, 3]:
    reading_frame = 1

export_format = settings.get("export_format", "PDF")
if export_format not in ["PDF", "CSV", "TXT"]:
    export_format = "PDF"


# ==================================================
# CUSTOM DESIGN — APPROVED LIGHT MODE CSS
# ==================================================

st.markdown("""
<style>

.stApp,
[data-testid="stAppViewContainer"],
[data-testid="stMain"] {
    background-color: #ffffff !important;
    color: #263445 !important;
    font-family: "Segoe UI", Arial, sans-serif;
}

[data-testid="stHeader"] {
    background-color: #ffffff !important;
}

.block-container {
    max-width: 1400px;
    padding-top: 1.5rem;
    padding-bottom: 2rem;
}

h1, h2, h3, h4 {
    color: #24364b !important;
}

p, label, li {
    color: #344256;
}

section[data-testid="stSidebar"] {
    background-color: #f4f6f9 !important;
    border-right: 1px solid #e1e6ed;
}

section[data-testid="stSidebar"] p,
section[data-testid="stSidebar"] span,
section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] div {
    color: #344256 !important;
}

section[data-testid="stSidebar"] h1,
section[data-testid="stSidebar"] h2,
section[data-testid="stSidebar"] h3 {
    color: #24364b !important;
}

section[data-testid="stSidebar"] .stButton button {
    background-color: #ffffff !important;
    color: #263445 !important;
    border: 1px solid #dce3eb !important;
    border-radius: 7px !important;
    min-height: 40px;
    text-align: left;
}

section[data-testid="stSidebar"] .stButton button p,
section[data-testid="stSidebar"] .stButton button span {
    color: #263445 !important;
}

section[data-testid="stSidebar"] .stButton button:hover {
    background-color: #eaf0f7 !important;
    border-color: #b8c8dc !important;
}

.stButton button {
    background-color: #315b85 !important;
    color: #ffffff !important;
    border: 1px solid #315b85 !important;
    border-radius: 7px !important;
    padding: 8px 16px !important;
    font-weight: 600 !important;
}

.stButton button p,
.stButton button span {
    color: #ffffff !important;
}

.stButton button:hover {
    background-color: #25496f !important;
}

.stDownloadButton button {
    background-color: #ffffff !important;
    color: #315b85 !important;
    border: 1px solid #315b85 !important;
    border-radius: 7px !important;
}

.stDownloadButton button p,
.stDownloadButton button span {
    color: #315b85 !important;
}

.stTextArea textarea,
.stTextInput input,
[data-baseweb="input"] input,
[data-baseweb="textarea"] textarea {
    background-color: #ffffff !important;
    color: #263445 !important;
    border: 1px solid #cbd5e1 !important;
    border-radius: 7px !important;
}

.stTextArea textarea::placeholder,
.stTextInput input::placeholder {
    color: #8491a3 !important;
}

[data-testid="stFileUploader"] {
    background-color: #ffffff !important;
    border-radius: 8px !important;
}

[data-testid="stFileUploader"] section {
    background-color: #f7f9fc !important;
    border: 1px dashed #b8c5d6 !important;
    border-radius: 8px !important;
}

[data-testid="stFileUploader"] p,
[data-testid="stFileUploader"] span,
[data-testid="stFileUploader"] small,
[data-testid="stFileUploader"] label,
[data-testid="stFileUploader"] div {
    color: #344256 !important;
}

[data-testid="stFileUploader"] button {
    background-color: #ffffff !important;
    border: 1px solid #cbd5e1 !important;
    color: #315b85 !important;
}

[data-testid="stFileUploader"] button p,
[data-testid="stFileUploader"] button span {
    color: #315b85 !important;
}

[data-testid="stCode"],
[data-testid="stCode"] pre,
[data-testid="stCode"] code,
.stCodeBlock,
pre,
code {
    background-color: #f1f4f8 !important;
    color: #263445 !important;
    border: 1px solid #e1e6ed !important;
    border-radius: 7px !important;
}

[data-testid="stCode"] *,
pre *,
code * {
    color: #263445 !important;
}

[data-testid="stAlert"] {
    background-color: #f1f5f9 !important;
    border: 1px solid #dce5ef !important;
    border-radius: 7px !important;
    color: #334155 !important;
}

[data-testid="stAlert"] p,
[data-testid="stAlert"] span,
[data-testid="stAlert"] div {
    background-color: transparent !important;
    color: #334155 !important;
}

[data-testid="stMetric"] {
    background-color: #f7f9fc !important;
    border: 1px solid #e2e8f0 !important;
    border-radius: 8px !important;
    padding: 12px !important;
}

[data-testid="stMetricLabel"] {
    color: #64748b !important;
}

[data-testid="stMetricValue"] {
    color: #24364b !important;
}

[data-testid="stExpander"] {
    background-color: #ffffff !important;
    border: 1px solid #e2e8f0 !important;
    border-radius: 7px !important;
}

[data-testid="stExpander"] summary {
    color: #263445 !important;
}

div[data-testid="stVerticalBlock"] {
    gap: 0.8rem;
}

hr {
    border-color: #e5eaf0 !important;
    margin-top: 0.8rem !important;
    margin-bottom: 0.8rem !important;
}

footer {
    visibility: hidden;
}

</style>
""", unsafe_allow_html=True)


# ==================================================
# DARK MODE CSS
# ==================================================

if current_theme == "Dark":
    st.markdown("""
    <style>

    .stApp,
    [data-testid="stAppViewContainer"],
    [data-testid="stMain"],
    [data-testid="stHeader"] {
        background-color: #111827 !important;
        color: #e5e7eb !important;
    }

    h1, h2, h3, h4 {
        color: #f1f5f9 !important;
    }

    p, label, li, span {
        color: #d1d5db;
    }

    section[data-testid="stSidebar"] {
        background-color: #1f2937 !important;
        border-right: 1px solid #374151 !important;
    }

    section[data-testid="stSidebar"] p,
    section[data-testid="stSidebar"] span,
    section[data-testid="stSidebar"] label,
    section[data-testid="stSidebar"] div,
    section[data-testid="stSidebar"] h1,
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3 {
        color: #e5e7eb !important;
    }

    section[data-testid="stSidebar"] .stButton button {
        background-color: #273548 !important;
        color: #f1f5f9 !important;
        border: 1px solid #475569 !important;
    }

    section[data-testid="stSidebar"] .stButton button p,
    section[data-testid="stSidebar"] .stButton button span {
        color: #f1f5f9 !important;
    }

    section[data-testid="stSidebar"] .stButton button:hover {
        background-color: #374b63 !important;
    }

    .stButton button {
        background-color: #315b85 !important;
        color: #ffffff !important;
        border: 1px solid #527ba5 !important;
    }

    .stDownloadButton button {
        background-color: #273548 !important;
        color: #93c5fd !important;
        border: 1px solid #527ba5 !important;
    }

    .stDownloadButton button p,
    .stDownloadButton button span {
        color: #93c5fd !important;
    }

    .stTextArea textarea,
    .stTextInput input,
    [data-baseweb="input"] input,
    [data-baseweb="textarea"] textarea {
        background-color: #1f2937 !important;
        color: #f1f5f9 !important;
        border: 1px solid #475569 !important;
    }

    .stTextArea textarea::placeholder,
    .stTextInput input::placeholder {
        color: #9ca3af !important;
    }

    [data-baseweb="select"] > div {
        background-color: #1f2937 !important;
        color: #f1f5f9 !important;
        border-color: #475569 !important;
    }

    [data-baseweb="select"] span {
        color: #f1f5f9 !important;
    }

    [role="radiogroup"] label {
        color: #e5e7eb !important;
    }

    [data-testid="stFileUploader"] {
        background-color: #1f2937 !important;
    }

    [data-testid="stFileUploader"] section {
        background-color: #273548 !important;
        border: 1px dashed #64748b !important;
    }

    [data-testid="stFileUploader"] p,
    [data-testid="stFileUploader"] span,
    [data-testid="stFileUploader"] small,
    [data-testid="stFileUploader"] label,
    [data-testid="stFileUploader"] div {
        color: #e5e7eb !important;
    }

    [data-testid="stFileUploader"] button {
        background-color: #1f2937 !important;
        border: 1px solid #64748b !important;
    }

    [data-testid="stCode"],
    [data-testid="stCode"] pre,
    [data-testid="stCode"] code,
    .stCodeBlock,
    pre,
    code {
        background-color: #1f2937 !important;
        color: #e5e7eb !important;
        border: 1px solid #374151 !important;
    }

    [data-testid="stCode"] *,
    pre *,
    code * {
        color: #e5e7eb !important;
    }

    [data-testid="stAlert"] {
        background-color: #273548 !important;
        border: 1px solid #475569 !important;
        color: #e5e7eb !important;
    }

    [data-testid="stAlert"] p,
    [data-testid="stAlert"] span,
    [data-testid="stAlert"] div {
        background-color: transparent !important;
        color: #e5e7eb !important;
    }

    [data-testid="stMetric"] {
        background-color: #1f2937 !important;
        border: 1px solid #374151 !important;
    }

    [data-testid="stMetricLabel"] {
        color: #9ca3af !important;
    }

    [data-testid="stMetricValue"] {
        color: #f1f5f9 !important;
    }

    [data-testid="stExpander"] {
        background-color: #1f2937 !important;
        border: 1px solid #374151 !important;
    }

    [data-testid="stExpander"] summary {
        color: #e5e7eb !important;
    }

    [data-testid="stExpander"] p,
    [data-testid="stExpander"] span {
        color: #e5e7eb !important;
    }

    hr {
        border-color: #374151 !important;
    }

    </style>
    """, unsafe_allow_html=True)


# ==================================================
# HELPER FUNCTIONS
# ==================================================

def set_sequence_and_open_dashboard(sequence: str) -> None:
    """Load a sequence into the dashboard input and open the dashboard."""
    st.session_state["dna_sequence"] = sequence
    st.session_state["pending_manual_input"] = sequence
    st.session_state["active_page"] = "Dashboard"


def safe_percent(value: Any) -> str:
    """Format percentage-like values without assuming a numeric type."""
    try:
        return f"{float(value):.2f}%"
    except (TypeError, ValueError):
        return str(value)


def dataframe_download_bytes(dataframe: pd.DataFrame) -> bytes:
    """Create a UTF-8 CSV file with correct quoting and escaping."""
    return dataframe.to_csv(index=False).encode("utf-8-sig")


def build_export_dataframe(
    sequence: str,
    length: int,
    composition: dict,
    ambiguous_count: int,
    gc_content: float,
    at_content: float,
    reverse_complement: str,
    rna_sequence: str,
    protein_sequence: str,
    reading_frame: int,
    six_frame_translations: Any,
    orfs: list,
    codon_statistics: dict,
    gc_statistics: dict,
) -> pd.DataFrame:
    """Flatten core analysis outputs into a safe, consistently quoted CSV."""
    rows = [
        ("Sequence", "Original DNA Sequence", sequence),
        ("Statistics", "Sequence Length (bp)", length),
        ("Statistics", "Ambiguous Nucleotide Count", ambiguous_count),
        ("Statistics", "GC Content (%)", gc_content),
        ("Statistics", "AT Content (%)", at_content),
        ("Composition", "Adenine (A)", composition.get("A", 0)),
        ("Composition", "Thymine (T)", composition.get("T", 0)),
        ("Composition", "Guanine (G)", composition.get("G", 0)),
        ("Composition", "Cytosine (C)", composition.get("C", 0)),
        ("Sequence", "Reverse Complement", reverse_complement),
        ("Sequence", "RNA Sequence", rna_sequence),
        ("Translation", f"Protein Translation (Frame {reading_frame})", protein_sequence),
        ("ORF Analysis", "Number of ORFs", len(orfs)),
    ]

    if isinstance(six_frame_translations, dict):
        for frame_name, translated in six_frame_translations.items():
            rows.append(("Six-Frame Translation", str(frame_name), str(translated)))
    elif isinstance(six_frame_translations, (list, tuple)):
        for index, translated in enumerate(six_frame_translations, start=1):
            rows.append(("Six-Frame Translation", f"Frame {index}", str(translated)))

    for key, value in codon_statistics.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            rows.append(("Codon Analysis", str(key), value))
        elif isinstance(value, dict):
            for subkey, subvalue in value.items():
                rows.append(("Codon Analysis", f"{key}: {subkey}", subvalue))

    for key, value in gc_statistics.items():
        if isinstance(value, (str, int, float, bool)) or value is None:
            rows.append(("Sliding-Window GC", str(key), value))

    for index, orf in enumerate(orfs, start=1):
        rows.append(("ORF", f"ORF {index} details", str(orf)))

    return pd.DataFrame(rows, columns=["Section", "Parameter", "Value"])


def build_text_report(
    sequence: str,
    length: int,
    composition: dict,
    ambiguous_count: int,
    gc_content: float,
    at_content: float,
    reverse_complement: str,
    rna_sequence: str,
    protein_sequence: str,
    reading_frame: int,
    six_frame_translations: Any,
    orfs: list,
    orf_statistics: dict,
    codon_statistics: dict,
    gc_window_results: list,
    gc_statistics: dict,
) -> str:
    """Build a readable plain-text report containing all analysis sections."""
    lines = [
        "GENOME SEQUENCE ANALYSIS TOOL",
        "DNA ANALYSIS REPORT",
        "=" * 60,
        "",
        "1. ORIGINAL DNA SEQUENCE",
        sequence,
        "",
        "2. SEQUENCE STATISTICS",
        f"Sequence Length: {length} bp",
        f"GC Content: {gc_content}%",
        f"AT Content: {at_content}%",
        f"Ambiguous Nucleotide Count: {ambiguous_count}",
        "GC/AT percentages use canonical A, T, G and C bases only.",
        "",
        "3. NUCLEOTIDE COMPOSITION",
        f"Adenine (A): {composition.get('A', 0)}",
        f"Thymine (T): {composition.get('T', 0)}",
        f"Guanine (G): {composition.get('G', 0)}",
        f"Cytosine (C): {composition.get('C', 0)}",
        "",
        "4. SEQUENCE TRANSFORMATIONS",
        "Reverse Complement:",
        reverse_complement,
        "",
        "RNA Transcription:",
        rna_sequence,
        "",
        "5. PROTEIN TRANSLATION",
        f"Selected Reading Frame: {reading_frame}",
        protein_sequence,
        "",
        "Six-Frame Translation:",
    ]

    if isinstance(six_frame_translations, dict):
        for frame_name, translated in six_frame_translations.items():
            lines.append(f"{frame_name}: {translated}")
    elif isinstance(six_frame_translations, (list, tuple)):
        for index, translated in enumerate(six_frame_translations, start=1):
            lines.append(f"Frame {index}: {translated}")
    else:
        lines.append(str(six_frame_translations))

    lines.extend([
        "",
        "6. OPEN READING FRAME (ORF) ANALYSIS",
        f"ORFs identified: {len(orfs)}",
        f"ORF statistics: {orf_statistics}",
    ])

    for index, orf in enumerate(orfs, start=1):
        lines.append(
            f"ORF {index}: strand={orf.get('strand')}, "
            f"frame={orf.get('reading_frame')}, "
            f"start={orf.get('start_position')}, "
            f"end={orf.get('end_position')}, "
            f"length={orf.get('length')} bp, "
            f"protein_length={orf.get('protein_length')} aa, "
            f"complete={orf.get('is_complete')}"
        )
        lines.append(f"DNA: {orf.get('sequence', '')}")
        lines.append(f"Protein: {orf.get('protein_sequence', '')}")

    lines.extend([
        "",
        "7. CODON ANALYSIS",
        str(codon_statistics),
        "",
        "8. SLIDING-WINDOW GC ANALYSIS",
        f"GC statistics: {gc_statistics}",
        f"Windows analysed: {len(gc_window_results)}",
        "",
    ])

    for window in gc_window_results:
        lines.append(
            f"Positions {window.get('start_position')}-"
            f"{window.get('end_position')}: "
            f"GC={window.get('gc_content')}%, "
            f"valid bases={window.get('valid_bases')}, "
            f"ambiguous bases={window.get('ambiguous_bases')}"
        )

    lines.extend([
        "",
        "=" * 60,
        "Generated by Genome Sequence Analysis Tool",
    ])

    return "\n".join(lines)


# ==================================================
# SIDEBAR NAVIGATION
# ==================================================

with st.sidebar:
    st.markdown("# 🧬 GenomeLab")
    st.caption("DNA SEQUENCE ANALYSIS")
    st.divider()
    st.markdown("### Navigation")

    navigation = [
        ("🏠", "Dashboard"),
        ("🔎", "Search"),
        ("🕘", "History"),
        ("📂", "Saved"),
        ("⚙️", "Settings"),
        ("ⓘ", "About"),
        ("📄", "Terms & Conditions"),
    ]

    for icon, label in navigation:
        if st.button(f"{icon}  {label}", key=f"nav_{label}"):
            st.session_state["active_page"] = label
            st.rerun()

    st.divider()
    st.caption("Genome Sequence Analysis Tool")
    st.caption("Scientific analysis workspace")


# ==================================================
# HEADER
# ==================================================

st.markdown(
    '<p class="main-title">Genome Sequence Analysis Tool</p>',
    unsafe_allow_html=True,
)
st.markdown(
    '<p class="subtitle">A Bioinformatics Platform for DNA Sequence Analysis</p>',
    unsafe_allow_html=True,
)
st.divider()

active_page = st.session_state["active_page"]


# ==================================================
# SEARCH
# ==================================================

if active_page == "Search":
    st.header("🔎 Search")
    st.write("Search functionality will be developed in a future update.")
    st.info(
        "This section will allow users to find sequences, "
        "analysis results and previous records."
    )


# ==================================================
# HISTORY
# ==================================================

elif active_page == "History":
    st.header("Analysis History")
    st.caption("View and manage your previous DNA sequence analyses.")

    history = get_history()

    if not history:
        st.info("No analysis history available yet.")
    else:
        metric_columns = st.columns(2)
        with metric_columns[0]:
            st.metric("Total Analyses", len(history))
        with metric_columns[1]:
            unique_sequences = len({
                record.get("sequence", "") for record in history
            })
            st.metric("Unique Sequences", unique_sequences)

        st.divider()
        col1, col2 = st.columns([3, 1])
        with col1:
            st.subheader("Previous Analyses")
        with col2:
            if st.button("Clear All History", key="clear_history"):
                st.session_state["confirm_clear_history"] = True

        if st.session_state["confirm_clear_history"]:
            st.warning(
                "Are you sure you want to delete your entire analysis history?"
            )
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Yes, Clear History", key="confirm_clear"):
                    clear_history()
                    st.session_state["confirm_clear_history"] = False
                    st.success("Analysis history cleared.")
                    st.rerun()
            with col2:
                if st.button("Cancel", key="cancel_clear"):
                    st.session_state["confirm_clear_history"] = False
                    st.rerun()

        st.divider()

        for index, record in enumerate(history):
            record_id = record.get("id", f"history_{index}")
            sequence_value = record.get("sequence", "")
            date = record.get("date", "Unknown date")
            source = record.get("source", "Unknown source")
            length_value = record.get("length", len(sequence_value))
            gc_value = record.get("gc_content", 0)

            with st.expander(
                f"Analysis {index + 1} | {length_value} bp | GC: {gc_value}%"
            ):
                st.caption(f"Date: {date}")
                st.caption(f"Source: {source}")

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Length", f"{length_value} bp")
                with col2:
                    st.metric("GC Content", f"{gc_value}%")
                with col3:
                    st.metric("AT Content", f"{record.get('at_content', 0)}%")

                st.write("**Nucleotide Composition**")
                st.json(record.get("composition", {}))
                st.write("**DNA Sequence**")
                st.code(sequence_value, language="text")

                col1, col2 = st.columns(2)
                with col1:
                    if st.button(
                        "Reanalyse Sequence",
                        key=f"reanalyse_{record_id}",
                    ):
                        set_sequence_and_open_dashboard(sequence_value)
                        st.rerun()
                with col2:
                    if st.button(
                        "Delete Record",
                        key=f"delete_{record_id}",
                    ):
                        delete_history(record_id)
                        st.rerun()


# ==================================================
# SAVED SEQUENCES
# ==================================================

elif active_page == "Saved":
    st.header("Saved Sequences")
    st.caption("Your bookmarked DNA sequences, stored for future analysis.")

    saved_sequences = get_saved_sequences()

    if not saved_sequences:
        st.info("You haven't saved any sequences yet.")
        st.write(
            "Analyse a DNA sequence on the Dashboard and save it here."
        )
    else:
        st.metric("Total Saved Sequences", len(saved_sequences))
        st.divider()

        for index, record in enumerate(saved_sequences):
            record_id = record.get("id", f"saved_{index}")
            name = record.get("name", f"Sequence {index + 1}")
            sequence_value = record.get("sequence", "")
            length_value = record.get("length", len(sequence_value))

            with st.expander(f"{name} | {length_value} bp"):
                st.write("**DNA Sequence**")
                st.code(sequence_value, language="text")

                col1, col2 = st.columns(2)
                with col1:
                    if st.button(
                        "Analyse Sequence",
                        key=f"analyse_saved_{record_id}",
                    ):
                        set_sequence_and_open_dashboard(sequence_value)
                        st.rerun()
                with col2:
                    if st.button(
                        "Delete Saved Sequence",
                        key=f"delete_saved_{record_id}",
                    ):
                        delete_saved_sequence(record_id)
                        st.rerun()


# ==================================================
# SETTINGS
# ==================================================

elif active_page == "Settings":
    st.header("⚙️ Settings")
    st.caption("Customise your Genome Sequence Analysis workspace.")

    settings_page = get_settings()

    saved_theme = settings_page.get("theme", "Light")
    if saved_theme not in ["Light", "Dark"]:
        saved_theme = "Light"

    saved_frame = settings_page.get("reading_frame", 1)
    if saved_frame not in [1, 2, 3]:
        saved_frame = 1

    saved_export = settings_page.get("export_format", "PDF")
    if saved_export not in ["PDF", "CSV", "TXT"]:
        saved_export = "PDF"

    st.subheader("Appearance")
    theme = st.radio(
        "Application Theme",
        options=["Light", "Dark"],
        index=["Light", "Dark"].index(saved_theme),
        horizontal=True,
        key="theme_setting",
    )
    st.caption("Choose the visual appearance of your application.")
    st.divider()

    st.subheader("Analysis Preferences")
    selected_frame = st.selectbox(
        "Default Translation Reading Frame",
        options=[1, 2, 3],
        index=saved_frame - 1,
        format_func=lambda value: f"Reading Frame {value}",
        key="reading_frame_setting",
    )
    st.caption(
        "Determines the starting nucleotide position used for "
        "the selected-frame protein translation."
    )
    st.divider()

    st.subheader("Report Preferences")
    selected_export = st.selectbox(
        "Default Export Format",
        options=["PDF", "CSV", "TXT"],
        index=["PDF", "CSV", "TXT"].index(saved_export),
        key="export_format_setting",
    )
    st.caption(
        "Your preferred format will appear first in the download section."
    )
    st.divider()

    new_settings = {
        "theme": theme,
        "reading_frame": selected_frame,
        "export_format": selected_export,
    }
    current_preferences = {
        "theme": saved_theme,
        "reading_frame": saved_frame,
        "export_format": saved_export,
    }

    if new_settings != current_preferences:
        update_settings(new_settings)
        st.success("Preferences saved successfully.")
        st.rerun()

    st.divider()
    st.subheader("Data Management")

    history = get_history()
    saved_sequences = get_saved_sequences()
    data_col1, data_col2 = st.columns(2)

    with data_col1:
        st.metric("Analysis Records", len(history))
    with data_col2:
        st.metric("Saved Sequences", len(saved_sequences))

    st.divider()
    st.markdown("#### Clear Analysis History")
    st.caption("Permanently delete all previous analysis records.")

    if st.button("Clear All Analysis History", key="settings_clear_history"):
        st.session_state["confirm_clear_history"] = True

    if st.session_state["confirm_clear_history"]:
        st.warning("This will permanently delete all analysis history.")
        col1, col2 = st.columns(2)
        with col1:
            if st.button(
                "Confirm Clear History",
                key="settings_confirm_clear",
            ):
                clear_history()
                st.session_state["confirm_clear_history"] = False
                st.success("Analysis history deleted.")
                st.rerun()
        with col2:
            if st.button(
                "Cancel",
                key="settings_cancel_clear",
            ):
                st.session_state["confirm_clear_history"] = False
                st.rerun()

    st.divider()
    st.markdown("#### Delete Saved Sequences")
    st.caption("Permanently delete all your bookmarked sequences.")

    if st.button("Delete All Saved Sequences", key="settings_delete_saved"):
        st.session_state["confirm_clear_saved"] = True

    if st.session_state["confirm_clear_saved"]:
        st.warning("This will permanently delete every saved sequence.")
        col1, col2 = st.columns(2)
        with col1:
            if st.button(
                "Confirm Delete All",
                key="settings_confirm_saved",
            ):
                for record in saved_sequences:
                    if record.get("id") is not None:
                        delete_saved_sequence(record["id"])
                st.session_state["confirm_clear_saved"] = False
                st.success("All saved sequences deleted.")
                st.rerun()
        with col2:
            if st.button(
                "Cancel",
                key="settings_cancel_saved",
            ):
                st.session_state["confirm_clear_saved"] = False
                st.rerun()

    st.divider()
    st.caption("Genome Sequence Analysis Tool | Application Settings")


# ==================================================
# ABOUT
# ==================================================

elif active_page == "About":
    st.header("ⓘ About the Project")
    st.write("""
    **Genome Sequence Analysis Tool** is a bioinformatics application
    designed to perform DNA sequence analysis.

    The application uses Python and Biopython for sequence processing
    and biological calculations.

    **Core capabilities**

    - DNA sequence validation and FASTA file parsing
    - Nucleotide composition and ambiguous nucleotide analysis
    - GC and AT content calculation
    - Reverse complement generation and DNA transcription
    - Protein translation in a selected reading frame and all six frames
    - Open reading frame (ORF) identification
    - Codon frequency and usage analysis
    - Sliding-window GC analysis
    - Interactive scientific visualisations
    - PDF, CSV and TXT report generation
    """)
    st.divider()
    st.subheader("Technology Stack")
    st.write(
        "Python · Streamlit · Biopython · Pandas · Plotly · ReportLab"
    )


# ==================================================
# TERMS & CONDITIONS
# ==================================================

elif active_page == "Terms & Conditions":
    st.header("📄 Terms & Conditions")
    st.write("""
    **1. Purpose**

    This application is intended for educational, research and
    exploratory bioinformatics purposes.

    **2. Scientific accuracy**

    Results depend on the quality and correctness of the input DNA sequence.

    **3. Interpretation**

    Generated results should be interpreted with appropriate biological
    knowledge and independently verified where necessary.

    **4. Limitations**

    Protein translation uses the standard genetic code. ORF detection
    identifies ATG start codons and in-frame stop codons. These computational
    results do not establish biological function or gene expression.

    **5. Responsibility**

    Users are responsible for verifying results before using them in
    scientific or clinical work. This application is not a clinical
    diagnostic tool.
    """)


# ==================================================
# MAIN DASHBOARD
# ==================================================

else:
    st.header("Dashboard")
    st.write(
        "Analyse DNA sequences, explore biological properties "
        "and generate downloadable reports."
    )
    st.divider()

    # ----------------------------------------------
    # DNA INPUT
    # ----------------------------------------------

    st.markdown(
        '<p class="section-heading">🧬 DNA Sequence Input</p>',
        unsafe_allow_html=True,
    )

    input_col, info_col = st.columns([3, 1])

    with input_col:
        if "pending_manual_input" in st.session_state:
            st.session_state["manual_sequence_input"] = st.session_state.pop(
                "pending_manual_input"
            )

        sequence_input = st.text_area(
            "Enter or paste your DNA sequence",
            height=150,
            placeholder="Example: ATGGCCATGCGCTAA",
            key="manual_sequence_input",
        )

        if st.button("✓ Validate DNA Sequence"):
            try:
                validated_sequence = validate_dna_sequence(sequence_input)
                st.session_state["dna_sequence"] = validated_sequence
                add_history(
                    sequence=validated_sequence,
                    source="Manual Input",
                )
                st.success(
                    "DNA sequence validated and saved to history!"
                )
                st.rerun()
            except ValueError as error:
                st.error(str(error))

    with info_col:
        st.markdown("#### Input Guidelines")
        st.write("Supported format:")
        st.code("A T G C")
        st.caption(
            "Ambiguous IUPAC nucleotide characters are also supported."
        )

    st.divider()

    # ----------------------------------------------
    # FASTA UPLOAD
    # ----------------------------------------------

    st.markdown(
        '<p class="section-heading">📁 FASTA File Upload</p>',
        unsafe_allow_html=True,
    )

    fasta_col1, fasta_col2 = st.columns([2, 1])

    with fasta_col1:
        uploaded_file = st.file_uploader(
            "Upload your FASTA file",
            type=["fasta", "fa", "fna"],
        )

        if uploaded_file is not None:
            if st.button("📥 Load FASTA File"):
                try:
                    uploaded_file.seek(0)
                    records = parse_fasta(uploaded_file)
                    st.session_state["fasta_records"] = records
                    st.success(
                        f"Successfully loaded {len(records)} sequence(s)!"
                    )
                except (ValueError, OSError, UnicodeError) as error:
                    st.error(str(error))

    with fasta_col2:
        st.markdown("#### FASTA Information")
        st.caption(
            "Upload FASTA files containing one or more DNA sequences."
        )

    if st.session_state["fasta_records"]:
        records = st.session_state["fasta_records"]
        st.markdown("##### Available Sequences")

        selected_index = st.selectbox(
            "Select a sequence",
            range(len(records)),
            format_func=lambda index: (
                f"{records[index]['id']} ({records[index]['length']} bp)"
            ),
            key="fasta_record_selection",
        )

        if st.button("🧬 Analyse Selected FASTA Sequence"):
            selected_record = records[selected_index]
            selected_sequence = selected_record["sequence"]
            st.session_state["dna_sequence"] = selected_sequence
            st.session_state["pending_manual_input"] = selected_sequence
            add_history(
                sequence=selected_sequence,
                source=f"FASTA: {selected_record['id']}",
            )
            st.success(
                f"Sequence '{selected_record['id']}' "
                "loaded and saved to history!"
            )
            st.rerun()

    st.divider()

    # ----------------------------------------------
    # ANALYSIS RESULTS
    # ----------------------------------------------

    sequence = st.session_state.get("dna_sequence", "")

    if sequence:
        st.markdown(
            '<p class="section-heading">📌 Save Sequence</p>',
            unsafe_allow_html=True,
        )

        with st.expander("Save This DNA Sequence"):
            with st.form("save_sequence_form"):
                sequence_name = st.text_input(
                    "Give your sequence a name",
                    placeholder="e.g. Human Gene Sample",
                )
                save_button = st.form_submit_button("Save Sequence")

                if save_button:
                    if not sequence_name.strip():
                        st.error("Please enter a name for this sequence.")
                    else:
                        save_sequence(
                            name=sequence_name.strip(),
                            sequence=sequence,
                        )
                        st.success("Sequence saved successfully!")
                        st.rerun()

        st.divider()

        st.markdown(
            '<p class="section-heading">📊 Sequence Analysis & Metrics</p>',
            unsafe_allow_html=True,
        )

        # Core sequence calculations
        length = calculate_sequence_length(sequence)
        composition = calculate_nucleotide_composition(sequence)
        ambiguous_count = calculate_ambiguous_nucleotide_count(sequence)
        gc_content = calculate_gc_content(sequence)
        at_content = calculate_at_content(sequence)
        reverse_complement = generate_reverse_complement(sequence)
        rna_sequence = transcribe_dna(sequence)
        protein_sequence = translate_dna(
            sequence,
            reading_frame=reading_frame,
        )

        # Extended scientific analyses
        six_frame_translations = translate_six_frames(sequence)
        orfs = find_all_orfs(sequence)
        orf_statistics = calculate_orf_statistics(orfs)
        codon_statistics = calculate_codon_statistics(sequence)

        # Allow sensible window settings for both short and long sequences.
        default_window = min(100, length)
        with st.expander("Sliding-Window GC Settings"):
            window_col1, window_col2 = st.columns(2)
            with window_col1:
                window_size = st.number_input(
                    "Window size (bp)",
                    min_value=1,
                    max_value=max(1, length),
                    value=max(1, default_window),
                    step=1,
                    key="gc_window_size",
                )
            with window_col2:
                step_size = st.number_input(
                    "Step size (bp)",
                    min_value=1,
                    max_value=max(1, length),
                    value=min(10, max(1, length)),
                    step=1,
                    key="gc_step_size",
                )

        gc_window_results = calculate_window_gc(
            sequence,
            window_size=int(window_size),
            step_size=int(step_size),
        )
        gc_statistics = calculate_gc_statistics(gc_window_results)

        # ------------------------------------------
        # METRIC CARDS
        # ------------------------------------------

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Sequence Length", f"{length} bp")
        col2.metric("GC Content", f"{gc_content}%")
        col3.metric("AT Content", f"{at_content}%")
        col4.metric(
            "Ambiguous Bases",
            ambiguous_count,
            help="IUPAC ambiguity symbols other than A, T, G and C.",
        )

        st.caption(
            "GC and AT percentages are calculated using canonical A, T, G "
            "and C bases only. Ambiguous nucleotides are excluded."
        )
        st.divider()

        # ------------------------------------------
        # NUCLEOTIDE COMPOSITION
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">🧪 Nucleotide Composition</p>',
            unsafe_allow_html=True,
        )

        n1, n2, n3, n4 = st.columns(4)
        n1.metric("Adenine (A)", composition.get("A", 0))
        n2.metric("Thymine (T)", composition.get("T", 0))
        n3.metric("Guanine (G)", composition.get("G", 0))
        n4.metric("Cytosine (C)", composition.get("C", 0))
        st.divider()

        # ------------------------------------------
        # SEQUENCE TRANSFORMATIONS
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">🧬 Sequence Results</p>',
            unsafe_allow_html=True,
        )

        with st.expander("Original DNA Sequence", expanded=True):
            st.code(sequence)

        with st.expander("Reverse Complement"):
            st.code(reverse_complement)

        with st.expander("DNA Transcription (RNA)"):
            st.code(rna_sequence)

        with st.expander("Protein Translation — Selected Frame"):
            st.code(protein_sequence)
            st.caption(
                "Translation uses the standard genetic code and "
                f"Reading Frame {reading_frame}."
            )

        with st.expander("Six-Frame Translation"):
            if isinstance(six_frame_translations, dict):
                for frame_name, translated in six_frame_translations.items():
                    st.markdown(f"**{frame_name}**")
                    st.code(str(translated))
            elif isinstance(six_frame_translations, (list, tuple)):
                for index, translated in enumerate(
                    six_frame_translations, start=1
                ):
                    st.markdown(f"**Frame {index}**")
                    st.code(str(translated))
            else:
                st.write(six_frame_translations)

        st.divider()

        # ------------------------------------------
        # ORF ANALYSIS
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">🧬 Open Reading Frame Analysis</p>',
            unsafe_allow_html=True,
        )

        orf_col1, orf_col2, orf_col3 = st.columns(3)
        orf_col1.metric("ORFs Identified", len(orfs))
        orf_col2.metric(
            "Longest ORF (bp)",
            max((orf.get("length", 0) for orf in orfs), default=0),
        )
        orf_col3.metric(
            "Complete ORFs",
            sum(bool(orf.get("is_complete")) for orf in orfs),
        )

        with st.expander("ORF Statistics"):
            st.json(orf_statistics)

        if orfs:
            orf_columns = [
                "strand",
                "reading_frame",
                "start_position",
                "end_position",
                "length",
                "protein_length",
                "start_codon",
                "stop_codon",
                "is_complete",
            ]
            orf_table = pd.DataFrame(orfs)
            visible_columns = [
                column for column in orf_columns
                if column in orf_table.columns
            ]
            st.dataframe(
                orf_table[visible_columns],
                use_container_width=True,
                hide_index=True,
            )

            with st.expander("ORF DNA and Protein Sequences"):
                for index, orf in enumerate(orfs, start=1):
                    st.markdown(
                        f"**ORF {index} — {orf.get('strand')} strand, "
                        f"frame {orf.get('reading_frame')}**"
                    )
                    st.caption(
                        f"Positions {orf.get('start_position')}–"
                        f"{orf.get('end_position')} | "
                        f"{orf.get('length')} bp | "
                        f"{orf.get('protein_length')} amino acids"
                    )
                    st.text("DNA:")
                    st.code(orf.get("sequence", ""))
                    st.text("Protein:")
                    st.code(orf.get("protein_sequence", ""))
        else:
            st.info(
                "No ATG-to-in-frame-stop ORFs were identified in this sequence."
            )

        st.divider()

        # ------------------------------------------
        # CODON ANALYSIS
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">🧪 Codon Analysis</p>',
            unsafe_allow_html=True,
        )

        st.caption(
            "Codon analysis uses non-overlapping triplets from the first "
            "base. Ambiguous codons and incomplete trailing bases are excluded."
        )

        with st.expander("Codon Statistics", expanded=True):
            st.json(codon_statistics)

        st.divider()

        # ------------------------------------------
        # SLIDING-WINDOW GC ANALYSIS
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">📈 Sliding-Window GC Analysis</p>',
            unsafe_allow_html=True,
        )

        st.caption(
            f"Window size: {int(window_size)} bp | "
            f"Step size: {int(step_size)} bp. "
            "Ambiguous bases are excluded from each window's GC denominator."
        )

        with st.expander("GC Window Statistics"):
            st.json(gc_statistics)

        if gc_window_results:
            gc_table = pd.DataFrame(gc_window_results)
            st.dataframe(
                gc_table,
                use_container_width=True,
                hide_index=True,
            )
        else:
            st.info("No GC windows were generated.")

        st.divider()

        # ------------------------------------------
        # VISUALISATIONS
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">📈 Visualisations</p>',
            unsafe_allow_html=True,
        )

        chart_template = (
            "plotly_dark" if current_theme == "Dark" else "plotly_white"
        )
        chart_background = (
            "#111827" if current_theme == "Dark" else "#ffffff"
        )

        chart_col1, chart_col2 = st.columns(2)

        with chart_col1:
            st.subheader("Nucleotide Distribution")
            nucleotide_fig = create_nucleotide_chart(composition)
            nucleotide_fig.update_layout(
                template=chart_template,
                paper_bgcolor=chart_background,
                plot_bgcolor=chart_background,
            )
            st.plotly_chart(
                nucleotide_fig,
                use_container_width=True,
            )

        with chart_col2:
            st.subheader("GC vs AT Distribution")
            if gc_content + at_content > 0:
                content_data = pd.DataFrame({
                    "Component": ["GC Content", "AT Content"],
                    "Percentage": [gc_content, at_content],
                })
                content_fig = px.pie(
                    content_data,
                    names="Component",
                    values="Percentage",
                    title="GC and AT Distribution",
                    hole=0.4,
                    template=chart_template,
                )
                content_fig.update_layout(
                    paper_bgcolor=chart_background,
                    plot_bgcolor=chart_background,
                )
                st.plotly_chart(
                    content_fig,
                    use_container_width=True,
                )
            else:
                st.info(
                    "GC/AT distribution cannot be plotted because the "
                    "sequence contains no canonical A, T, G or C bases."
                )

        st.divider()

        # Codon chart
        st.subheader("Codon Frequency")
        codon_fig = create_codon_frequency_chart(codon_statistics)
        if codon_fig is not None:
            codon_fig.update_layout(
                template=chart_template,
                paper_bgcolor=chart_background,
                plot_bgcolor=chart_background,
            )
            st.plotly_chart(codon_fig, use_container_width=True)
        else:
            st.info("No unambiguous codons are available to visualise.")

        # GC window chart
        st.subheader("GC Content Across Sequence")
        if gc_window_results:
            gc_fig = create_gc_window_chart(gc_window_results)
            if gc_fig is not None:
                gc_fig.update_layout(
                    template=chart_template,
                    paper_bgcolor=chart_background,
                    plot_bgcolor=chart_background,
                )
                st.plotly_chart(gc_fig, use_container_width=True)
        else:
            st.info("No GC-window data is available to visualise.")

        # ORF chart
        st.subheader("ORF Locations")
        if orfs:
            orf_fig = create_orf_visualisation(orfs)
            if orf_fig is not None:
                orf_fig.update_layout(
                    template=chart_template,
                    paper_bgcolor=chart_background,
                    plot_bgcolor=chart_background,
                )
                st.plotly_chart(orf_fig, use_container_width=True)
        else:
            st.info("No ORFs are available to visualise.")

        st.divider()

        # ------------------------------------------
        # DOWNLOAD REPORTS
        # ------------------------------------------

        st.markdown(
            '<p class="section-heading">📥 Download Reports</p>',
            unsafe_allow_html=True,
        )
        st.write(
            "Export your sequence analysis in your preferred format."
        )

        text_report = build_text_report(
            sequence=sequence,
            length=length,
            composition=composition,
            ambiguous_count=ambiguous_count,
            gc_content=gc_content,
            at_content=at_content,
            reverse_complement=reverse_complement,
            rna_sequence=rna_sequence,
            protein_sequence=protein_sequence,
            reading_frame=reading_frame,
            six_frame_translations=six_frame_translations,
            orfs=orfs,
            orf_statistics=orf_statistics,
            codon_statistics=codon_statistics,
            gc_window_results=gc_window_results,
            gc_statistics=gc_statistics,
        )

        pdf_data = generate_pdf_report(
            sequence=sequence,
            length=length,
            composition=composition,
            gc_content=gc_content,
            at_content=at_content,
            reverse_complement=reverse_complement,
            rna_sequence=rna_sequence,
            protein_sequence=protein_sequence,
            ambiguous_count=ambiguous_count,
            six_frame_translations=six_frame_translations,
            orfs=orfs,
            codon_statistics=codon_statistics,
            gc_window_results=gc_window_results,
            gc_statistics=gc_statistics,
        )

        csv_dataframe = build_export_dataframe(
            sequence=sequence,
            length=length,
            composition=composition,
            ambiguous_count=ambiguous_count,
            gc_content=gc_content,
            at_content=at_content,
            reverse_complement=reverse_complement,
            rna_sequence=rna_sequence,
            protein_sequence=protein_sequence,
            reading_frame=reading_frame,
            six_frame_translations=six_frame_translations,
            orfs=orfs,
            codon_statistics=codon_statistics,
            gc_statistics=gc_statistics,
        )
        csv_data = dataframe_download_bytes(csv_dataframe)

        download_options = {
            "PDF": {
                "label": "📕 Download PDF",
                "data": pdf_data,
                "filename": "DNA_Analysis_Report.pdf",
                "mime": "application/pdf",
            },
            "CSV": {
                "label": "📊 Download CSV",
                "data": csv_data,
                "filename": "DNA_Analysis_Results.csv",
                "mime": "text/csv",
            },
            "TXT": {
                "label": "📄 Download TXT",
                "data": text_report,
                "filename": "DNA_Analysis_Report.txt",
                "mime": "text/plain",
            },
        }

        ordered_formats = [export_format] + [
            fmt for fmt in ["PDF", "CSV", "TXT"] if fmt != export_format
        ]
        download_columns = st.columns(3)

        for column, file_format in zip(download_columns, ordered_formats):
            option = download_options[file_format]
            with column:
                button_label = option["label"]
                if file_format == export_format:
                    button_label += " (Default)"

                st.download_button(
                    button_label,
                    data=option["data"],
                    file_name=option["filename"],
                    mime=option["mime"],
                    use_container_width=True,
                    key=f"download_{file_format}",
                )

    else:
        st.info(
            "Enter and validate a DNA sequence or upload a FASTA file "
            "to begin your analysis."
        )


# ==================================================
# FOOTER
# ==================================================

st.divider()
st.markdown(
    "<center style='color: #999999; font-size: 12px;'>"
    "Genome Sequence Analysis Tool | Bioinformatics Workspace"
    "</center>",
    unsafe_allow_html=True,
)
