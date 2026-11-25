import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def create_nucleotide_chart(composition: dict):
    """Create an interactive nucleotide distribution chart."""

    df = pd.DataFrame({
        "Nucleotide": list(composition.keys()),
        "Count": list(composition.values())
    })

    fig = px.bar(
        df,
        x="Nucleotide",
        y="Count",
        title="Nucleotide Composition",
        labels={
            "Nucleotide": "DNA Nucleotide",
            "Count": "Number of Bases"
        },
        text="Count"
    )

    fig.update_layout(
        template="plotly_white",
        title_x=0.5
    )

    return fig


def create_codon_frequency_chart(codon_statistics: dict):
    """Create a chart showing codon frequencies."""

    counts = codon_statistics.get("codon_counts", {})

    if not counts:
        fig = go.Figure()
        fig.add_annotation(
            text="No valid codons available for analysis.",
            showarrow=False
        )
        fig.update_layout(
            title="Codon Frequency",
            template="plotly_white"
        )
        return fig

    df = pd.DataFrame({
        "Codon": list(counts.keys()),
        "Frequency": list(counts.values())
    })

    df = df.sort_values(
        by="Frequency",
        ascending=False
    )

    fig = px.bar(
        df,
        x="Codon",
        y="Frequency",
        title="Codon Frequency Distribution",
        text="Frequency",
        labels={
            "Codon": "DNA Codon",
            "Frequency": "Occurrences"
        }
    )

    fig.update_layout(
        template="plotly_white",
        title_x=0.5
    )

    return fig


def create_gc_window_chart(window_results: list):
    """Visualise GC content across sliding DNA windows."""

    if not window_results:
        fig = go.Figure()
        fig.add_annotation(
            text="No GC window data available.",
            showarrow=False
        )
        fig.update_layout(
            title="GC Content Distribution",
            template="plotly_white"
        )
        return fig

    df = pd.DataFrame(window_results)

    fig = px.line(
        df,
        x="start_position",
        y="gc_content",
        markers=True,
        title="GC Content Across DNA Sequence",
        labels={
            "start_position": "Sequence Position",
            "gc_content": "GC Content (%)"
        },
        hover_data=[
            "end_position",
            "valid_bases",
            "ambiguous_bases"
        ]
    )

    fig.update_layout(
        template="plotly_white",
        title_x=0.5,
        yaxis_range=[0, 100]
    )

    return fig


def create_orf_visualisation(orfs: list):
    """Visualise ORF locations on forward and reverse strands."""

    if not orfs:
        fig = go.Figure()
        fig.add_annotation(
            text="No complete ORFs detected.",
            showarrow=False
        )
        fig.update_layout(
            title="Open Reading Frame Distribution",
            template="plotly_white"
        )
        return fig

    df = pd.DataFrame(orfs)

    fig = go.Figure()

    for strand in ("+", "-"):

        strand_data = df[df["strand"] == strand]

        if strand_data.empty:
            continue

        fig.add_trace(
            go.Bar(
                x=(
                    strand_data["end_position"]
                    - strand_data["start_position"]
                    + 1
                ),
                base=strand_data["start_position"],
                y=strand_data["frame"],
                orientation="h",
                name=f"{strand} Strand",
                text=[
                    f"ORF {i + 1}"
                    for i in range(len(strand_data))
                ],
                hovertemplate=(
                    "Frame: %{y}<br>"
                    "Start: %{base}<br>"
                    "Length: %{x} bp"
                    "<extra></extra>"
                )
            )
        )

    fig.update_layout(
        title="Open Reading Frame Distribution",
        xaxis_title="DNA Sequence Position",
        yaxis_title="Reading Frame",
        template="plotly_white",
        title_x=0.5,
        barmode="overlay"
    )

    return fig