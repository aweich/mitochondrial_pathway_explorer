"""Interactive mitochondrial pathway dashboard backed only by Parquet/JSON."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st


DATA_DIR = Path(__file__).resolve().parent / "dashboard_data"


@st.cache_data
def load_data(data_dir: str):
    root = Path(data_dir)
    mean = pd.read_parquet(root / "gene_mean_expression.parquet")
    pct = pd.read_parquet(root / "gene_pct_expression.parquet")
    zscore = pd.read_parquet(root / "gene_mean_zscore.parquet")
    ucell = pd.read_parquet(root / "ucell_by_celltype.parquet")
    metadata = pd.read_parquet(root / "pathway_metadata.parquet")
    gene_sets = json.loads((root / "gene_sets.json").read_text(encoding="utf-8"))
    return mean, pct, zscore, ucell, metadata, gene_sets


def _path_label(path: str, metadata: pd.DataFrame) -> str:
    row = metadata.loc[metadata["full_path"] == path].iloc[0]
    return f"{'  ' * int(row['depth'])}{row['display_name']}  [{row['n_genes']} genes]"


def _download_controls(fig, frame: pd.DataFrame, filename: str, key: str) -> None:
    st.download_button("Download CSV", frame.to_csv(index=False).encode("utf-8"),
                       f"{filename}.csv", "text/csv", key=f"csv-{key}")
    st.download_button("Download HTML", fig.to_html(include_plotlyjs="cdn").encode("utf-8"),
                       f"{filename}.html", "text/html", key=f"html-{key}")
    for extension, mime in (("png", "image/png"), ("svg", "image/svg+xml")):
        try:
            image = fig.to_image(format=extension)
        except Exception:
            st.caption(f"{extension.upper()} download requires the Plotly Kaleido package.")
        else:
            st.download_button(f"Download {extension.upper()}", image,
                               f"{filename}.{extension}", mime, key=f"{extension}-{key}")


def _ucell_panel(path: str, ucell: pd.DataFrame) -> None:
    frame = ucell[[path]].rename(columns={path: "ucell_score"}).reset_index(names="cell_type")
    figure = px.imshow(frame.set_index("cell_type").T, aspect="auto",
                        color_continuous_scale="RdBu_r",
                        labels={"x": "Cell type", "y": "Pathway", "color": "UCell score"})
    figure.update_layout(height=350, margin=dict(l=20, r=20, t=40, b=20))
    st.plotly_chart(figure, use_container_width=True)
    _download_controls(figure, frame, f"{path}_ucell_heatmap", "ucell")


def _gene_panel(path: str, genes: list[str], mean: pd.DataFrame,
                pct: pd.DataFrame, zscore: pd.DataFrame) -> None:
    display_genes = genes[:180]
    if len(genes) > len(display_genes):
        st.warning("Showing the first 180 genes for plot readability; all genes remain in the precomputed data.")
    mode = st.radio("Expression color", ["Raw mean", "Z-score"], horizontal=True,
                    key="expression-mode")
    values = mean if mode == "Raw mean" else zscore
    color_name = "mean_expression" if mode == "Raw mean" else "mean_zscore"
    frame = values.loc[display_genes].rename_axis("gene").reset_index().melt(
        id_vars="gene", var_name="cell_type", value_name=color_name)
    fraction = pct.loc[display_genes].rename_axis("gene").reset_index().melt(
        id_vars="gene", var_name="cell_type", value_name="detection_fraction_percent")
    frame = frame.merge(fraction, on=["gene", "cell_type"])
    figure = px.scatter(frame, x="cell_type", y="gene", color=color_name,
                        size="detection_fraction_percent", color_continuous_scale="RdBu_r",
                        hover_data=["gene", "cell_type", color_name,
                                    "detection_fraction_percent"])
    figure.update_layout(height=max(600, len(display_genes) * 22), margin=dict(l=20, r=20, t=40, b=80))
    st.plotly_chart(figure, use_container_width=True)
    _download_controls(figure, frame, f"{path}_dotplot", "dotplot")

    heatmap = px.imshow(values.loc[display_genes], aspect="auto", color_continuous_scale="RdBu_r",
                        labels={"x": "Cell type", "y": "Gene", "color": color_name})
    heatmap.update_layout(height=max(450, len(display_genes) * 18), margin=dict(l=20, r=20, t=40, b=80))
    st.plotly_chart(heatmap, use_container_width=True)
    _download_controls(heatmap, values.loc[display_genes].rename_axis("gene").reset_index(),
                       f"{path}_gene_heatmap", "gene-heatmap")


st.set_page_config(page_title="Mitochondrial pathways", layout="wide")
st.title("Mitochondrial pathway explorer")
st.caption("Pre-computed expression and UCell summaries across CD8 T cell subtypes. Based on MitoCarta 3.0 Pathways.")

try:
    gene_mean, gene_pct, gene_zscore, ucell, metadata, gene_sets = load_data(str(DATA_DIR))
except FileNotFoundError:
    st.error(f"Dashboard data not found at {DATA_DIR}. Run dashboard/build_dashboard_data.py first.")
    st.stop()

path = st.selectbox("Pathway", metadata["full_path"].tolist(),
                    format_func=lambda item: _path_label(item, metadata))
genes = gene_sets[path]
st.subheader(_path_label(path, metadata).strip())
st.write(f"{len(genes)} genes in this pathway; all are shown below.")

tab_ucell, tab_genes = st.tabs(["UCell heatmap", "Gene expression"])
with tab_ucell:
    _ucell_panel(path, ucell)
with tab_genes:
    _gene_panel(path, genes, gene_mean, gene_pct, gene_zscore)
