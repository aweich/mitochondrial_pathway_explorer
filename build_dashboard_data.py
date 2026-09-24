"""Precompute compact tables for the mitochondrial pathway dashboard."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", type=Path, required=True)
    parser.add_argument("--pathways", type=Path,
                        default=Path("../mt_variants/data/mitocarta/MitoPathways3.0.gmx"))
    parser.add_argument("--output", type=Path, default=Path("dashboard_data"))
    parser.add_argument("--cell-type-column", default="cell.type")
    parser.add_argument("--cell-chunk-size", type=int, default=512,
                        help="Number of cells to process at once; lower this if memory is tight")
    return parser.parse_args()


def _zscore_by_row(frame: pd.DataFrame) -> pd.DataFrame:
    means = frame.mean(axis=1)
    stds = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(means, axis=0).div(stds, axis=0).fillna(0.0)


def _aggregate_expression(adata, genes: list[str], cell_types: pd.Index,
                           cell_type_values: np.ndarray,
                           chunk_size: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate selected genes while keeping each matrix slice small."""
    gene_indices = adata.var_names.get_indexer(genes)
    if (gene_indices < 0).any():
        raise ValueError("Some requested genes are missing from AnnData")

    sums = np.zeros((len(cell_types), len(genes)), dtype=np.float64)
    detected = np.zeros_like(sums)
    type_positions = {cell_type: position
                      for position, cell_type in enumerate(cell_types)}

    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        values = adata.X[start:stop, gene_indices]
        labels = cell_type_values[start:stop]
        if hasattr(values, "tocsr"):
            values = values.tocsr()
        else:
            values = np.asarray(values)

        for cell_type, position in type_positions.items():
            mask = labels == cell_type
            if not mask.any():
                continue
            selected = values[mask]
            if hasattr(selected, "getnnz"):
                sums[position] += np.asarray(selected.sum(axis=0)).ravel()
                detected[position] += selected.getnnz(axis=0)
            else:
                selected = np.asarray(selected)
                sums[position] += selected.sum(axis=0)
                detected[position] += np.count_nonzero(selected, axis=0)

    counts = np.asarray([(cell_type_values == cell_type).sum()
                         for cell_type in cell_types], dtype=np.float64)
    mean = pd.DataFrame(sums / counts[:, None], index=cell_types, columns=genes).T
    pct = pd.DataFrame(100 * detected / counts[:, None],
                       index=cell_types, columns=genes).T
    return mean, pct


def build_dashboard_data(adata_path: Path, pathway_path: Path,
                         output_dir: Path,
                         cell_type_column: str = "cell.type",
                         cell_chunk_size: int = 512) -> None:
    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit(
            "anndata is required to build the tables. Install dashboard/requirements.txt."
        ) from exc

    scripts_dir = Path(__file__).resolve().parent.parent / "mt_variants" / "scripts"
    sys.path.insert(0, str(scripts_dir))
    from mitochondrial_pathways import load_mitochondrial_pathways

    if cell_chunk_size < 1:
        raise ValueError("cell_chunk_size must be at least 1")
    adata = ad.read_h5ad(adata_path, backed="r")
    if cell_type_column not in adata.obs:
        raise ValueError(f"Missing required AnnData observation column: {cell_type_column}")
    if not any(str(column).endswith("_UCell") for column in adata.obs.columns):
        raise ValueError(
            "No *_UCell columns found. Provide an AnnData file after upstream "
            "UCell scoring; this pipeline does not score pathways."
        )

    pathways = load_mitochondrial_pathways(pathway_path, available_genes=adata.var_names)
    gene_sets = {path: list(pathways.available_signatures[path])
                 for path in pathways.signature_paths()}
    hierarchy = pathways.hierarchy.set_index("path")
    cell_types = pd.Index(sorted(adata.obs[cell_type_column].dropna().unique()),
                          name=cell_type_column)

    genes = sorted({gene for gene_set in gene_sets.values() for gene in gene_set})
    cell_type_values = adata.obs[cell_type_column].to_numpy()
    gene_mean, gene_pct = _aggregate_expression(
        adata, genes, cell_types, cell_type_values, cell_chunk_size
    )
    gene_zscore = _zscore_by_row(gene_mean)

    ucell_by_celltype = pd.DataFrame(index=cell_types)
    for path in gene_sets:
        column = f"{path}_UCell"
        if column in adata.obs:
            ucell_by_celltype[path] = (adata.obs[column]
                                       .groupby(adata.obs[cell_type_column])
                                       .mean().reindex(cell_types))
        else:
            ucell_by_celltype[path] = np.nan

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_mean.to_parquet(output_dir / "gene_mean_expression.parquet")
    gene_pct.to_parquet(output_dir / "gene_pct_expression.parquet")
    gene_zscore.to_parquet(output_dir / "gene_mean_zscore.parquet")
    ucell_by_celltype.to_parquet(output_dir / "ucell_by_celltype.parquet")

    metadata = hierarchy[["name", "parent", "depth"]].copy()
    metadata = metadata.rename(columns={"name": "display_name"})
    metadata["full_path"] = metadata.index
    metadata["genes"] = metadata.index.map(lambda path: json.dumps(gene_sets[path]))
    metadata["n_genes"] = metadata["full_path"].map(lambda path: len(gene_sets[path]))
    metadata = metadata.reset_index(drop=True)[
        ["full_path", "display_name", "parent", "depth", "n_genes", "genes"]
    ]
    metadata.to_parquet(output_dir / "pathway_metadata.parquet", index=False)
    (output_dir / "gene_sets.json").write_text(json.dumps(gene_sets, indent=2), encoding="utf-8")
    adata.file.close()
    print(f"Wrote {len(genes)} genes and {len(gene_sets)} pathways to {output_dir}")


if __name__ == "__main__":
    args = _parse_args()
    build_dashboard_data(args.adata, args.pathways, args.output,
                         args.cell_type_column, args.cell_chunk_size)
