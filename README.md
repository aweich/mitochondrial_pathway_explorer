# Mitochondrial pathway dashboard

Streamlit app for exploring precomputed CD8 T-cell mitochondrial pathway
expression and UCell scores. Runtime data is in `dashboard_data/`.

## Run locally

```bash
python -m pip install -r requirements.txt
streamlit run app.py
```

## Rebuild data

The input `.h5ad` must contain `adata.obs["cell.type"]` and existing
`*_UCell` columns. The builder does not calculate UCell scores.

```bash
python build_dashboard_data.py \
  --adata ../mt_variants/data/CD8T/CD8_scored_with_signatures.h5ad \
  --pathways ../mt_variants/data/mitocarta/MitoPathways3.0.gmx \
  --output dashboard_data \
  --cell-chunk-size 128
```

For Streamlit Community Cloud, select `app.py` as the app entrypoint
and commit the generated `dashboard_data/` files.
