# SDC3b results website

Static site showing the SDC3b team posteriors: interactive corner plots, credible intervals, comparison metrics and data downloads. It is deployed to GitHub Pages by [`.github/workflows/pages.yml`](../.github/workflows/pages.yml) on every push to `main` that changes this folder.

## Files

| Path | Content |
|------|---------|
| `index.html`, `assets/` | The page (Plotly.js from a CDN, no build step) |
| `data/<case>.js` | Marginals and metrics per dataset, written by `build_data.py` |
| `data/index.js` | List of datasets and the generation date |
| `build_data.py` | Regenerates `data/` from the posterior FITS files |

## Updating the data

The data files are committed, so the workflow only publishes the folder. After the posteriors or `analysis/config/*.yml` change:

```bash
python analysis/scripts/fetch_posteriors.py   # or set SDC3B_DATA to a local copy
python website/build_data.py
git add website/data && git commit -m "Update website data"
```

The metrics are computed with `analysis/sdc3b_tools`, the same code used for the paper.

## Local preview

Open `website/index.html` directly in a browser. The data files are plain scripts, so no server is needed.

Views can be linked directly, e.g. `#PS3/corner`, `#PS1/metrics`, `#IM1/data`.
