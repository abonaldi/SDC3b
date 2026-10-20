# SDC3b

Scripts used to generate the simulated data of the third SKA Data Challenge, part b (SDC3b), and to analyse the results submitted by the participating teams.

Results website: <https://abonaldi.github.io/SDC3b/>

## Contents

| Folder | Content |
|--------|---------|
| [`SDC3b_repro/`](SDC3b_repro) | Pipeline that generated the 21-cm signal of the challenge datasets. See its README for the full calling sequence. |
| [`analysis/`](analysis) | Python tools and an example notebook to load the team posteriors, make corner plots and credible-interval plots, and compute the comparison diagnostics used in the paper. |
| [`website/`](website) | Source of the results website (interactive posterior comparison), deployed with GitHub Pages. |
| [`data/truth/`](data/truth) | True neutral fraction of the input simulations: band-centre values, ranges and $x_{\rm HI}(z)$ histories. |

## Team posteriors

The posteriors submitted for PS1, PS2, PS3 and IM1 are attached to the [GitHub releases](https://github.com/abonaldi/SDC3b/releases) of this repository. Download and unpack them with

```bash
pip install -r analysis/requirements.txt
python analysis/scripts/fetch_posteriors.py
```

Then open [`analysis/notebooks/sdc3b_posteriors_example.ipynb`](analysis/notebooks/sdc3b_posteriors_example.ipynb) for a worked example. The file format is described in [`analysis/README.md`](analysis/README.md).

The posteriors are shown for all teams that agreed to the public release of their submissions. The comparison diagnostics describe how each posterior relates to the input values; they are not a ranking of the teams.

## Citation

If you use this repository or the posteriors, please cite the SDC3b paper (reference to be added on publication).
