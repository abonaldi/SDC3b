# SDC3b results analysis

Code to compare the posterior distributions of the neutral fraction $x_{\rm HI}$ submitted by the SDC3b teams, and to reproduce the results figures and tables of the SDC3b paper.

## Setup

```bash
pip install -r analysis/requirements.txt
python analysis/scripts/fetch_posteriors.py        # downloads all cases (~80 MB)
```

The posteriors are attached to the GitHub release of this repository. They are unpacked into `data/posteriors/<case>/`, which is not tracked by git. Set `SDC3B_DATA` to use another location.

## Data

Each submission is a FITS file `<team>_<case>.fits` holding a $100^3$ grid $P(x_1, x_2, x_3)$, one axis per frequency band:

| Band | Frequency (MHz) | Redshift |
|------|-----------------|----------|
| $x_1$ | 151.0 – 165.9 | 8.46 – 7.57 |
| $x_2$ | 166.0 – 180.9 | 7.57 – 6.92 |
| $x_3$ | 181.0 – 195.9 | 6.92 – 6.31 |

- Cell centres are $x_{\rm HI} = 0.005, 0.015, \dots, 0.995$, as given by the FITS WCS keywords.
- Read with astropy, numpy axis 0 is $x_1$.
- All grids in the release are normalised to unit sum. Submissions that were not normalised have been normalised.

The cases are:
- **PS1, PS2, PS3:** power spectra.
- **IM1:** the image cube whose power spectrum is PS3, so IM1 shares PS3's true values.

## Layout

| Path | Content |
|------|---------|
| `config/teams.yml` | Teams, file names and display labels as used in the paper |
| `config/truth.yml` | True band-centre values and ranges for each case |
| `sdc3b_tools/` | Loading (`io.py`), corner plots and metrics (`plotting.py`) |
| `scripts/fetch_posteriors.py` | Download the release data |
| `notebooks/sdc3b_posteriors_example.ipynb` | Worked example: input histories, corner plots, average posterior, credible intervals, diagnostics |
| `scripts/build_release.py` | Maintainers: build the release zips from the collected submissions |
| `../data/truth/` | Simulation outputs: true values and $x_{\rm HI}(z)$ histories |

## Example

```python
import sys; sys.path.insert(0, 'analysis')
import sdc3b_tools as st

posteriors, coords = st.load_case('PS1')
truth = st.load_truth('PS1')

comp = st.ComparisonPlotter(coords, true_values=truth['centre'], true_ranges=truth['ranges'])
for name, grid in posteriors:
    comp.add_posterior(grid, name)
print(comp.get_metrics_table())   # teams in teams.yml order
comp.plot_boxplot()
```
