"""
Loading of SDC3b configuration, true values and team posteriors.

Posterior files are 3-D FITS grids P(x1, x2, x3) over the neutral fraction in
the three frequency bands. When read with astropy, numpy axis 0 is band z1
(151-166 MHz), axis 1 is z2 and axis 2 is z3. The x_HI value of each cell
centre is taken from the FITS WCS keywords (CRVAL, CDELT, CRPIX); for SDC3b
this is 0.005, 0.015, ..., 0.995 on all axes.
"""

import os
from pathlib import Path

import numpy as np
import yaml
from astropy.io import fits

REPO_DIR   = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_DIR / 'analysis' / 'config'
DATA_DIR   = Path(os.environ.get('SDC3B_DATA', REPO_DIR / 'data' / 'posteriors'))

CASES  = ['PS1', 'PS2', 'PS3', 'IM1']
LABELS = [r'x_{\rm HI,1}', r'x_{\rm HI,2}', r'x_{\rm HI,3}']


def _read_yaml(name):
    with open(CONFIG_DIR / name) as f:
        return yaml.safe_load(f)


def load_teams(case, include_private=False):
    """
    List of team entries (dicts with 'name', 'stem', ...) that submitted `case`.

    Entries with `public: false` are skipped unless include_private=True.
    """
    teams = _read_yaml('teams.yml')['teams']
    return [t for t in teams
            if case in t['cases'] and (include_private or t.get('public', True))]


def load_truth(case):
    """
    True values for `case`.

    Returns a dict with
        centre : list of x_HI at the band centres [x'_1, x'_2, x'_3]
        ranges : list of [min, max] x_HI over each band (the truth cuboid)
        model  : name of the input simulation
    """
    cases = _read_yaml('truth.yml')['cases']
    entry = cases[case]
    if 'same_as' in entry:
        entry = cases[entry['same_as']]
    return {
        'centre': entry['xhi_centre'],
        'ranges': entry['xhi_range'],
        'model':  entry['model'],
    }


def posterior_path(stem, case, data_dir=None):
    return Path(data_dir or DATA_DIR) / case / f'{stem}_{case}.fits'


def coords_from_header(header):
    """Cell-centre coordinates for each numpy axis (axis 0 first)."""
    naxis = header['NAXIS']
    coords = []
    for n in range(naxis, 0, -1):    # FITS axis NAXIS is numpy axis 0
        pix = np.arange(header[f'NAXIS{n}']) + 1
        coords.append(header[f'CRVAL{n}'] + (pix - header[f'CRPIX{n}']) * header[f'CDELT{n}'])
    return coords


def load_posterior(stem, case, data_dir=None):
    """
    Read one posterior grid.

    Returns (grid, coords): grid normalised to unit sum and the shared 1-D
    cell-centre coordinate array. Submissions that were not normalised are
    normalised here.
    """
    with fits.open(posterior_path(stem, case, data_dir)) as hdul:
        grid   = hdul[0].data.astype(np.float64)
        coords = coords_from_header(hdul[0].header)
    for c in coords[1:]:
        if not np.allclose(c, coords[0]):
            raise ValueError(f'{stem}_{case}: axes have different coordinates')
    grid = np.nan_to_num(grid, nan=0.0)
    grid /= grid.sum()
    return grid, coords[0]


def load_case(case, include_private=False, data_dir=None):
    """
    Load all posteriors for `case`.

    Returns (posteriors, coords) where posteriors is a list of
    (team_name, grid) in the order of teams.yml.
    """
    posteriors, coords = [], None
    for team in load_teams(case, include_private):
        grid, c = load_posterior(team['stem'], case, data_dir)
        if coords is None:
            coords = c
        elif not np.allclose(c, coords):
            raise ValueError(f"{team['name']} {case}: grid differs from other teams")
        posteriors.append((team['name'], grid))
    return posteriors, coords
