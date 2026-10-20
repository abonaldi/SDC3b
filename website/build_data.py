"""
Pre-compute the data shown on the SDC3b results website.

For every case it writes website/data/<case>.js with, per team,
    marg1d : 1-D marginals, scaled to unit peak
    marg2d : 2-D credible-level maps, base64 uint8 (see encode_levels)
    stats  : credible intervals and comparison metrics (as in the paper)
and website/data/index.js listing the cases. The data is wrapped in JavaScript
(rather than plain JSON) so the page also works when opened from disk, where
browsers block fetch() of local files.

Requires the posteriors in data/posteriors/ (or $SDC3B_DATA); fetch them with
analysis/scripts/fetch_posteriors.py.

Usage:
    python website/build_data.py
"""

import base64
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from sdc3b_tools import CASES, ComparisonPlotter, load_case, load_truth  # noqa: E402
from sdc3b_tools.io import _read_yaml  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent / 'data'
PAIRS   = [(0, 1), (0, 2), (1, 2)]

# Matches the colour cycle of GridCornerPlotter (tab20 + first 10 of tab20b).
PALETTE = [
    '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
    '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
    '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
    '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5',
    '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939',
    '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31', '#bd9e39',
]

# Credible-level maps store s = -log10(1 - c), where c is the posterior mass
# enclosed by the highest-density region reaching that cell. s is clipped to
# [0, S_MAX] and quantised to one byte; the 68/95/99.7% contours sit at
# s = 0.49, 1.30, 2.52.
S_MAX = 4.0

CASE_INFO = {
    'PS1': 'Foreground-free cylindrical power spectra: EoR signal, 100 h noise and calibration errors.',
    'PS2': 'Foreground-free cylindrical power spectra: EoR signal, 100 h noise and calibration errors.',
    'PS3': 'Foreground-contaminated cylindrical power spectra (1000 h noise), computed from IM1.',
    'IM1': 'Foreground-contaminated image cube (8×8 deg, 151–196 MHz) from which PS3 was computed.',
}

CASE_SIMULATION = {
    'PS1': 'Simulator 1 (pyC2ray)',
    'PS2': 'Simulator 2 (21cmFAST)',
    'PS3': 'Simulator 2 (21cmFAST)',
    'IM1': 'Simulator 2 (21cmFAST)',
}


def team_colours():
    """Colour per team name, fixed by position in teams.yml across all cases."""
    teams = _read_yaml('teams.yml')['teams']
    return {t['name']: PALETTE[i % len(PALETTE)] for i, t in enumerate(teams)}


def encode_levels(m2d):
    """2-D marginal (axes: x_i, x_j) -> base64 uint8 level map, row-major in (x_j, x_i)."""
    flat  = m2d.flatten() / m2d.sum()
    order = np.argsort(flat)[::-1]
    enclosed = np.empty_like(flat)
    enclosed[order] = np.cumsum(flat[order])
    s = -np.log10(np.clip(1.0 - enclosed, 10**-S_MAX, 1.0))
    s[flat <= 0] = S_MAX
    q = np.round(s / S_MAX * 255).astype(np.uint8).reshape(m2d.shape)
    return base64.b64encode(q.T.tobytes()).decode('ascii')


def finite(x, ndigits=4):
    x = float(x)
    return round(x, ndigits) if np.isfinite(x) else None


def team_stats(stats):
    cis = stats['cis']
    return {
        'median':      [finite(c['median']) for c in cis],
        'ci68':        [[finite(c['lo1']), finite(c['hi1'])] for c in cis],
        'ci95':        [[finite(c['lo2']), finite(c['hi2'])] for c in cis],
        'mean':        [finite(v) for v in stats['means']],
        'sigma':       [finite(v) for v in stats['sigmas']],
        'z':           [finite(v, 3) for v in stats['z_scores']],
        'pit':         [finite(v, 3) for v in stats['pit']],
        'mahalanobis': finite(stats['mahalanobis'], 3),
        'score':       finite(stats['sdc3b_score']),
    }


def build_case(case, colours):
    posteriors, coords = load_case(case)
    truth = load_truth(case)
    comp = ComparisonPlotter(coords, true_values=truth['centre'], true_ranges=truth['ranges'])

    teams = []
    for name, grid in posteriors:
        comp.add_posterior(grid, name)
        marg1d = []
        for axis in range(3):
            m = grid.sum(axis=tuple(k for k in range(3) if k != axis))
            marg1d.append([round(float(v), 4) for v in m / m.max()])
        marg2d = {f'{i}_{j}': encode_levels(grid.sum(axis=3 - i - j)) for i, j in PAIRS}
        teams.append({
            'name':   name,
            'colour': colours[name],
            'marg1d': marg1d,
            'marg2d': marg2d,
            'stats':  team_stats(comp.posteriors[-1]['stats']),
        })

    return {
        'case':   case,
        'info':   CASE_INFO[case],
        'simulation': CASE_SIMULATION[case],
        'model':  truth['model'],
        'coords': [round(float(c), 4) for c in coords],
        'truth':  {'centre': truth['centre'], 'ranges': truth['ranges']},
        's_max':  S_MAX,
        'teams':  teams,
    }


def write_js(path, key, obj):
    path.write_text('(window.SDC3B_DATA = window.SDC3B_DATA || {})'
                    f'[{json.dumps(key)}] = {json.dumps(obj, separators=(",", ":"))};\n')


def main():
    OUT_DIR.mkdir(exist_ok=True)
    colours = team_colours()
    index = {'generated': date.today().isoformat(), 'cases': []}
    for case in CASES:
        data = build_case(case, colours)
        path = OUT_DIR / f'{case}.js'
        write_js(path, case, data)
        index['cases'].append({'case': case, 'info': CASE_INFO[case], 'simulation': CASE_SIMULATION[case],
                               'model': data['model'], 'n_teams': len(data['teams'])})
        print(f'{path.name}: {len(data["teams"])} teams, {path.stat().st_size / 1e6:.2f} MB')
    write_js(OUT_DIR / 'index.js', 'index', index)


if __name__ == '__main__':
    main()
