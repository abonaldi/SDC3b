"""
Assemble the public SDC3b posterior release from the collected submissions.

Reads <source>/<case>/<file>.fits, keeps only public entries in teams.yml,
normalises each grid to unit sum, and writes

    <out>/<case>/<stem>_<case>.fits
    <out>/SDC3b_posteriors_<case>.zip
    <out>/SHA256SUMS

The zips and SHA256SUMS are the files to attach to the GitHub release.

Usage:
    python analysis/scripts/build_release.py ../paper_results dist
"""

import argparse
import hashlib
import sys
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sdc3b_tools.io import CASES, load_teams  # noqa: E402

# Release files whose source file name differs from <stem>_<case>.fits.
# SegU-Net SEarCH was submitted twice under the Modern/Traditional SEarCH names
# (identical files); it is released once under its own name.
SOURCE_OVERRIDES = {
    ('SegU-Net_SEarCH', 'IM1'): 'Modern_SEarCH_IM1.fits',
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def build_case(case, source, out):
    case_dir = out / case
    case_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for team in load_teams(case):
        stem = team['stem']
        src = source / case / SOURCE_OVERRIDES.get((stem, case), f'{stem}_{case}.fits')
        dst = case_dir / f'{stem}_{case}.fits'
        with fits.open(src) as hdul:
            header = hdul[0].header
            grid = np.nan_to_num(hdul[0].data.astype(np.float64), nan=0.0)
        grid /= grid.sum()
        fits.PrimaryHDU(grid, header=header).writeto(dst, overwrite=True)
        written.append(dst)
        print(f'  {case}: {team["name"]:30s} <- {src.name}')

    zip_path = out / f'SDC3b_posteriors_{case}.zip'
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in written:
            zf.write(path, arcname=f'{case}/{path.name}')
    return zip_path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('source', type=Path, help='directory with the collected submissions')
    parser.add_argument('out', type=Path, help='output directory')
    args = parser.parse_args()

    zips = [build_case(case, args.source, args.out) for case in CASES]

    with open(args.out / 'SHA256SUMS', 'w') as f:
        for path in zips:
            f.write(f'{sha256(path)}  {path.name}\n')
    for path in zips:
        print(f'{path.name}: {path.stat().st_size / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
