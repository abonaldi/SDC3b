"""
Download the SDC3b team posteriors from the GitHub release and unpack them
into data/posteriors/<case>/ (or $SDC3B_DATA).

Usage:
    python analysis/scripts/fetch_posteriors.py            # all cases
    python analysis/scripts/fetch_posteriors.py PS1 IM1    # selected cases

From a notebook:
    from fetch_posteriors import fetch
    fetch(['PS1'])
"""

import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sdc3b_tools.io import CASES, DATA_DIR  # noqa: E402

RELEASE_TAG = 'v0.1'
RELEASE_URL = f'https://github.com/abonaldi/SDC3b/releases/download/{RELEASE_TAG}'


def _download(name, dest):
    url = f'{RELEASE_URL}/{name}'
    print(f'Downloading {url}')
    urllib.request.urlretrieve(url, dest)


def _checksums(cache):
    path = cache / 'SHA256SUMS'
    if not path.exists():
        _download('SHA256SUMS', path)
    sums = {}
    for line in path.read_text().splitlines():
        digest, name = line.split()
        sums[name] = digest
    return sums


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fetch(cases=None, data_dir=None, force=False):
    """Download and unpack the posteriors for `cases` (default: all)."""
    data_dir = Path(data_dir or DATA_DIR)
    missing = []
    for case in cases or CASES:
        if (data_dir / case).is_dir() and not force:
            print(f'{case}: already present in {data_dir / case}')
        else:
            missing.append(case)
    if not missing:
        return

    cache = data_dir / '_downloads'
    cache.mkdir(parents=True, exist_ok=True)
    sums = _checksums(cache)

    for case in missing:
        name = f'SDC3b_posteriors_{case}.zip'
        zip_path = cache / name
        if not zip_path.exists() or _sha256(zip_path) != sums[name]:
            _download(name, zip_path)
            if _sha256(zip_path) != sums[name]:
                raise RuntimeError(f'Checksum mismatch for {name}')
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(data_dir)
        print(f'{case}: unpacked to {data_dir / case}')


if __name__ == '__main__':
    fetch(sys.argv[1:] or None)
