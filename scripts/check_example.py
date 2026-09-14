#!/usr/bin/env python3
"""Check the committed example's files; optionally decode its media with FFmpeg."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess


def check(folder, *, decode=False):
    folder = folder.resolve()
    manifest = json.loads((folder / 'manifest.json').read_text())
    expected = manifest['files']
    actual = {
        str(p.relative_to(folder))
        for p in folder.rglob('*')
        if p.is_file() and p.name not in {'manifest.json', '.DS_Store'}
    }
    if actual != set(expected):
        raise ValueError(f'File list differs: missing={set(expected)-actual}, extra={actual-set(expected)}')
    for relative, digest in expected.items():
        path = folder / relative
        if not path.resolve().is_relative_to(folder):
            raise ValueError('Example path leaves its directory: ' + relative)
        with path.open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != digest:
                raise ValueError('Checksum mismatch: ' + relative)
    media = [folder / p for p in expected if Path(p).suffix in {'.png', '.gif', '.mp4'}]
    counts = dict(Counter(p.suffix[1:] for p in media))
    if counts != manifest['media_counts']:
        raise ValueError('Media counts differ from the manifest')
    if decode:
        for path in media:
            subprocess.run(
                ['ffmpeg', '-v', 'error', '-xerror', '-i', str(path), '-f', 'null', '-'],
                check=True, stdout=subprocess.DEVNULL,
            )
    return dict(files=len(expected), media=counts, checksums='PASS',
                decoding='PASS' if decode else 'NOT_CHECKED')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'examples/central_fans')
    parser.add_argument('--decode', action='store_true', help='Decode all images and movies with FFmpeg')
    args = parser.parse_args()
    print(json.dumps(check(args.folder, decode=args.decode), indent=2))


if __name__ == '__main__':
    main()
