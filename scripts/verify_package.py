"""Verify all release payloads and the byte-identical archival source copies."""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

def safe_path(root, relative):
    """Accept regular files inside the release, without symlink traversal."""
    root = Path(root).resolve()
    if not isinstance(relative, str) or not relative or '\\' in relative:
        raise ValueError('unsafe manifest path')
    path = Path(relative)
    if path.is_absolute() or '..' in path.parts or '.' in relative.split('/') or str(path) != relative:
        raise ValueError('unsafe manifest path')
    current = root
    for part in path.parts:
        current /= part
        if current.is_symlink():
            raise ValueError('symlink in release payload: ' + relative)
    if not current.resolve().is_relative_to(root) or not current.is_file():
        raise ValueError('payload missing or outside release: ' + relative)
    return current

def verify(root=ROOT, quiet=False):
    root = Path(root).resolve()
    manifest = json.loads(safe_path(root, 'MANIFEST.json').read_text())
    payloads = {}
    for row in manifest['files']:
        path = safe_path(root, row['path'])
        if row['path'] in payloads or not re.fullmatch('[0-9a-f]{64}', row['sha256']) or type(row['bytes']) is not int or row['bytes'] < 0:
            raise ValueError('invalid or duplicate payload entry')
        raw = path.read_bytes()
        if len(raw) != row['bytes'] or hashlib.sha256(raw).hexdigest() != row['sha256']:
            raise ValueError('payload mismatch: ' + row['path'])
        payloads[row['path']] = row
    for directory in ['code', 'scripts', 'tests']:
        for path in (root / directory).rglob('*.py'):
            relative = str(path.relative_to(root))
            if relative not in payloads:
                raise ValueError('unlisted executable source: ' + relative)
    sources = json.loads(safe_path(root, 'provenance/sources.json').read_text())
    for row in sources['files']:
        path = safe_path(root, row['path'])
        if row['path'] not in payloads or payloads[row['path']]['sha256'] != row['sha256'] or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('archival source changed: ' + row['path'])
    if not quiet:
        print(f"Verified {len(manifest['files'])} release files and {len(sources['files'])} original source copies.")

if __name__ == '__main__':
    verify()
