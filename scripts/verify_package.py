"""Verify all release payloads and the byte-identical archival source copies."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def verify():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    for row in manifest['files']:
        path = ROOT / row['path']
        if not path.is_relative_to(ROOT) or '..' in Path(row['path']).parts:
            raise ValueError('unsafe manifest path')
        raw = path.read_bytes()
        if len(raw) != row['bytes'] or hashlib.sha256(raw).hexdigest() != row['sha256']:
            raise ValueError('payload mismatch: ' + row['path'])
    sources = json.loads((ROOT / 'provenance/sources.json').read_text())
    for row in sources['files']:
        if hashlib.sha256((ROOT / row['path']).read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('archival source changed: ' + row['path'])
    print(f"Verified {len(manifest['files'])} release files and {len(sources['files'])} original source copies.")

if __name__ == '__main__':
    verify()
