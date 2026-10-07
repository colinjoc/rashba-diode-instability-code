"""Regression checks for release integrity and filesystem confinement."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from verify_package import safe_path, verify

class IntegrityTests(unittest.TestCase):
    def test_path_traversal_and_absolute_paths_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            for relative in ['../outside', '/etc/passwd', 'a/../../outside', './data.json', 'a\\b']:
                with self.subTest(relative=relative), self.assertRaises(ValueError):
                    safe_path(Path(temporary), relative)

    def test_symlink_payload_and_parent_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'real').mkdir()
            (root / 'real/data').write_text('example')
            (root / 'link').symlink_to(root / 'real', target_is_directory=True)
            (root / 'linked-file').symlink_to(root / 'real/data')
            for relative in ['link/data', 'linked-file']:
                with self.subTest(relative=relative), self.assertRaises(ValueError):
                    safe_path(root, relative)

    def test_tampering_and_unlisted_python_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'provenance').mkdir()
            (root / 'scripts').mkdir()
            (root / 'payload').write_bytes(b'original')
            (root / 'provenance/sources.json').write_text('{"files": []}')
            files = []
            for name in ['payload', 'provenance/sources.json']:
                raw = (root / name).read_bytes()
                files.append(dict(path=name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
            (root / 'MANIFEST.json').write_text(json.dumps(dict(files=files)))
            verify(root, quiet=True)
            (root / 'payload').write_bytes(b'tampered')
            with self.assertRaisesRegex(ValueError, 'payload mismatch'):
                verify(root, quiet=True)
            (root / 'payload').write_bytes(b'original')
            (root / 'scripts/unlisted.py').write_text('raise RuntimeError("unexpected")')
            with self.assertRaisesRegex(ValueError, 'unlisted executable'):
                verify(root, quiet=True)

    def test_provenance_paths_are_also_confined(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'provenance').mkdir()
            raw = json.dumps(dict(files=[dict(path='/etc/passwd', sha256='0'*64)])).encode()
            (root / 'provenance/sources.json').write_bytes(raw)
            (root / 'MANIFEST.json').write_text(json.dumps(dict(files=[dict(path='provenance/sources.json', bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())])))
            with self.assertRaisesRegex(ValueError, 'unsafe manifest path'):
                verify(root, quiet=True)

if __name__ == '__main__':
    unittest.main()
