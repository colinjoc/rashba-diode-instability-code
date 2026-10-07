"""Bounded integration checks of the released original scientific routines."""
import importlib.util
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

class ReproductionTests(unittest.TestCase):
    def test_primary_at_published_witness(self):
        """A fresh original-kernel evaluation reproduces the saved primary value."""
        sys.path.insert(0, str(ROOT / 'code/endpoint'))
        primary = load('release_primary', 'code/endpoint/primary.py')
        frozen = json.loads((ROOT / 'scientific_data.json').read_text())['perturbation']
        v = primary.decode(frozen['v'])
        result = primary.Primary().evaluate(frozen['x'], frozen['K'], nt=256, nw=512, nested=False)
        matrix = primary.decode(result['table']['256,512']['matrix'])
        coefficient = float(np.vdot(v, matrix @ v).real)
        self.assertAlmostEqual(coefficient, frozen['cP'], delta=2e-13)
        self.assertLess(np.linalg.norm(matrix - matrix.conj().T), 1e-14)

    def test_nonlinear_first_order_matches_independent_linear_response(self):
        """Both independent formulae agree on physical and grazing sample nodes."""
        linear = load('release_linear', 'code/linear/accepted271/linear_response.py')
        nonlinear = load('release_nonlinear', 'code/nonlinear/shifted330/nonlinear_response.py')
        frozen = json.loads((ROOT / 'scientific_data.json').read_text())['perturbation']
        T, D, q, h = frozen['x']
        v = tuple(complex(a, b) for a, b in zip(frozen['v']['real'], frozen['v']['imag']))
        eta = np.array([(v[0] - 1j * v[1]).conjugate(), 0, v[0] + 1j * v[1]])
        for n in [0, 32]:
            for helicity in [-1, 1]:
                for direction in [0., .47, 1.]:
                    with self.subTest(n=n, helicity=helicity, direction=direction):
                        z = complex((2*n + 1)*np.pi*T, (q/2 - helicity*h)*direction)
                        projection = frozen['K'][0]*direction
                        O, a, rO, ra = nonlinear.background(z, D, 0, 0, 1e-12)
                        options = dict(rO=rO, ra=ra, rD=0, ruK=0, floor=1e-12, guard=lambda: None)
                        U, eU = nonlinear.sector(eta, a, O, D, projection, 2, **options)
                        V, eV = nonlinear.sector(np.conj(eta[::-1]), a, O, D, -projection, 2, **options)
                        _, _, coeff, identity, _ = nonlinear.force_series(U, V, a, eta, 2, eU=eU, eV=eV, ra=ra, floor=1e-12, guard=lambda: None)
                        _, H = linear.force_matrix(z, D, projection)
                        expected = linear.contract(H, v)
                        self.assertAlmostEqual(coeff[0], expected, delta=1e-12)
                        self.assertLess(max(nonlinear.norm(x) for x in identity[:3]), 1e-12)

    def test_accepted_linear_reassembly_matches_saved_evidence(self):
        with tempfile.TemporaryDirectory(prefix='rashba-linear-test-') as tmp:
            output = Path(tmp) / 'result'
            subprocess.run([sys.executable, str(ROOT / 'scripts/reproduce.py'), 'linear-reassembly', '--output', str(output)], check=True, capture_output=True, text=True)
            report = json.loads((output / 'QUALIFICATION_RESULT.json').read_text())
        saved = json.loads((ROOT / 'data/evidence/LINEAR271.json').read_text())
        for key in ['independent_coefficient', 'complete_physical_EI', 'Ectrl']:
            self.assertAlmostEqual(report[key], saved[key], delta=1e-16)
        self.assertEqual(report['checks'], saved['checks'])

if __name__ == '__main__':
    unittest.main()
