"""Run archival numerical routines in an isolated local output directory."""
import argparse
import importlib.util
import math
import os
from pathlib import Path
import sys
import time
from verify_package import verify

for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(key, '1')
ROOT = Path(__file__).resolve().parents[1]
PACKAGES = {
    'linear-reassembly': 'linear/accepted271',
    'fresh-linear': 'linear/fresh264',
    'nonlinear-baseline': 'nonlinear/baseline308',
    'nonlinear-shifted': 'nonlinear/shifted330',
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=PACKAGES)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wall-seconds', type=float, default=86400)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or args.wall_seconds <= 0:
        parser.error('--wall-seconds must be positive and finite')
    output = args.output.resolve()
    if output.exists():
        parser.error('output already exists; choose a new directory')
    if output.is_relative_to(ROOT / 'code') or output.is_relative_to(ROOT / 'data') or output.is_relative_to(ROOT / 'paper') or output.is_relative_to(ROOT / 'provenance'):
        parser.error('output must not be inside frozen source, evidence, or manuscript directories')
    if output.is_relative_to(ROOT) and not output.is_relative_to(ROOT / 'outputs'):
        parser.error('outputs inside this repository must be under outputs/')
    verify(quiet=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    package = ROOT / 'code' / PACKAGES[args.mode]
    sys.path.insert(0, str(package))
    spec = importlib.util.spec_from_file_location('release_qualification_runner', package / 'qualification_runner.py')
    runner = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = runner
    spec.loader.exec_module(runner)
    deadline = time.monotonic() + args.wall_seconds
    def guard():
        if time.monotonic() >= deadline:
            raise TimeoutError('standalone reproduction wall deadline reached')
    report = runner.run(package, output, guard)
    passed = report.get('whole_qualification_pass', report.get('qualification_pass', False))
    print(f"{args.mode}: {report.get('outcome')}; qualification_pass={passed}; result={output / 'QUALIFICATION_RESULT.json'}")
    print('Numerical reproduction of inspected inputs; no new held-out confirmation or research-harness authority.')
    if not passed:
        raise SystemExit(1)

if __name__ == '__main__':
    main()
