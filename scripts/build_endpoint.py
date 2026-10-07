"""Build the original compensated summation library from portable C++ source."""
from pathlib import Path
import os
import subprocess

ROOT = Path(__file__).resolve().parents[1]
source = ROOT / 'code/endpoint/sums.cpp'
output = source.with_suffix('.so')
subprocess.run([os.environ.get('CXX', 'g++'), '-std=c++11', '-O2', '-fPIC', '-shared', str(source), '-o', str(output)], check=True)
print(output)
