"""Build the original compensated summation library from portable C++ source."""
from pathlib import Path
import os
import subprocess
import tempfile
from verify_package import verify, safe_path

ROOT = Path(__file__).resolve().parents[1]
def main():
    verify(quiet=True)
    source = safe_path(ROOT, 'code/endpoint/sums.cpp')
    output = source.with_suffix('.so')
    if output.is_symlink():
        raise ValueError('refusing a symlink at the compiled library output')
    with tempfile.TemporaryDirectory(prefix='.sums-build-', dir=source.parent) as temporary:
        built = Path(temporary) / 'sums.so'
        subprocess.run([os.environ.get('CXX', 'g++'), '-std=c++11', '-O2', '-fPIC', '-shared', str(source), '-o', str(built)], check=True)
        os.replace(built, output)
    print(output)

if __name__ == '__main__':
    main()
