"""Write BUILD_MANIFEST.json for the frozen GMRES build.

Records source and binary SHA-256 hashes, toolchain, flags, Python and
package versions, BLAS/LAPACK link, CPU, kernel and CPU frequency governor.

usage: python tests/write_build_manifest.py
"""
import hashlib
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

CODE = Path(__file__).resolve().parents[1]
SRC = CODE / 'msspec_gmres' / 'src'
FORTRAN = SRC / 'msspec' / 'spec' / 'fortran' / 'phd_gm_noso_nosp_nosym'
OUT = CODE / 'BUILD_MANIFEST.json'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout.strip()


def rel(p):
    return str(p.relative_to(CODE))


def source_hashes():
    files = sorted(list(SRC.rglob('*.f')) + list(SRC.rglob('*.py'))
                   + [SRC / 'options.mk'] + [CODE / 'check_release.sh']
                   + sorted((CODE / 'tests').glob('*.py')))
    files = [f for f in files if '_venv' not in f.parts and 'build' not in f.parts]
    entries = {rel(f): sha256(f) for f in files if f.is_file()}
    combined = hashlib.sha256()
    for name in sorted(entries):
        combined.update((name + ' ' + entries[name] + '\n').encode())
    return entries, combined.hexdigest()


def binary_hashes():
    bins = sorted(SRC.rglob('*.so'))
    return {rel(b): sha256(b) for b in bins}


def blas_link():
    so = SRC / 'msspec' / 'spec' / 'fortran' / '_phd_gm_noso_nosp_nosym.so'
    lines = run(['ldd', str(so)]).splitlines() if so.exists() else []
    keep = [l.strip() for l in lines if re.search(r'lapack|blas|openblas|gfortran', l)]
    return keep


def cpu_info():
    model = run(['lscpu']).split('Model name:')[-1].splitlines()[0].strip() if 'Model name:' in run(['lscpu']) else ''
    gov_path = Path('/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor')
    gov = gov_path.read_text().strip() if gov_path.exists() else 'unknown'
    return {'model': model, 'logical_cpus': run(['nproc']), 'cpu_governor': gov}


def main():
    src_entries, src_combined = source_hashes()
    try:
        import numpy, scipy
        numpy_v, scipy_v = numpy.__version__, scipy.__version__
    except ImportError:
        numpy_v = scipy_v = 'not installed'
    try:
        import ase
        ase_v = ase.__version__
    except ImportError:
        ase_v = 'not installed'
    makefile = (SRC / 'options.mk').read_text()
    fflags = re.search(r'GFORTRAN_FFLAGS\s*=\s*(.*)', makefile).group(1).strip()
    f2py_flags = re.search(r'F2PYFLAGS\s*=\s*(.*)', makefile).group(1).strip()
    manifest = {
        'generated_utc': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'source_tree_sha256': src_combined,
        'source_files_sha256': src_entries,
        'binary_files_sha256': binary_hashes(),
        'toolchain': {
            'gfortran': run(['gfortran', '--version']).splitlines()[0],
            'f2py': run(['f2py', '-v']),
            'gfortran_flags': fflags,
            'f2py_flags': f2py_flags,
            'build_command': 'cd msspec_gmres/src/msspec/spec/fortran && make all BUILDDIR=$PWD/build VERBOSE=1',
        },
        'python': {
            'version': platform.python_version(),
            'executable_note': 'code/_venv',
            'numpy': numpy_v,
            'scipy': scipy_v,
            'ase': ase_v,
        },
        'blas_lapack_link': blas_link(),
        'platform': {
            'kernel': platform.release(),
            'os': run(['uname', '-srm']),
            **cpu_info(),
        },
        'notes': [
            'Timing numbers in the material check were taken under the powersave governor.',
            'The material check report is in tests/material_check_report.json.',
        ],
    }
    OUT.write_text(json.dumps(manifest, indent=2) + '\n')
    print('wrote', OUT)
    print('source_tree_sha256', src_combined)
    print('binaries', len(manifest['binary_files_sha256']))


if __name__ == '__main__':
    sys.exit(main())
