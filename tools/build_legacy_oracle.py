"""Build old kernels in an isolated directory and produce FEM oracle fixtures.

This is an explicit development tool. Normal builds never read legacy sources.
It compiles the original paper-case implementation and headers, with only the
two mesh constructors extracted from model.cpp to avoid the old LOD dependency.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def extract_function(text, signature):
    start = text.index(signature)
    end = text.index('{', start) + 1
    depth = 1
    while depth:
        if text[end] == '{':
            depth += 1
        elif text[end] == '}':
            depth -= 1
        end += 1
    return text[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--work', type=Path, default=ROOT/'_local/oracle')
    parser.add_argument('--fixtures', type=Path, default=ROOT/'tests/fixtures')
    args = parser.parse_args()
    src = args.source.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src/'include', work/'include', dirs_exist_ok=True)
    for header in ['problems.hpp', 'meshes.hpp']:
        (work/'include/alod').mkdir(exist_ok=True)
        shutil.copyfile(ROOT/'include/alod'/header, work/'include/alod'/header)
    units = ['src/mesh/refine.cpp', 'src/mesh/edges.cpp', 'src/helmholtz/boundary.cpp',
             'src/helmholtz/quadrature.cpp', 'src/helmholtz/operators.cpp',
             'src/helmholtz/manufactured.cpp', 'src/helmholtz/benchmarks/paper_cases.cpp']
    hashes = {}
    for unit in units:
        shutil.copyfile(src/unit, work/Path(unit).name)
        hashes[unit] = hashlib.sha256((src/unit).read_bytes()).hexdigest()
    model = (src/'src/helmholtz/model.cpp').read_text(encoding='utf-8-sig')
    constructors = '\n'.join(extract_function(model, 'TriMesh '+name+'()') for name in
                            ['make_helmholtz_unit_square_mesh', 'make_helmholtz_l_shape_mesh'])
    (work/'initial.cpp').write_text('#include "helmholtz/model.h"\n#include "helmholtz/boundary.h"\n'
        '#include <cmath>\nnamespace lod2d::helmholtz {\n'+constructors+'\n}\n', encoding='utf-8')
    hashes['src/helmholtz/model.cpp'] = hashlib.sha256((src/'src/helmholtz/model.cpp').read_bytes()).hexdigest()
    for path in ['apps/fem_smoke.cpp', 'src/problems/problems.cpp']:
        shutil.copyfile(ROOT/path, work/Path(path).name)
    cmake = '''cmake_minimum_required(VERSION 3.20)
project(LegacyKernelOracle LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_path(UMFPACK_INCLUDE_DIR umfpack.h PATH_SUFFIXES suitesparse REQUIRED)
find_library(UMFPACK_LIBRARY NAMES umfpack REQUIRED)
add_executable(legacy_oracle fem_smoke.cpp problems.cpp initial.cpp SOURCES)
target_include_directories(legacy_oracle PRIVATE include ${UMFPACK_INCLUDE_DIR})
target_compile_definitions(legacy_oracle PRIVATE LOD2D_HAVE_UMFPACK=1)
target_link_libraries(legacy_oracle PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX ${UMFPACK_LIBRARY})
'''.replace('SOURCES', ' '.join(Path(x).name for x in units))
    (work/'CMakeLists.txt').write_text(cmake, encoding='utf-8')
    subprocess.run(['cmake', '-S', str(work), '-B', str(work/'build'), '-DCMAKE_BUILD_TYPE=Release'], check=True)
    subprocess.run(['cmake', '--build', str(work/'build'), '-j', '2'], check=True)
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
    args.fixtures.mkdir(parents=True, exist_ok=True)
    for problem in ['E1', 'E2']:
        result = subprocess.run([str(work/'build/legacy_oracle'), problem, '6'],
                                check=True, capture_output=True, text=True, env=env, timeout=90)
        fixture = {'provenance': 'Original E2 production-snapshot kernels; nominal FEM only, not an adaptive replay.',
                   'source_sha256': hashes, 'output': json.loads(result.stdout)}
        (args.fixtures/f'fem_{problem}.json').write_text(json.dumps(fixture, indent=2)+'\n', encoding='utf-8')
        print(f'Wrote {problem} legacy oracle')


if __name__ == '__main__':
    main()
