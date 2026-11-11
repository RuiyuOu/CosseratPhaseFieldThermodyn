#!/bin/bash
#SBATCH -J ff-test
#SBATCH -p small
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 0:30:00
#SBATCH -o logs/ff_test_%j.log
# Test of the FreeFEM + PETSc installation: versions, plugins, ff-c++ (user plugin), MPI launcher, PETSc solvers (GAMG, MUMPS, hypre) on 4 ranks.
# usage (repo root):  mkdir -p logs; sbatch cluster/test_install.sbatch      env: FF_MODULE (freefem/4.15), MODDIR ($HOME/modulefiles)
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module use ${MODDIR:-$HOME/modulefiles}; module load ${FF_MODULE:-freefem/4.15} || { echo "FAILED: module load ${FF_MODULE:-freefem/4.15}"; exit 1; }
FAIL=0
echo "== executables"; command -v FreeFem++ FreeFem++-mpi ff-mpirun ff-c++ mpiexec
ME=$(grep -o 'mpirun="[^"]*"' "$(command -v ff-mpirun)" | tail -1 | cut -d'"' -f2); echo "ff-mpirun launcher: ${ME:-?}"
echo "== plugins (serial)"; FreeFem++ -nw -v 0 cluster/test_plugins.edp 2>&1 | grep -E "plugins|3D P2|TEST|rror" || true
FreeFem++ -nw -v 0 cluster/test_plugins.edp 2>&1 | grep -q "PLUGIN TEST OK" || { echo "PLUGIN TEST FAILED"; FAIL=1; }
echo "== optional plugins (adaptive meshing: mmg serial, parmmg MPI)"
W=$(mktemp -d); echo 'load "mmg"' > $W/m.edp; echo 'cout << "mmg ok" << endl;' >> $W/m.edp
FreeFem++ -nw -v 0 $W/m.edp 2>&1 | grep -q "mmg ok" && echo "mmg: available" || echo "mmg: not available (optional; needed only for adaptive remeshing)"
rm -rf "$W"
echo "== ff-c++ (user plugin)"
W=$(mktemp -d)
cat > $W/hello.cpp <<'CPP'
#include "ff++.hpp"
double ffhello(const long &n) { return 2.0*n; }
static void Load_Init() { Global.Add("ffhello", "(", new OneOperator1_<double, long>(ffhello)); }
LOADFUNC(Load_Init)
CPP
echo 'load "hello"' > $W/t.edp; echo 'cout << "ffhello(21) = " << ffhello(21) << endl;' >> $W/t.edp
( cd $W && ff-c++ hello.cpp "CXXFLAGS=-fopenmp -O3" -lgomp > build.log 2>&1; FreeFem++ -nw -v 0 t.edp 2>&1 | grep ffhello ) || true
( cd $W && FreeFem++ -nw -v 0 t.edp 2>&1 | grep -q "ffhello(21) = 42" ) && echo "FF-C++ TEST OK" || { echo "FF-C++ TEST FAILED"; tail -20 $W/build.log; FAIL=1; }
rm -rf "$W"
echo "== PETSc on $SLURM_NTASKS ranks"
L=${ME:-mpiexec}
$L -n $SLURM_NTASKS --bind-to none FreeFem++-mpi -nw -v 0 cluster/test_petsc.edp 2>&1 | grep -E "test_petsc|PETSC TEST|rror" | tee $SLURM_SUBMIT_DIR/logs/.ff_test_petsc.txt
grep -q "PETSC TEST OK" $SLURM_SUBMIT_DIR/logs/.ff_test_petsc.txt || FAIL=1; rm -f $SLURM_SUBMIT_DIR/logs/.ff_test_petsc.txt
W=$(mktemp -d); printf 'load "PETSc"\nload "parmmg"\nif (mpirank == 0) cout << "parmmg ok" << endl;\n' > $W/p.edp
$L -n 2 --bind-to none FreeFem++-mpi -nw -v 0 $W/p.edp 2>&1 | grep -q "parmmg ok" && echo "parmmg: available" || echo "parmmg: not available (optional; needed only for parallel adaptive remeshing)"
rm -rf "$W"
[ $FAIL = 0 ] && echo "ALL TESTS OK" || echo "SOME TESTS FAILED (see above)"
