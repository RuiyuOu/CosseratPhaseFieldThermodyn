#!/bin/bash
# build_plugin.sh -- builds the solver's plugin kwc2d.so (element-map derivatives by forward-mode AD: kwc2d.cpp, dual.hpp) into the
# repository root, where FreeFem++ loads it, plus a copy in mpi/ (FreeFem++-mpi looks there first), and regenerates cosserat_kwc_mpi.edp
# (and the production v2 solver cosserat_kwc_v2.edp, cosserat_kwc_v2_mpi.edp) when missing or older than their sources.  Rebuilds only what is stale; serialised by a lock.
# Called by run.sh and the submit scripts; safe to run by hand (needs the freefem module: ff-c++).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
exec 8> plugin/.build.lock; flock 8
if [ ! -s kwc2d.so ] || [ plugin/kwc2d.cpp -nt kwc2d.so ] || [ plugin/dual.hpp -nt kwc2d.so ]; then
  command -v ff-c++ >/dev/null || { echo "FAILED: ff-c++ not on PATH (load the freefem module)"; exit 1; }
  # built in a scratch directory: ff-c++ also tries an MPI build into ../mpi when such a directory exists next to the sources
  W=$(mktemp -d); cp plugin/kwc2d.cpp plugin/dual.hpp $W/
  (cd $W && ff-c++ kwc2d.cpp "CXXFLAGS=-fopenmp -O3" -lgomp > build.log 2>&1) || true; cp $W/build.log plugin/build.log
  [ -s $W/kwc2d.so ] || { echo "FAILED: plugin build, see plugin/build.log"; tail -20 plugin/build.log; rm -rf $W; exit 1; }
  cp $W/kwc2d.so kwc2d.so; rm -rf $W; echo "built kwc2d.so"
fi
# the kernels use OpenMP only (no MPI calls), so the same .so serves FreeFem++-mpi
mkdir -p mpi; cmp -s kwc2d.so mpi/kwc2d.so || cp kwc2d.so mpi/kwc2d.so
if [ ! -s cosserat_kwc_mpi.edp ] || [ cosserat_kwc.edp -nt cosserat_kwc_mpi.edp ] || [ mpi/make_mpi.py -nt cosserat_kwc_mpi.edp ]; then
  python3 mpi/make_mpi.py cosserat_kwc.edp cosserat_kwc_mpi.edp
fi
# production v2 solver (v2/README.md): generated from the base solver, then its MPI version
if [ -f v2/make_v2.py ]; then
  if [ ! -s cosserat_kwc_v2.edp ] || [ cosserat_kwc.edp -nt cosserat_kwc_v2.edp ] || [ v2/make_v2.py -nt cosserat_kwc_v2.edp ]; then
    python3 v2/make_v2.py cosserat_kwc.edp cosserat_kwc_v2.edp
  fi
  if [ ! -s cosserat_kwc_v2_mpi.edp ] || [ cosserat_kwc_v2.edp -nt cosserat_kwc_v2_mpi.edp ] || [ mpi/make_mpi.py -nt cosserat_kwc_v2_mpi.edp ]; then
    python3 mpi/make_mpi.py cosserat_kwc_v2.edp cosserat_kwc_v2_mpi.edp
  fi
fi
