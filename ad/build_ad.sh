#!/bin/bash
# build_ad.sh -- AD variant of the 2D solver: builds ad/kwc2d.so (copied to the repo root, where FreeFem loads it) and generates
# cosserat_kwc_ad.edp (make_ad.py) and cosserat_kwc_ad_mpi.edp (mpi/make_mpi.py) when missing or older than their sources.
# Called by run.sh when TANGENT=ad; safe to run by hand (needs the freefem module: ff-c++).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
exec 8> ad/.build.lock; flock 8
if [ ! -s kwc2d.so ] || [ ad/kwc2d.cpp -nt kwc2d.so ] || [ ad/dual.hpp -nt kwc2d.so ]; then
  command -v ff-c++ >/dev/null || { echo "FAILED: ff-c++ not on PATH (load the freefem module)"; exit 1; }
  # built in a scratch directory: ff-c++ also tries an MPI build into ../mpi when such a directory exists next to the sources
  W=$(mktemp -d); cp ad/kwc2d.cpp ad/dual.hpp $W/
  (cd $W && ff-c++ kwc2d.cpp "CXXFLAGS=-fopenmp -O3" -lgomp > build.log 2>&1) || true; cp $W/build.log ad/build.log
  [ -s $W/kwc2d.so ] || { echo "FAILED: see ad/build.log"; tail -20 ad/build.log; rm -rf $W; exit 1; }
  cp $W/kwc2d.so kwc2d.so; rm -rf $W; echo "built kwc2d.so"
fi
# FreeFem++-mpi looks for plugins in ./mpi/ first (not ./): keep a copy there (the kernels use OpenMP only, no MPI calls, so the same .so works)
mkdir -p mpi; cmp -s kwc2d.so mpi/kwc2d.so || cp kwc2d.so mpi/kwc2d.so
if [ ! -s cosserat_kwc_ad.edp ] || [ cosserat_kwc.edp -nt cosserat_kwc_ad.edp ] || [ ad/make_ad.py -nt cosserat_kwc_ad.edp ]; then
  python3 ad/make_ad.py cosserat_kwc.edp cosserat_kwc_ad.edp
fi
if [ ! -s cosserat_kwc_ad_mpi.edp ] || [ cosserat_kwc_ad.edp -nt cosserat_kwc_ad_mpi.edp ] || [ mpi/make_mpi.py -nt cosserat_kwc_ad_mpi.edp ]; then
  python3 mpi/make_mpi.py cosserat_kwc_ad.edp cosserat_kwc_ad_mpi.edp
fi
