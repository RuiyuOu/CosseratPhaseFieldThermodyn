#!/bin/bash
# install_freefem_petsc.sh -- build FreeFEM (with PETSc/SLEPc, MUMPS, hypre, HPDDM, Mmg/ParMmg) from source on a UNC cluster and make it a module.
#
# usage:  bash cluster/install_freefem_petsc.sh <stage> [<stage> ...]
#   stages, in order (each one is restartable and writes $BUILD/log_<stage>.txt):
#     fetch      clone FreeFEM-sources at v$FF_VERSION, autoreconf             (login node: needs internet)
#     configure  ./configure --enable-download --enable-optim --prefix=$PREFIX
#     getall     ./3rdparty/getall -a: download the third-party tarballs       (login node: needs internet)
#     petsc      3rdparty/ff-petsc: make petsc-slepc (PETSc downloads MUMPS, hypre, ScaLAPACK, METIS, ParMETIS, Mmg, ParMmg, ...), ./reconfigure
#     build      make -j$JOBS
#     install    make install
#     module     write the Lmod module file $MODDIR/freefem/$FF_VERSION.lua
#     all        every stage above
# env:  FF_VERSION (4.15)  PREFIX ($HOME/sw/freefem/$FF_VERSION)  BUILD (/work/users/<o>/<n>/<onyen>/build/freefem-$FF_VERSION, else $HOME/build/...)
#       MODULES (Sycamore "gcc/13.3.0 openmpi/5.0.6/gcc_13.3.0", Longleaf "gcc/12.2.0 openmpi/5.0.6-gcc_12.2.0 autoconf/2.71")
#       JOBS (number of cores; 8 on a login node)  MODDIR ($HOME/modulefiles)
# Known FreeFEM 4.15 pitfalls handled here: the Boost archive of bemtool fails its checksum (-> --disable-bemtool, not downloaded);
# CMake 4 rejects old cmake_minimum_required lines of bundled packages (-> CMAKE_POLICY_VERSION_MINIMUM=3.5).
#
# Heavy stages (petsc, build) take 1-2 h: run them on a compute node (cluster/install_freefem_petsc.sbatch, or an interactive session).
# The compiler and MPI chosen here must be the ones the jobs load: the module file loads them.
set -uo pipefail
FF_VERSION=${FF_VERSION:-4.15}
PREFIX=${PREFIX:-$HOME/sw/freefem/$FF_VERSION}
WORKDIR=/work/users/${USER:0:1}/${USER:1:1}/$USER
[ -d "$WORKDIR" ] || WORKDIR=$HOME
BUILD=${BUILD:-$WORKDIR/build/freefem-$FF_VERSION}
case "$(hostname -f 2>/dev/null || hostname)" in
  *longleaf*) DEFMOD="gcc/12.2.0 openmpi/5.0.6-gcc_12.2.0 autoconf/2.71" ;;
  *) DEFMOD="gcc/13.3.0 openmpi/5.0.6/gcc_13.3.0" ;;
esac
MODULES=${MODULES:-$DEFMOD}
MODDIR=${MODDIR:-$HOME/modulefiles}
if [ -z "${JOBS:-}" ]; then if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then JOBS=$SLURM_CPUS_PER_TASK; else JOBS=8; fi; fi
SRC=$BUILD/FreeFem-sources
[ $# -ge 1 ] || { sed -n '2,20p' "$0"; exit 1; }

module purge >/dev/null 2>&1
module load $MODULES || { echo "FAILED: module load $MODULES"; exit 1; }
export CMAKE_POLICY_VERSION_MINIMUM=3.5
# PETSc's METIS/ParMETIS/Mmg downloads need a WORKING cmake. On Sycamore /usr/bin/cmake exists but dies under the gcc/13.3.0 module
# ("libstdc++.so.6: version CXXABI_1.3.15 not found": it was built with a newer GCC and the module's libstdc++ shadows the system one).
# So test it, then try a cmake module, then a pip-installed (self-contained) cmake in a venv inside the build folder.
# CMake 4 is fine because of CMAKE_POLICY_VERSION_MINIMUM above. None of this goes into the FreeFEM module file.
cmake_ok() { cmake --version >/dev/null 2>&1; }
cmake_ok || module load ${CMAKE_MODULE:-cmake/4.3.3} >/dev/null 2>&1
if ! cmake_ok; then
  CV=$BUILD/cmake-venv
  [ -x "$CV/bin/cmake" ] || { mkdir -p "$BUILD" && python3 -m venv "$CV" && "$CV/bin/pip" install -q cmake; } >&2
  export PATH="$CV/bin:$PATH"
fi
cmake_ok || { echo "FAILED: no working cmake (system, module and pip all failed)"; exit 1; }
echo "  cmake: $(command -v cmake)  $(cmake --version | head -1)"
for t in cmake gcc g++ gfortran mpicc mpicxx mpif90 autoreconf automake libtool make patch m4 bison flex git wget; do command -v $t >/dev/null || { echo "FAILED: $t not found after 'module load $MODULES'"; exit 1; }; done
export CC=gcc CXX=g++ FC=gfortran F77=gfortran MPICC=mpicc MPICXX=mpicxx MPIFC=mpif90
mkdir -p "$BUILD"
echo "FreeFEM $FF_VERSION -> $PREFIX   build: $BUILD   modules: $MODULES   jobs: $JOBS"
echo "  $(gcc --version | head -1);  $(mpicc --version 2>/dev/null | head -1);  mpiexec: $(command -v mpiexec)"

run() {   # run <stage> <command...>: log, and stop on failure with the tail of the log
  local st=$1; shift; local lg=$BUILD/log_$st.txt
  echo "[$(date +%T)] $st ... (log: $lg)"
  if ! ( "$@" ) > "$lg" 2>&1; then echo "FAILED in stage $st; last lines of $lg:"; tail -n 40 "$lg"; exit 1; fi
  echo "[$(date +%T)] $st done"
}
st_fetch() {
  if [ ! -d "$SRC/.git" ]; then git clone --branch v$FF_VERSION --depth 1 https://github.com/FreeFem/FreeFem-sources.git "$SRC"; fi
  cd "$SRC" && autoreconf -i
}
CONF="--enable-download --enable-optim --disable-bemtool --prefix=$PREFIX"
# every getall package except bemtool/Boost (checksum failure in 4.15)
GETALL=${GETALL:-ARPACK,BLAS,OpenBLAS,FFTW,freeYams,Gmm++,Ipopt,METIS,ParMETIS,MMG3D,mshmet,MUMPS,NLopt,ScaLAPACK,Scotch,SuiteSparse,SuperLU,TetGen,mmg,parmmg,PETSc,htool,hpddm,libpthread-google}
st_configure() { cd "$SRC" && ./configure $CONF; }
st_getall() { cd "$SRC" && ./3rdparty/getall -a -o "$GETALL" && ./configure $CONF; }   # configure again so the downloads are picked up
st_petsc() { cd "$SRC/3rdparty/ff-petsc" && make petsc-slepc && cd "$SRC" && ./reconfigure; }
st_build() { cd "$SRC" && make -j"$JOBS"; }
st_install() { cd "$SRC" && make install; }
st_module() { MODULES="$MODULES" MODDIR="$MODDIR" FF_VERSION="$FF_VERSION" bash "$(dirname "${BASH_SOURCE[0]}")/write_modulefile.sh" "$PREFIX"; }
for stage in "$@"; do
  case $stage in
    all) for s in fetch configure getall petsc build install module; do run $s st_$s; done ;;
    fetch|configure|getall|petsc|build|install) run $stage st_$stage ;;
    module) st_module ;;
    *) echo "unknown stage $stage"; exit 1 ;;
  esac
done
case " $* " in *" module "*|*" all "*) echo; echo "use it:  module use $MODDIR; module load freefem/$FF_VERSION      then: sbatch cluster/test_install.sbatch";; esac
