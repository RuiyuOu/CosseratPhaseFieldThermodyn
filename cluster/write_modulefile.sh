#!/bin/bash
# write_modulefile.sh -- (re)write the Lmod module file of an existing FreeFEM installation, e.g. after ~/modulefiles was deleted.
# usage:  bash cluster/write_modulefile.sh [PREFIX]
#   PREFIX: the installation (folder with bin/FreeFem++); default: the first of $HOME/sw/freefem/<v>, /hickory|/work/users/<o>/<n>/<onyen>/software/freefem/<v>
# env:    MODULES (compiler + MPI modules it was built with; default by cluster: Longleaf "gcc/12.2.0 openmpi/5.0.6-gcc_12.2.0",
#         Sycamore "gcc/13.3.0 openmpi/5.0.6/gcc_13.3.0"), MODDIR ($HOME/modulefiles), FF_VERSION (from the folder name)
# The module loads those modules and sets PATH, LD_LIBRARY_PATH, FF_LOADPATH, FF_INCLUDEPATH, CPATH, FREEFEM_ROOT (and PETSC_DIR if PETSc is there).
set -euo pipefail
case "$(hostname -f 2>/dev/null || hostname)" in
  *longleaf*) DEFMOD="gcc/12.2.0 openmpi/5.0.6-gcc_12.2.0" ;;
  *sycamore*) DEFMOD="gcc/13.3.0 openmpi/5.0.6/gcc_13.3.0" ;;
  *) DEFMOD="" ;;
esac
MODULES=${MODULES:-$DEFMOD}
[ -n "$MODULES" ] || { echo "set MODULES (the compiler and MPI modules FreeFEM was built with) on this cluster"; exit 1; }
MODDIR=${MODDIR:-$HOME/modulefiles}
PREFIX=${1:-}
if [ -z "$PREFIX" ]; then
  for c in $HOME/sw/freefem/* /hickory/users/${USER:0:1}/${USER:1:1}/$USER/software/freefem/* /work/users/${USER:0:1}/${USER:1:1}/$USER/software/freefem/*; do
    [ -x "$c/bin/FreeFem++" ] && { PREFIX=$c; break; }
  done
fi
[ -n "$PREFIX" ] && [ -x "$PREFIX/bin/FreeFem++" ] || { echo "no FreeFEM installation found (give its folder: bash cluster/write_modulefile.sh /path/to/freefem/4.15)"; exit 1; }
PREFIX=$(cd "$PREFIX" && pwd)
FF_VERSION=${FF_VERSION:-$(basename "$PREFIX")}
LIBD=$(ls -d "$PREFIX"/lib/ff++/*/lib 2>/dev/null | head -1 || true)
IDPD=$(ls -d "$PREFIX"/lib/ff++/*/idp 2>/dev/null | head -1 || true)
INCD=$(ls -d "$PREFIX"/lib/ff++/*/include 2>/dev/null | head -1 || true)
PETD=$( [ -d "$PREFIX/ff-petsc/r/lib" ] && echo "$PREFIX/ff-petsc/r" || true)
mkdir -p "$MODDIR/freefem"
F=$MODDIR/freefem/$FF_VERSION.lua
{
  echo "-- FreeFEM $FF_VERSION module (written by cluster/write_modulefile.sh on $(hostname -s), $(date +%F))"
  echo "help([[FreeFEM $FF_VERSION in $PREFIX.  Commands: FreeFem++, FreeFem++-mpi, ff-mpirun, ff-c++.  Built with: $MODULES]])"
  echo "whatis(\"Name: FreeFEM\")"
  echo "whatis(\"Version: $FF_VERSION\")"
  for m in $MODULES; do echo "load(\"$m\")"; done
  echo "local root = \"$PREFIX\""
  echo "prepend_path(\"PATH\", pathJoin(root, \"bin\"))"
  echo "prepend_path(\"LD_LIBRARY_PATH\", pathJoin(root, \"lib\"))"
  echo "setenv(\"FREEFEM_ROOT\", root)"
  [ -n "$LIBD" ] && echo "prepend_path(\"FF_LOADPATH\", \"$LIBD\")"
  [ -n "$IDPD" ] && echo "prepend_path(\"FF_INCLUDEPATH\", \"$IDPD\")"
  [ -n "$INCD" ] && echo "prepend_path(\"CPATH\", \"$INCD\")"
  [ -n "$PETD" ] && { echo "prepend_path(\"LD_LIBRARY_PATH\", \"$PETD/lib\")"; echo "setenv(\"PETSC_DIR\", \"$PETD\")"; }
} > "$F"
echo "wrote $F:"; cat "$F"
echo
echo "PETSc in this installation: $( [ -n "$PETD" ] && echo "yes ($PETD)" || echo "NO (MPI runs of the solvers need it: build with cluster/install_freefem_petsc.sh)")"
echo "MPI plugins: $(ls "$LIBD"/mpi/PETSc*.so 2>/dev/null | xargs -n1 basename 2>/dev/null | tr '\n' ' ' || echo none)"
echo "use:  module use $MODDIR && module load freefem/$FF_VERSION"
