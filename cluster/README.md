# FreeFEM + PETSc on the UNC clusters

How to build FreeFEM (FreeFem++, FreeFem++-mpi, ff-c++) with PETSc/SLEPc from source on UNC Research Computing clusters, as a module you can
`module load`, and test it. Written for Sycamore (MPI cluster) and usable on Longleaf (set `MODULES` to Longleaf's compiler/MPI modules).
Nothing here needs root. The build follows FreeFEM's own procedure (configure --enable-download, 3rdparty/getall, ff-petsc `make petsc-slepc`).

| file | what |
|---|---|
| install_freefem_petsc.sh | the build, in restartable stages, each logged to `$BUILD/log_<stage>.txt` |
| install_freefem_petsc.sbatch | runs the heavy stages (petsc, build, install, module) on a compute node |
| test_install.sbatch | checks the installation: plugins, ff-c++ (user plugins), MPI launcher, PETSc solvers on 4 ranks |
| test_plugins.edp, test_petsc.edp | the FreeFEM scripts used by the test |
| write_modulefile.sh | (re)writes the module file of an existing installation (also the last install stage) |

## Quick start (both clusters tested 2026-11-11, `ALL TESTS OK`)

Run from the folder that contains `cluster/` (e.g. the repo root). The script picks the compiler/MPI modules from the host name.

**Sycamore** (gcc 13.3 + openmpi 5.0.6; compute nodes reach the internet; whole compute part ~25 min on 32 cores):
```bash
mkdir -p logs
export BUILD=/hickory/users/<o>/<n>/<onyen>/build/freefem-4.15     # optional; default /work/users/... (or $HOME/build if no /work)
bash cluster/install_freefem_petsc.sh fetch configure getall        # login node, ~3 min
sbatch cluster/install_freefem_petsc.sbatch                         # partition small: petsc, build, install, module file
sbatch cluster/test_install.sbatch                                  # last line of logs/ff_test_<job>.log: ALL TESTS OK
```

**Longleaf** (gcc 12.2 + openmpi 5.0.6-gcc_12.2.0 + autoconf 2.71; partition `general`). The current Longleaf installation was built earlier
with the same recipe by another script and passes `test_install.sbatch`; this script's full build has so far been run on Sycamore only:
```bash
mkdir -p logs
bash cluster/install_freefem_petsc.sh fetch configure getall
sbatch -p general cluster/install_freefem_petsc.sbatch
sbatch -p general cluster/test_install.sbatch
```

**Then, in every shell or job script:** `module use $HOME/modulefiles && module load freefem/4.15`
(installation in `$HOME/sw/freefem/4.15`, about 3 GB; set `PREFIX=...` for another place, `MODDIR=...` for the module folder).
Only the module file lost (installation still there)? `bash cluster/write_modulefile.sh [installation folder]`.

## Module file lost? (the installation itself is elsewhere)

`~/modulefiles/freefem/<v>.lua` only points to the installation. If it was deleted, the installation still works: recreate the file with
`bash cluster/write_modulefile.sh` (it finds `$HOME/sw/freefem/*` or `/hickory|/work/users/<o>/<n>/<onyen>/software/freefem/*`, or give the
folder) and `module use $HOME/modulefiles && module load freefem/4.15`. It also reports whether that installation has PETSc.

## Choices (Sycamore, RHEL 10.2, AMD EPYC 9654/9684X)

* Compiler + MPI: `gcc/13.3.0` + `openmpi/5.0.6/gcc_13.3.0`. GCC 15 is also installed, but its C23 default breaks some older third-party
  sources that PETSc downloads; GCC 13 is the safe choice. Every job must load the same pair: the module file does it.
* CMake: PETSc needs a working cmake to build METIS/ParMETIS/Mmg. On Sycamore `/usr/bin/cmake` exists but does not run under the
  gcc/13.3.0 module (`libstdc++.so.6: version CXXABI_1.3.15 not found`), and PETSc then reports "CMake not found, needed to build METIS".
  The script therefore tests `cmake --version`, then tries the `cmake/4.3.3` module (override: `CMAKE_MODULE=...`), then falls back to a
  pip-installed cmake in `$BUILD/cmake-venv` (self-contained, needs internet). `CMAKE_POLICY_VERSION_MINIMUM=3.5` lets CMake 4 accept the
  older bundled projects. None of this is in the FreeFEM module file.
* Locations: install `PREFIX=$HOME/sw/freefem/4.15` (about 3 GB; home is backed up), build tree in `/work/users/<o>/<n>/<onyen>/build`
  (about 15 GB; can be deleted after the install), module file `$HOME/modulefiles/freefem/4.15.lua`.
* PETSc 3.22.2 (pinned by FreeFEM 4.15, installed in `$PREFIX/ff-petsc/r` and `/c`) brings MUMPS, ScaLAPACK, hypre, METIS/ParMETIS,
  SuperLU(_dist), SuiteSparse, HPDDM, SLEPc and Mmg/ParMmg (adaptive remeshing).
* FreeFEM 4.15 pitfalls handled by the script: the Boost archive needed by bemtool fails its checksum (`--disable-bemtool`, Boost not
  downloaded; boundary elements are not used); CMake 4 rejects old `cmake_minimum_required` lines of bundled packages
  (`CMAKE_POLICY_VERSION_MINIMUM=3.5`).
* Longleaf (rhel9): `MODULES="gcc/12.2.0 openmpi/5.0.6-gcc_12.2.0 autoconf/2.71"` (the script picks it on a Longleaf host), partition `general`:
  `sbatch -p general cluster/install_freefem_petsc.sbatch`, `sbatch -p general cluster/test_install.sbatch`.

## Steps

```bash
# 0. on the login node, in the repo
cd ~/CosseratPhaseFieldThermodyn && mkdir -p logs
srun -p inter -n 1 -t 0:05:00 curl -sI https://web.cels.anl.gov/projects/petsc/download/ | head -1   # do compute nodes reach the internet? (HTTP/2 200)

# 1. login node (needs internet; light work): source, configure, third-party tarballs        ~15 min
bash cluster/install_freefem_petsc.sh fetch configure getall

# 2. compute node: PETSc and its packages, FreeFEM, install, module file                     ~1.5-2 h
sbatch cluster/install_freefem_petsc.sbatch           # log: logs/ff_install_<job>.log (stage logs in the build folder)
#    (if compute nodes have no internet: run step 2 on the login node instead, gently: JOBS=8 bash cluster/install_freefem_petsc.sh petsc build install module)

# 3. test                                                                                     ~5 min
sbatch cluster/test_install.sbatch                     # log: logs/ff_test_<job>.log, last line ALL TESTS OK

# 4. use
module use $HOME/modulefiles && module load freefem/4.15
```

A failed stage prints the last 40 lines of its log; fix the cause and rerun that stage (and the ones after it), e.g.
`bash cluster/install_freefem_petsc.sh configure getall`. To start over: delete the build folder (and `$PREFIX`).

## Using it in job scripts

```bash
module purge
module use $HOME/modulefiles        # or the shared location below
module load freefem/4.15            # loads gcc/13.3.0 + openmpi/5.0.6 as well
mpiexec -n $SLURM_NTASKS --bind-to none FreeFem++-mpi -nw script.edp     # use the mpiexec of the loaded OpenMPI
```
The launcher must be the MPI FreeFem++-mpi was built with; with the wrong one every rank starts as its own 1-rank job (test_install.sbatch
catches this: "only 1 rank").

## Sharing one installation with a group

Install into a group directory and point others to its module folder, e.g.
`PREFIX=/proj/<group>/sw/freefem/4.15 MODDIR=/proj/<group>/modulefiles bash cluster/install_freefem_petsc.sh ...`; users then
`module use /proj/<group>/modulefiles`. (/proj space is given to PIs on request.)

## Other clusters

Set `MODULES` to that cluster's compiler and MPI modules (one compiler, the MPI built with it), e.g. on Longleaf
`MODULES="gcc/<ver> openmpi/<ver>" bash cluster/install_freefem_petsc.sh all`. Check with `module avail gcc openmpi`.

## Log of the Sycamore installation (2026-11)

* Environment survey: RHEL 10.2, 2 x 96-core EPYC per node, Lmod; gcc 13.3/15.2, openmpi 5.0.6/5.0.10, mvapich 4.1, intel 2024/2026, aocc 5.1;
  system cmake 3.31, autoconf 2.71, python 3; login node reaches GitHub and the PETSc server.
* Longleaf, earlier build (build_freefem.sbatch of another session): gcc 12.2 + openmpi 5.0.6, installed in /hickory/.../software/freefem/4.15;
  needed --disable-bemtool and CMAKE_POLICY_VERSION_MINIMUM=3.5 (cmake/4.0.3). Its module file was lost and is recreated by write_modulefile.sh.
  test_install.sbatch on Longleaf (2026-11-11): ALL TESTS OK (plugins, mmg/parmmg, ff-c++, PETSc GAMG/MUMPS/hypre on 4 ranks, max u = 0.0736673).
* Sycamore build (2026-11-11, job on `small`, 32 cores, gcc 13.3 + openmpi 5.0.6, build tree in /hickory/users/r/y/ryo7/build, install in ~/sw/freefem/4.15):
  compute nodes reach the internet; fetch/configure/getall on the login node 3 min; petsc 13 min; build 9 min; install 24 s (the whole job
  about 23 min, far less than the 1.5-2 h estimated). First attempts failed because the compute node has no cmake in PATH, then because
  `/usr/bin/cmake` does not run under the gcc/13.3.0 module (CXXABI_1.3.15); fixed by the cmake_ok check (module, then pip venv cmake 4.4.4).
  `PETSc.so` and `PETSc-complex.so` installed. test_install.sbatch on Sycamore (2026-11-11): ALL TESTS OK (plugins, mmg/parmmg, ff-c++, PETSc GAMG/MUMPS/hypre on 4 ranks, max u = 0.0736673, differences 1.3e-9).
