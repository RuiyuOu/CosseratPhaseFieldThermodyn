#!/bin/bash
#SBATCH -J a1mesh
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --array=0-4
#SBATCH -o logs/a1_%A_%a.log
# A1 (Longleaf, 2D repository): KWC mesh study.  The production set-up (same 8 grains, same physics, ANGLE / UPMS) on near-uniform meshes
# of element size h = 0.2, 0.14, 0.1, 0.05, 0.035 um (array task 0-4; the production mesh.msh is ~0.06 um and is the run you already have),
# zero-load relaxation (NRELAX steps) then the plate, up to TEND ns.  Then grain_post (boundary lengths per misorientation bin over time)
# and remesh/b1_curvature.py (rotation front, refinement indicator) on each run.  Compare with remesh/a1_compare.py.
# usage (repository root, after 'bash plugin/build_plugin.sh'):  sbatch remesh/a1.sbatch     env: ANGLE (90), UPMS (400), TEND (5), NRELAX (100)
# output: remesh/a1/h<h>/ (mesh, solver copy, out/ = run, grains/, b1/)
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module use "$HOME/modulefiles"; module load freefem/4.15 sqlite/3.45.0
HS=(0.2 0.14 0.1 0.05 0.035); ID=${SLURM_ARRAY_TASK_ID:-0}; H=${HS[$ID]}
export ANGLE=${ANGLE:-90}; UPMS=${UPMS:-400}; TEND=${TEND:-5}; MODE=mpi; source ./tag.sh || exit 1
PHYS=$(echo "$PHYS" | sed -E "s/-nrelax [^ ]+/-nrelax ${NRELAX:-100}/")
W=remesh/a1/h$H; mkdir -p $W/out
module load python3 2>/dev/null || true; [ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"   # numpy etc. (setup_py.sh)
# mesh: the production recipe (mesh/make_mesh.py) with segment length h on grain boundaries and outer sides (buildmesh fills the inside alike)
if [ ! -s $W/mesh.msh ]; then
  python3 mesh/make_mesh.py mesh/MPTD_2D.geo mesh/MPTD_2D.ori --hgb $H --hout $H --out $W/meshbuild > /dev/null || exit 1
  (cd $W/meshbuild && FreeFem++ -nw -v 0 make_mesh.edp > /dev/null) && cp $W/meshbuild/mesh.msh $W/ || { echo "FAILED: mesh"; exit 1; }
fi
# the solver reads mesh.msh, grain_angles_rad.txt and kwc2d.so from its own folder (FreeFem++ -cd: a solver option, also FreeFem's chdir)
cp cosserat_kwc_mpi.edp grain_angles_rad.txt grain_angles_deg.txt kwc2d.so $W/; mkdir -p $W/mpi; cp kwc2d.so $W/mpi/
echo "A1 h=$H: $(grep -c . $W/mesh.msh) lines in mesh.msh; tag $TAG, $SLURM_NTASKS ranks"
FFM=$(command -v ff-mpirun); ME=$(grep -o 'mpirun="[^"]*"' "$FFM" | tail -1 | cut -d'"' -f2); [ -x "$ME" ] || ME=mpiexec
O=$W/out/a1
if ! tail -1 ${O}_log.csv 2>/dev/null | awk -F, -v e=$TEND '{exit !($2+0 >= e - 0.011)}'; then
  $ME -n $SLURM_NTASKS --bind-to none FreeFem++-mpi -nw -v 0 $W/cosserat_kwc_mpi.edp -pistonup $UN -pistonut $UT -ptr 0.5 -pth 2.0 -tend $TEND -dt 0.01 \
     -vtkevery 10 -ckevery 50 -tstop 1400 $PHYS -npexpect $SLURM_NTASKS -out $(realpath $W/out)/a1 > $W/out/a1.out 2>&1
  echo "solver exit $?; last log line: $(tail -1 ${O}_log.csv)"
fi
# post-processing
(cd post && BASE=a1 NP=2 python3 grain_post.py ../$W/out ../$W/grains ../$W/anim) || echo "grain_post failed"
python3 remesh/b1_curvature.py $W/out --hmin 0.1 --hmax 0.4 --out $W/b1 || echo "b1 failed"
