#!/bin/bash
# one-time: create a python venv with the post-processing dependencies (run on a login node; needs pip access)
# usage: bash setup_py.sh   (optionally `module load python3` first; python >= 3.9)
set -e
python3 -m venv $HOME/dynpost_venv
source $HOME/dynpost_venv/bin/activate
pip install --upgrade pip
pip install numpy scipy matplotlib meshio pillow
python3 -c "import numpy, scipy, matplotlib, meshio, PIL; print('post-processing deps OK')"
which zip || echo "WARNING: 'zip' not found; package.sh needs it (module load zip, or ask admins)"
