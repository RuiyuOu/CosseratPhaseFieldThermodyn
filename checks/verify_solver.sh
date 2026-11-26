#!/bin/bash
# verify_solver.sh -- the frozen solver files (version 8.0) against SOLVER_SHA256; exit 0 if all match.
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
sha256sum -c SOLVER_SHA256 && echo "solver frozen: all files match SOLVER_SHA256"
