# Tschopp Cu tilt-boundary energies -> KWC fit

1. Download (login node; compute nodes usually have no internet), standard library only:

       python3 tschopp_grab.py --out tschopp_data

   Writes `tschopp_data/gb_energies.csv` (file, metal, axis, sigma, angle_deg, energy) and `tschopp_data/info_heads.txt`
   (first lines of every info file). If the CSV is empty or the columns look wrong, send `info_heads.txt`.

2. Fit (needs numpy + scipy; matplotlib optional). One process per boundary:

       sbatch -J kwcfit -p general -N 1 -n 1 -c 8 --mem=8G -t 02:00:00 -o kwcfit_%j.log --wrap 'source ~/dynpost_venv/bin/activate; python3 kwc_fit_tschopp.py tschopp_data/gb_energies.csv --metal Cu --axis 110 --max-angle 30 --delta 1e-3 --np 8 --out fit_cu110'

   If scipy is missing from the venv: `source ~/dynpost_venv/bin/activate; pip install scipy`.
   The pre-2026-11-10 production energy had eps = 0.004 um, h = 1 and delta = 0.005 (run scripts in the repository history).
   To fit in that setting with eps held fixed (fits f0 and s only, much faster), add `--eps 0.004 --delta 0.005 --h 1`.
   Without `--eps`, eps is fitted too (eps and s are nearly interchangeable for boundary energies alone).

Output: `fit_cu110/fit_summary.txt` (eps = 0 fit, eps fit, production values, gamma at 2-45 deg) and `fit_cu110/kwc_fit.png`.
Model: gamma = f0 min int[(1-phi)^2/2 + a^2 phi'^2/2 + s g(phi)(sqrt(theta'^2+delta^2)-delta) + eps^2/2 h(phi) theta'^2] dx, h = 1 (default; `--h phi2` for phi^2),
g = -2[phi + ln(1-phi)]; a is held at 0.10 um (it sets the width; gamma depends on f0 a, s/a, eps/a).
