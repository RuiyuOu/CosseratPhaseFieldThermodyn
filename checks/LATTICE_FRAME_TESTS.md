# Lattice-frame tests: stress and slip follow the lattice

`checks/test_stress_follows_lattice.py` and `checks/test_slip_follows_lattice.py` check that rotating the elasticity tensor ℂ and the slip systems once is correct. They are rotated by the grain angle θ₀ at setup and never updated, yet the stress and the slip still follow the current lattice orientation R_t = R R₀. This holds however the lattice turned: rigid rotation, slip, the grain-boundary spin w*, or any mix.

Both tests are pointwise, at a single material point, and need only Python ≥ 3.8 and numpy; matplotlib is optional (figures). They run in a few seconds and do not need FreeFem, the mesh or a checkpoint. The proof is in `checks/lattice_frame_proof.pdf` (source `lattice_frame_proof.md`).

```
python3 checks/test_stress_follows_lattice.py            # both elastic laws, figures in checks/out/
python3 checks/test_slip_follows_lattice.py
python3 checks/test_slip_follows_lattice.py --law prod --no-fig --tol 1e-10
```

Exit code 0 means PASS, 1 means FAIL. Run both after any change to the kinematics (F = R Ū F_p), the elastic law, the slip-system set-up, the rotation of ℂ, or the F_p update.

## What is computed

The model at one point is built in `checks/lattice_frame_common.py`:

- **Kinematics.** F = R Ū F_p with the exact constraint R = polar(F F_p⁻¹). The multiplier λ is set by the torque balance.
- **Stresses.** P = R (Π + 2λ ∂a/∂Ū) F_p⁻ᵀ and τ = P Fᵀ. The Mandel stress is M = ŪᵀΠ, and the resolved shear stress is τ_α = M:(s_α ⊗ n_α).
- **Derivatives.** Taken by complex step, so they are exact to round-off.
- **Elastic laws.**
  - `svk`: plain St Venant–Kirchhoff, ½E:ℂ:E with E = ½(ŪᵀŪ − I), plus 2μ_c a².
  - `prod`: the production law, Vinet W(J) + ½Ẽ:ℂ:Ẽ (isochoric Green strain) + 2μ_c a².
- **Constants.** Cu (C11, C12, C44 = 168.4, 121.4, 75.4 GPa) in the solver's lattice frame x = [110], y = [001], z = [1̄10]; μ_c = 754 GPa.
- **Slip systems.** The solver's three: s = (1, 0), (cos a_R, ± sin a_R), a_R = atan √2, n = (−s_y, s_x).
- **Comparison with the solver.** At start-up both tests read `cosserat_kwc_mpi.edp` (or `cosserat_kwc.edp`) from the repository root. They compare `Clxx`, `Clyy`, `Clxy`, `Clss`, `Clxz`, `Clyz` and the slip-direction definition with the values used here. A mismatch is a FAIL; a source that cannot be found is a WARN.

**Load set.** 35 stretch states V, each at 12 lattice angles α between −100° and 170°:

- uniaxial x and y, and equibiaxial, at 0.02, 0.1, 0.25 and −0.2;
- three pure shears;
- 20 random large stretches, J ≈ 0.75–1.3.

The orientation changes the stress by up to 116 % (prod) and 138 % (svk), so agreement is not trivial.

**Reference.** A crystal that starts at orientation α: ℂ and the slip systems rotated by α, F = V, F_p = I.

**Test cases.** ℂ and the slip systems are kept at orientation 0, and the lattice is turned by α during the history.

| Case | History | Expected |
|---|---|---|
| T1 | rigid rotation: F = V R(α), F_p = I | = reference |
| T2 | grain-boundary spin: F = V, F_p = R(−α) (material not rotated) | = reference |
| T3 | slips on the three solver systems (γ = 0.5, −0.2, 0.3), then R(−α); F = V R(α) F_p | = reference |
| T4 | 200 steps of F_p ← exp(Δw* Ŵ) F_p (exact, as in `trialFp`) to w* = −73°, while F goes from I to V | = reference at every 20th step |
| W | T2 with ℂ (and, in the slip test, s_α, n_α) also rotated by α | must differ (> 0.1 / > 0.05), and must equal a crystal at 2α |
| Ws | slip test only: T2 with only s_α, n_α rotated | must differ (> 0.05) |

**What is compared.**

- The stress test compares τ = P Fᵀ: the maximum over the components of |τ − τ_ref| / max|τ_ref|.
- The slip test compares:
  - τ_α for the three systems, relative to max|τ_ref|;
  - the current slip directions F_e s_α and normals F_e⁻ᵀ n_α;
  - the push-forward identity F_e L_p F_e⁻¹ = Σ γ̇_α (F_e s_α) ⊗ (F_e⁻ᵀ n_α).
- Tolerance: `--tol`, default 10⁻¹⁰.

## Expected output (17–18 Nov 2026)

Stress, maximum relative error:

| law | T1 | T2 | T3 | T4 | W | W vs 2α |
|---|---|---|---|---|---|---|
| svk | 1.4e-14 | 1.9e-14 | 2.2e-14 | 2.4e-13 | 1.70 | 2.3e-14 |
| prod | 2.8e-14 | 1.6e-14 | 3.1e-14 | 2.7e-13 | 1.16 | 1.4e-14 |

Slip, maximum error in τ_α (relative) / current vectors (absolute):

| law | T1 | T2 | T3 | T4 | W | Ws | W vs 2α | push-forward |
|---|---|---|---|---|---|---|---|---|
| svk | 3.3e-15 / 4.4e-16 | 3.5e-15 / 6.7e-16 | 5.8e-15 / 8.9e-16 | 4.4e-14 / 4.2e-15 | 2.5 | 2.1 | 4.8e-15 | 8.9e-16 |
| prod | 2.9e-15 / 4.4e-16 | 3.7e-15 / 6.7e-16 | 6.6e-15 / 8.9e-16 | 3.7e-14 / 4.2e-15 | 2.7 | 2.3 | 5.7e-15 | 8.9e-16 |

**Figures** (`checks/out/`):

- `stress_follows_lattice_<law>.png`: τ₁₁, τ₂₂, τ₁₂ against α for V = diag(1.10, 0.95); lines are the reference, markers the tests. Includes the T4 time history, the W case against the 2α curve, and a summary panel.
- `slip_follows_lattice_<law>.png`: the same for τ_α of the three systems, plus summary panels for τ_α and the current slip vectors.

## Why W and Ws must fail

ℂ, s_α and n_α are fixed in the isoclinic intermediate configuration, and the lattice rotation is carried by F_e (and, through the constraint, by R). Rotating ℂ or the slip systems as well when the lattice turns applies the rotation twice: W reproduces a crystal at 2α exactly. These two cases show that the tests can detect a wrong implementation.

## Scope

- **Pointwise only.** The tests check the constitutive kinematics at a point. They do not cover finite elements, curvature (KWC) torques or the balance laws.
- **Constraint enforcement.** The constraint is exact here (a = 0). In the solver it is exact in the element mean (P0 multiplier), and the 2μ_c a² penalty controls the rest.
- **Coverage elsewhere.** The AD element map, the tangent and the energy derivatives are covered by `checks/run_tangent_check.sh`.
