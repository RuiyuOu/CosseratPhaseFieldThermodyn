# `-baek 1`: duality-based Cosserat–KWC coupling (experimental code only)

File: `cosserat2d_mob_baek.edp` (2026-11-03): the stable code + `-baek 1`, kept separate from the KWC crystal-symmetry code
(`cosserat2d_mob_exp.edp` is back to its pre-Baek state; the development version with both is `cosserat2d_mob_exp_baekdev.edp`). Serial only so far.
Defaults of the existing code are unchanged; everything below is active only with `-baek 1`.

## Model
* Elastic law θ-free (μ_c = 0): ψ(C^e).
  * The lattice rotation θ_e = atan((F^e_21 − F^e_12)/(F^e_11 + F^e_22)) comes from the polar decomposition of F^e = F F^p⁻¹.
    It is an output of the element map, with dθ_e/dF^e = R^e Ŵ / tr U^e.
* Orientation field ϑ ∈ HCT (C¹), an absolute angle. It carries the KWC curvature energy
  W(∇ϑ) = f0 (s g(φ)(√(|κ|²+δ²) − δ) + ε²|κ|²/2).
* Discrete Lagrangian per step:
  L = Σ_e |e| ψ_e + ∫[W(κ) − m·κ + m·∇θ0s] − Σ_e θ_e ∫_{∂e∖∂Ω} m·n + s·2μ_c ∫(ϑ − θ0s − θ_e)² − (h²·proj·Ā/8) Σ_e ‖Δϑ − ⟨Δϑ⟩_e‖²
  * m = W′(κ).
  * The edge sum is the distributional ⟨m, ∇θ_e⟩ for element-wise constant θ_e. This is the IBP form of Baek's Eq. 6.28 (stage 1).
* Element-map input is the torque z_e = (1/|e|)(∫_{∂e∖∂Ω} m·n + s·4μ_c ∫_e(ϑ − θ0s)), i.e. the element mean of div m = 2 × skew stress.
  * ψ_e gets −z θ_e + s·2μ_c θ_e².
  * It reaches the slip and GB-spin drives through the Mandel stress F^eᵀ ∂ψ/∂F^e.

## Robustness choices: each one was needed on the bicrystal
1. **`-baekls −1` (concave least squares).** The dual part −W*(m) + m·∇θ_L is maximised by ϑ = θ_L, so the ϑ block is negative definite.
   Baek's +1 gives −div(A∇) + 4μ_c. That is indefinite and near-singular wherever A/ℓ² ≈ 4μ_c, which the degenerate KWC A
   (tiny in the GB cores) always produces. With +1 the ϑ initialisation stalled (30 Newton its, |r| only 50× down); with −1 it took 7.
   Default `-baekc 1e-4` (μ_c = 4.8 MPa); Baek reports results insensitive for c = 1e-5 … 1e-2.
2. **`-baeklag 1` (lagged diffusivity, default).** Within a step, W is replaced by its quadratic majorant ½A_n|κ|² with
   A_n = W′(|κ_n|)/|κ_n| from the step start.
   * W(√t) is concave in t, so this is the classical MM/lagged-diffusivity scheme for TV-like energies: energy-decreasing, first order in time.
   * Without it, the exact dual Lagrangian is concave only within ~δ of ∇ϑ = ∇θ_L. Step 2 of the relaxation then stalled:
     line-search steps 1/64, |R| 19 → 3.4 after 12 its. Gauss-Newton (`-baeklag 0 -baekgn 1`) was slow.
   * With lagging, the ϑ part of every step is linear. S_z, the ϑ Hessian and the constant residual part are assembled once per step,
     and the per-iteration HCT assembly disappears.
3. **`-baekspin 1` (GB spin dw* as a global unknown, default).**
   * The element equation is r_w = μ|e|(k_spin τ* − dw*), solved monolithically with (u, θ, ϑ).
   * With the θ-free law nothing restrains the spin locally: at fixed z, k_spin τ* reached ~7 rad in the GB cores.
     That hit the spin cap and 2π branches, local residuals stayed at 6.3, and step 3 stalled.
   * The restraint is the curvature response of the whole ϑ system, which only the monolithic solve sees.
     The local problem is now slip only.
4. **Uncapped flow law.** `flow()` caps x = |τ|/(f τ_c) at 300 as an overflow guard. In Baek mode it now uses the overflow-free
   form gd = x/(cc + x^{1−m}/g0), with no cap and an exact derivative. The Baek torque drives ~40 elements beyond x = 300 in the
   first relaxation steps; the cap made the tangent wrong (1.8e-4) and the residual non-smooth.

## Newton system
Unknowns are (u, θ, ϑ, dw*). θ (P1 Cosserat) is a dummy that stays 0. Blocks:

```
[ K_uu            B_l S_z                      B_w            ]
[ −S_zᵀ G         K_ϑϑ − S_zᵀ C S_z            −S_zᵀ g_w      ]
[ −μ k_s T_k      −μ k_s t_z S_z               μ(1 − k_s t_w) ]
```

* K_uu is the existing consistent element tangent.
* B_l, B_w are d(stress)/dz and d(stress)/d dw*.
* G, C, g_w are dθ_e/dk, dθ_e/dz and dθ_e/d dw*.
* T_k, t_z, t_w are the derivatives of τ*.
* All include the implicit slip response. Solved with UMFPACK.

ϑ is initialised by an L2 projection of θ0s, then equilibrated with θ_e = 0 (lagged: 6 outer updates of A_n).
The checkpoint side file is `<ck>.thb`. Quadrature is `tripleQF(qf5pT)` on the HCT sub-triangles: on the stage-1 manufactured
solution it is as accurate as qf7pT and half the cost.
The KWC phase field, energy and output use the element means of ∇ϑ. VTK fields:
* `dtheta_deg` = θ_e
* `theta_kwc_deg` = ϑ (folded)
* `m_x, m_y` = element means of m
* `w_cosserat` = θ_e (rad)

## Options
* `-baek 1`
* `-baekc` (1e-4)
* `-baekproj` (4)
* `-baekls` (−1)
* `-baeklag` (1)
* `-baekspin` (1)
* `-baekgn` (0)
* `-baekinit` (30)
* `-tceps` (FD step of `-tancheck`, 1e-7)

Requires `-variant A -tangent consistent -newton full -linsolver umfpack` and `-symfold 0 -kwcfield 0 -mucdual 0 -kwcreg 0`.

## Checks (bicrystal, serial)
* Tangent vs central FD (`-tancheck 1 -tceps 1e-5`): all blocks 1e-6 … 1e-8. The exception is u-rows along d(dw*), at 7e-2 of
  |FD| = 1.6e-4, i.e. FD round-off. The ϑ-direction errors fall with the FD step (1e-4 → 1e-7 for ε = 1e-7 → 1e-4): round-off.
* KWC off (`-kwcs 0 -kwce 0`), 2 relaxation and 3 piston steps: `-baek 1` reproduces `-muc 0`. All logged quantities agree to
  < 1e-6 relative, including the initial residual of every step. This is the stage-2 checkpoint (run before items 2–4).
* KWC on, relaxation: every step converges with full Newton steps, quadratically at the end (6–8 iterations).
* Cost ≈ 70 s per relaxation step serial on the 9600-element bicrystal, against ≈ 4 s for the plain μ_c = 0 step.
  Most of it is the UMFPACK solve of the 83k-dof block system (≈ 6 s per iteration); the element map and tangent take ≈ 2 s.

## Finding: KWC relaxation with exact compatibility (bicrystal, 20 relaxation steps)
* Default (η*₀ = 1e-3): the GB-core elements (φ ≈ 0.025, ~610 elements along x = 5) spin without bound.
  * θ_e reached 79° after 19 steps; dw* is ~1.5 rad per step; Newton needs 6 → 35 iterations; it stalled at step 20 (`BGS2`).
  * The penalty run on the same schedule (`PEN20`, μ_c = 75.4) stays at |Δθ| ≤ 0.7° and |w*| ≤ 0.012.
  * The rotated elements are exactly the GB core (all 610 with φ < 0.5).
* Cause: the dual is the μ_c → ∞ limit, so ϑ = θ_lattice is enforced.
  * The smoothed initial orientation θ0s keeps a KWC gradient in the grain interiors next to the GB. Its couple-stress flux
    m ≈ f0 s g(φ) sign(κ) enters the core elements as a torque z.
  * Only the interior lattice could remove that gradient, and it cannot rotate: k_spin ≈ 0 where φ ≈ 1, and elastic rotation is
    resisted. The core can spin freely: k_spin = 5 at φ → 0, and the θ-free law gives no elastic resistance. The core torque
    therefore never relaxes.
  * With the penalty, the Cosserat θ in the interior deviates from the lattice by a ≈ m/(4μ_c h). That flattens the KWC gradient,
    and the skew stress carries the rest.
* `-baekls +1 -baekc 1e-2` (0.48 GPa least squares) does not help: it stalls at step 5.
* `-etastar 0.1` (100× slower GB spin) relaxes smoothly over 8 steps (5–6 Newton its each), with |θ_e| ≤ 0.43° and |w*| ≤ 0.004.
  The core torque is still there, though, so the rotation only grows more slowly.
* Candidate remedies (open):
  * (a) start from the sharp grain orientation θ0 (P0, no smoothing), which the edge form handles, so the interior carries no
    KWC gradient;
  * (b) a GB-spin law that is not driven by the full KWC torque in the core;
  * (c) keep a finite compatibility (penalty) in the relaxation.

## Update 2026-11-03: sharp initial orientation (option 1)
* `-baeksharp 1` (default): θ0 = raw grain angle (P0), entering the dual edge sum with θ_e; `-baeksharp 0` = smoothed th0s.
* Least-squares reference is now the equilibrated initial ϑ_0 (not θ0), so the LS term only fixes the constant mode and is 0 at t = 0.
  KWC off: reproduces `-muc 0` again (all logged quantities < 1e-7 relative).
* Bicrystal relaxation (`SH2`): the GB-core spin per step drops from ~1.5 rad to 0.15–0.37 rad, but does not stop; steps 7 and 9
  failed their first Newton attempt (retry converged at 7; run stopped during the step-9 retry). Option 1 alone is not enough.

## Not done yet
* MPI (`hct_ddm.idp` createMatHCT; the `int1d` calls already use the physical labels), `make_mpi.py` patches, run scripts.
* The production polycrystal test.
* The KG regularisation (`-kwcreg 1`) in the ϑ terms.
* A check that the lagged scheme's first-order error in W is acceptable at dt = 0.01 (compare dt/2).
