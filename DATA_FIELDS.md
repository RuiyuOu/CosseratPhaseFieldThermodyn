# Data fields

Every field written by the solver (`cosserat_kwc.edp`) and the post-processing (`post/`), with its unit and the formula that produces it.
How to run the code is in `README.md`; where each formula sits in the solver is in `CODE_MAP.md`.

## Conventions

| | |
|---|---|
| Units | length µm, time ns, stress and energy density GPa, temperature K, density ρ in 10³ kg/m³. Velocities in µm/ns = km/s. Dislocation densities in µm⁻² (× 10¹² for m⁻²). Viscosities in GPa·ns = Pa·s. |
| Energies | per unit thickness (plane strain): GPa·µm² = 10⁻³ J/m |
| Domain | 10 × 10 µm square, 8 grains (triangle `Label` 1–8). Boundary labels: 1 left (clamped), 2 bottom, 3 right (plate), 4 top. "Domain means" divide by the area 100 µm². |
| Discretisation | displacement u: P2; microrotation θ, φ, T: P1 (nodes); history and stresses: P0 (one value per triangle) |
| Element-mean state | P0 stresses and driving forces are evaluated at the element-mean kinematics: the centroid value of the P2 gradient, written F̄. The energy and internal forces in the solve use 3 Gauss points instead. |
| Steps | `step` n counts every step, including the NRELAX zero-load relaxation steps. Time t = (n − NRELAX)·Δt; output starts at n = NRELAX (the relaxed state, t = 0). |
| Increments | Δx = x(n+1) − x(n) over one step; quantities marked "n+1" are end-of-step values. |

### Symbols used in the formulas

| Symbol | Meaning |
|---|---|
| θ | Cosserat microrotation relative to the initial grain orientation (rad); `dtheta_deg` |
| θ₀ | grain orientation, smoothed: (1 − ℓ²∇²)θ₀ = θ_grain with ℓ = 0.04 µm |
| θ_t = θ₀ + θ | total lattice orientation (`theta_lattice_deg`); ∇θ_t is its element gradient |
| F = I + ∇u = R(θ) Ū F_p | multiplicative split; R(θ) the rotation by θ, F_p the plastic part (det F_p = 1) |
| J = det F = det Ū | volume ratio |
| C̄ = ŪᵀŪ, Ē = (C̄ − I)/2 | elastic right Cauchy–Green tensor and Green strain |
| Ẽ = (J^(−2/3) C̄ − I)/2 | isochoric Green strain (with Ẽ₃₃ = (J^(−2/3) − 1)/2 in plane strain) |
| a = (Ū₂₁ − Ū₁₂)/2 | skew part of Ū: the mismatch between lattice and microrotation |
| ψ | elastic energy density: ψ = W_V(J) + ½ Ẽ:C:Ẽ + 2μ_c a² |
| W_V | Vinet: W_V = 9K₀/η² [1 − (1 − η(1 − x)) e^(η(1 − x))], x = J^(1/3), η = 1.5(K′ − 1); pressure P = 3K₀ (1 − x)/x² e^(η(1 − x)) |
| Π = ∂ψ/∂Ū, M = ŪᵀΠ | stress conjugate to Ū and the Mandel stress |
| s_a, n_a | slip direction and normal of system a (n = s rotated by +90°); lattice directions 0°, ±54.74° from the lattice x axis, rotated by the grain angle |
| f_a | Rice planar triple-slip weights (√3, 2/√3, 2/√3) for systems 1, 2, 3 |
| Ŵ | the 2-D rotation generator [[0, −1], [1, 0]] |
| τ_a = M : (s_a ⊗ n_a) | resolved shear stress of system a |
| τ* = M₂₁ − M₁₂ | driving force of the grain-boundary spin w* |
| φ | KWC order parameter (1 in the grain, lower in the boundary core) |
| g(φ) = −2[φ + ln(1 − φ)] | KWC coupling; ḡ is its element mean |
| N(∇θ_t) = √(\|∇θ_t\|² + δ²) − δ | regularised orientation-gradient norm |
| ρ_S, ρ_G, ρ_M | statistically stored, geometrically necessary and evolved mobile dislocation densities (µm⁻²) |
| ρ_m = max(min(ρ_M, ρ_S + ρ_G), 10⁻³) | mobile density used in the drag |
| Γ = Σ_a f_a \|Δγ_a\| | weighted slip increment (drives storage and ρ_M) |
| Z = ½ χ μ b² φ | stored-energy factor per unit density |
| Parameters | μ = 48 GPa, b = 2.556×10⁻⁴ µm, τ₀ = 0.001 GPa, λ_T = 0.3, ρ₀ = 10 µm⁻², γ̇₀ = 10⁻⁶ ns⁻¹, m = 50, B₀ = 2×10⁻⁵ GPa·ns, K = 0.05, d = 5×10⁻⁴ µm, χ = 1, T₀ = 300 K, α_μ = 3.8×10⁻⁴ K⁻¹, ρC = 3.45×10⁻³ GPa/K, k = 4×10⁻⁴ GPa·µm²/(ns·K), ρ = 8.96. The rest (f₀, a, s, δ, η_φ, Q_φ, η*₀, C_D, C_A, v_c, ρ_Msat, C₂, ρ_Fcut, ρ_Mmin, μ_c, j, K₀, K′, softmin) are the solver options in `README.md` and in `<tag>_meta.json`. |

### Constitutive updates referenced below

- **Strength:** τ_c = (τ₀ + λ_T μ b √(ρ_S + ρ_G)) · max(1 − α_μ (T − T₀), softmin). It is evaluated once at the start of each step, from the end state of the previous step (T at the element centroid). The same holds for ρ_m, η_d and the relativistic rate limit below.
- **Drag:** η_d = B₀ (T/T₀)/(b² ρ_m). The relativistic effective drag is η_eff = η_d √(1 + s²) with s = \|τ\|/(η_d ρ_m b v_c), so \|γ̇\| ≤ ρ_m b v_c.
- **Slip rate:** thermal activation in series with drag. γ̇_a = sign(τ_a) x^m / (1/γ̇₀ + η_eff x^(m−1)/(f_a τ_c)), where x = min(\|τ_a\|/(f_a τ_c), 300). Backward Euler: Δγ_a = Δt γ̇_a(τ_a,n+1).
- **Grain-boundary spin:** Δw* = Δt τ*/η_s with η_s = 2η*₀/(1 − φ)², and \|Δw*\| ≤ `spinbranch` per step.
- **Plastic deformation:** F_p,n+1 = exp(Σ_a Δγ_a s_a ⊗ n_a + Δw* Ŵ) F_p,n, using the exact exponential.
- **Stored dislocations:** Δρ_S = (K/b)√ρ_S Γ − (2d/b) ρ_S Γ − C_D max(Δφ, 0) tanh(C_A² \|∇θ_t\|²) ρ_S.
  - The first term is storage, the second dynamic recovery, the third recovery by a sweeping boundary.
  - Storage is capped so that Z Δρ_stored ≤ Σ_a \|τ_a Δγ_a\| (Clausius–Duhem), and ρ_S ≥ 10⁻³.
- **Mobile density:** ρ_M* = max(ρ_Msat (1 − √((ρ_S + ρ_G)/ρ_Fcut)), ρ_Mmin). Then ρ_M,n+1 = max((ρ_M,n + C₂ Γ ρ_M*)/(1 + C₂ Γ), ρ_Mmin).
- **GND (Cermelli–Gurtin):** G = F_p (Curl F_p).
  - The curl is c_i = ∂₁F_p,i2 − ∂₂F_p,i1, computed as a weak nodal curl (lumped P1 mass), then averaged over the element.
  - ρ_G = \|G − G_R\|/b, where G_R is G at the end of the relaxation.
- **Dissipation per element and step** (energy density, GPa):
  - D = Σ_a τ_a,n+1 Δγ_a + τ* Δw* + η_φ A(T) (Δφ)²/Δt − Z Δρ_S;
  - A(T) = max(exp(Q_φ/R (1/T − 1/T₀)), `etaphimin`);
  - the first two terms are ≥ 0 for a converged backward-Euler step.
- **Heat:** ρC (T − T_n)/Δt − k ∇²T = D/Δt, with natural (insulated) boundaries.
- **Phase field:** a backward-Euler gradient flow of the KWC energy plus the stored energy, with 0 ≤ φ ≤ 1 − 10⁻⁷ enforced exactly.

## Files at a glance

| File (in `out_<tag>/`) | Written | Content |
|---|---|---|
| `<tag>_s<k>.vtu` | every VTKEVERY loaded steps (k = n − NRELAX) | full field snapshot (below) |
| `<tag>_t<k>.vtu` | at t = 2, 5, 7, 10 ns and t_end (k = 100·t) | same fields |
| `<tag>_relaxed.vtu` | end of the relaxation (t = 0) | same fields; the reference frame for the animations |
| `<tag>_stall<n>.vtu`, `<tag>_tstop<n>.vtu` | at a non-converged step / a TSTOP | same fields, last committed state |
| `<tag>.pvd` | with every VTU | ParaView time series: `timestep` = t (ns) and the file name |
| `<tag>_log.csv` | every step from n = NRELAX | energies, energy audit, extremes |
| `<tag>_thermo.csv` | every step from n = NRELAX | dissipation budget of the step |
| `<tag>_diag.csv` | every step, 3 rows | the worst elements of the step |
| `<tag>_meta.json` | at the start | every parameter, units, lattice, grain angles |
| `<tag>.args` | by `run.sh` | the option list the run was started with (resume guard) |
| `<tag>_ck.txt` (+ `.manifest`, `.gnd`, `.plate`, `.rhoM`), `<tag>_ckroll{0,1}.txt`, `<tag>_relaxed_ck.txt` | checkpoints | restart state (not for analysis) |
| `<tag>_TSTOP.txt`, `_NOCONV.txt`, `_NAN.txt`, `STALLED` | on a stop | one line: why and where the run stopped |
| `<tag>_stage.txt` | always (overwritten) | the solver stage in progress (left behind by a hard crash) |

## VTU snapshots

### Point data (P1, one value per node)

| Field | Unit | Definition |
|---|---|---|
| `u` (3 comp., z = 0) | µm | displacement u at the mesh vertices |
| `v` (3 comp., z = 0) | µm/ns | velocity u̇ (Newmark) |
| `dtheta_deg` | deg | microrotation θ relative to the initial grain orientation, including the rotation during relaxation |
| `theta_lattice_deg` | deg | lattice orientation θ_t = θ₀ + θ |
| `theta_dot` | rad/ns | microrotation rate θ̇ |
| `T` | K | temperature |
| `phi` | – | KWC order parameter φ |

### Cell data (P0, one value per triangle)

All stresses are at the element-mean state F̄ of the committed step, with zero increments.

| Field | Unit | Definition |
|---|---|---|
| `Label` | – | grain number 1–8 on triangles; boundary number 1–4 on the boundary line cells (added by FreeFem) |
| `p1`, `p2`, `p3` | – | accumulated signed slip γ_a = Σ_steps Δγ_a |
| `slip_path_eff` | – | Σ_steps Σ_a \|Δγ_a\| (total slip path, unweighted) |
| `slip_path_phys` | – | Σ_steps Γ = Σ_steps Σ_a f_a \|Δγ_a\| (slip path with the planar-triple-slip weights) |
| `rhoS` | µm⁻² | ρ_S |
| `rhoG` | µm⁻² | ρ_G = \|G − G_R\|/b (Cermelli–Gurtin, G = F_p Curl F_p); 0 at the end of the relaxation |
| `kp_x`, `kp_y` | – | always 0 (kept for the checkpoint and VTU layout) |
| `wstar` | rad | accumulated grain-boundary spin w* = Σ_steps Δw* |
| `tau_c` | GPa | τ_c used in the last step: from ρ_S, ρ_G and T of the previous step, i.e. one step behind `rhoS`, `rhoG` and `T` in the same file |
| `tau_1`, `tau_2`, `tau_3` | GPa | τ_a = M : (s_a ⊗ n_a), signed |
| `sig_11`, `sig_22`, `sig_12`, `sig_21` | GPa | Cauchy stress σ = P Fᵀ/J with P = ∂ψ/∂F = R Π F_p⁻ᵀ (σ_ij: i = force direction; σ₁₂ ≠ σ₂₁ in a Cosserat medium) |
| `sig_33` | GPa | out-of-plane Cauchy stress S₃₃/J (S₃₃ from the Vinet + isochoric law, F₃₃ = 1) |
| `m_x`, `m_y` | GPa·µm | KWC couple stress m = f₀ s ḡ ∇θ_t / √(\|∇θ_t\|² + δ²) |
| `D_rate` | GPa/ns | D/Δt of the last step (power density; × 10¹⁸ for W/m³) |
| `Fp_11`, `Fp_12`, `Fp_21`, `Fp_22` | – | F_p (global reference-lattice frame) |
| `w_cosserat` | rad | a = (Ū₂₁ − Ū₁₂)/2 |
| `J` | – | det F̄ |
| `psi_elastic` | GPa | ψ = W_V(J) + ½ Ẽ:C:Ẽ + 2μ_c a² at F̄ (KWC and stored energies not included) |
| `rho_mob` | µm⁻² | ρ_m = max(min(ρ_M, ρ_S + ρ_G), 10⁻³) used in the last step (previous-step state, like `tau_c`); the evolved ρ_M itself is only in the checkpoint `.rhoM` file |
| `eta_drag` | GPa·ns | low-stress drag η_d = B₀ (T/T₀)/(b² ρ_m) used in the last step (previous-step state; before the relativistic factor) |

## `<tag>_log.csv`

One row per step from n = NRELAX. After a restart, `run.sh` trims rows past the checkpoint, so no step appears twice.

| Column | Unit | Definition |
|---|---|---|
| `step` | – | step number n (relaxation included) |
| `t_ns` | ns | t = (n − NRELAX)·Δt |
| `load_GPa` | GPa | always 0 (no traction loading; kept for the column layout) |
| `iters` | – | global Newton iterations of the accepted solve (of the retry, if there was one) |
| `locres` | – | largest local residual \|Δγ_a − Δt γ̇_a\| or spin residual over all elements, last local iteration |
| `KE_trans` | GPa·µm² | ∫ ½ ρ \|u̇\|² dA |
| `KE_rot` | GPa·µm² | ∫ ½ ρ j θ̇² dA |
| `E_elastic` | GPa·µm² | Gauss-point elastic energy Σ_e A_e Σ_q w_q ψ(F_q) plus the KWC energy ∫ f₀ [½(1 − φ)² + ½ a² \|∇φ\|² + s g(φ) N(∇θ_t)] dA |
| `E_stored` | GPa·µm² | ∫ ½ χ μ b² φ ρ_S dA |
| `Q_heat` | GPa·µm² | ∫ ρC (T − T₀) dA |
| `W_ext` | GPa·µm² | always 0 (no traction loading); the plate work is `W_plate` |
| `D_cum` | GPa·µm² | ∫ Σ_steps D dA since the end of the relaxation |
| `audit_residual` | GPa·µm² | W_ext − (E − E₀) = −(E − E₀), E = KE_trans + KE_rot + E_elastic + E_stored + Q_heat. E₀ is E evaluated during the reset at the end of the relaxation, just before the state of the first row is logged, so the first row can show a small offset (≈ 10⁻⁶) that carries through. |
| `max_abs_slip` | – | max over elements of Σ_a \|γ_a\| |
| `mean_abs_slip` | – | (1/100 µm²) ∫ Σ_a \|γ_a\| dA |
| `max_dT` | K | max T − T₀ |
| `max_abs_dtheta_deg` | deg | max over nodes of \|θ\| |
| `mean_rhoS_um2` | µm⁻² | (1/100 µm²) ∫ ρ_S dA |
| `max_rhoG_um2` | µm⁻² | max ρ_G |
| `min_phi` | – | min φ |
| `neg_D_count` | – | element-steps with D < −10⁻¹⁴ since the end of the relaxation (cumulative; should stay 0) |
| `W_plate` | GPa·µm² | work done on the body through the prescribed dofs: Σ_steps ½(f_c,n + f_c,n+1)·(U_c,n+1 − U_c,n), f_c = −R (reaction of the discrete equations, inertia included). Only the plate moves, so this is the plate work. |
| `audit_full` | GPa·µm² | energy balance: audit_residual − R₀ + W_plate = W_plate − (E − E₀) − R₀. Dissipation becomes heat, which is part of E, so an exactly conserving scheme would keep it at 0. The backward-Euler slip and phase-field updates dissipate slightly more than the D they report, so it drifts slowly; watch its size relative to W_plate. R₀ = 0, except after a restart without a `.plate` file, where the ledger restarts from the first step (that row shows `nan`). |

## `<tag>_thermo.csv`

The dissipation budget of each step. Sums run over the elements e, weighted by area A_e.

| Column | Unit | Definition |
|---|---|---|
| `step`, `t_ns` | – , ns | as in the log |
| `D_heat` | GPa·µm² | Σ_e D_e A_e: the heat source of the step |
| `D_implicit_slip` | GPa·µm² | Σ_e Σ_a τ_a,n+1 Δγ_a A_e |
| `D_midpoint_slip` | GPa·µm² | Σ_e Σ_a ½(τ_a,n + τ_a,n+1) Δγ_a A_e (the value consistent with the trapezoidal energy identity) |
| `alg_mid_minus_imp_cum` | GPa·µm² | cumulative Σ_steps (D_midpoint_slip − D_implicit_slip): the algorithmic dissipation, O(Δt²) per step, accumulated from the start of the run, relaxation included (from 0 again after a restart) |
| `neg_D_elems` | – | elements with D < −10⁻¹⁴ in this step |
| `min_D` | GPa | min_e D_e |
| `stored_energy_rate` | GPa·µm²/ns | Σ_e Z Δρ_S A_e / Δt |
| `cdstore` | – | always 1: storage is Clausius–Duhem capped |

## `<tag>_diag.csv`

Three rows per step, one per "worst" element, selected by the value in `tag`. All values are at the end of the step.

| Column | Unit | Definition |
|---|---|---|
| `step`, `t_ns` | – , ns | as in the log |
| `tag` | – | `maxslipinc` (largest Σ_a \|Δγ_a\|), `maxlocres` (largest local residual), `maxD` (largest D) |
| `elem` | – | triangle index |
| `grain` | – | its grain label 1–8 |
| `x`, `y` | µm | centroid |
| `dist_edge` | µm | distance of the centroid to the nearest domain edge |
| `phi` | – | φ at the centroid |
| `grad_theta_tot` | rad/µm | \|∇θ_t\| |
| `rhoG`, `rhoS` | µm⁻² | ρ_G, ρ_S |
| `slip_inc` | – | Σ_a \|Δγ_a\| in this step |
| `locres` | – | max_a \|Δγ_a − Δt γ̇_a\| |
| `tau_star` | GPa | τ* |
| `dw_star` | rad | Δw* in this step |
| `tau_c` | GPa | τ_c used in this step (previous-step state) |
| `eta_drag` | GPa·ns | η_d |
| `D` | GPa | D of the element in this step |
| `max_dT` | K | max T − T₀ (domain) |
| `rho_mob` | µm⁻² | ρ_m |
| `eta_eff` | GPa·ns | max_a η_eff,a (relativistic drag actually used) |
| `tau_abs_max` | GPa | max_a \|τ_a\| |
| `v_kms` | km/s | slip_inc/(Δt ρ_m b): mean dislocation speed from Orowan's relation |

## `<tag>_meta.json`

Written once at the start; post-processing reads parameters from it.

| Key | Content |
|---|---|
| `units` | the unit system above |
| `mesh`, `j_um2`, `load`, `plate_un_um_ns`, `plate_ut_um_ns`, `ramp_ns`, `hold_ns`, `tend_ns`, `dt_ns`, `nrelax` | run settings: micro-inertia j (µm²), plate velocity components (normal −x, tangential +y, µm/ns), quintic ramp and hold times |
| `mu_c`, `rho`, `C_lattice`, `lattice_frame` | Cosserat modulus, density, cubic constants in the lattice frame x = [110], y = [001], z = [−110] |
| `slip_systems_lattice` | s, n and weight f for each system |
| `grain_theta0_rad` | the 8 grain angles (from `grain_angles_rad.txt`) |
| `hardening`, `rate`, `mobile_drag`, `kwc`, `thermal`, `elastic_law` | the parameters in the formulas above |
| `safeguards`, `robustness`, `tangent`, `newton` | solver settings |
| `kinematics`, `fields_note` | text description of the kinematics and of the VTU fields |

## Post-processing outputs

### Animations (`anim_<tag>/<field>/png/`, `<field>.gif`; packaged into `deliver_<tag>/`)

| Field | Unit | Definition |
|---|---|---|
| `phi` | – | φ |
| `pressure` | MPa | −(σ₁₁ + σ₂₂)/2 (in-plane mean; the tau statistics below use −(σ₁₁ + σ₂₂ + σ₃₃)/3) |
| `shear` | MPa | (σ₁₂ + σ₂₁)/2 |
| `stress_frob` | MPa | √(σ₁₁² + σ₂₂² + σ₁₂² + σ₂₁² + σ₃₃²) |
| `divv`, `curlv` | 1/ns | ∇·v and ∂ₓv_y − ∂ᵧv_x (constant per triangle, from the P1 nodal velocity) |
| `microrot` | rad/ns | θ̇ |
| `temperature` | K | T − 300 |
| `deformation` | nm | \|u\|, drawn on the mesh displaced by a magnification shown in the title |
| `slip_acc` | – | `slip_path_eff` (log scale) |
| `slip_sys1..3` | – | `p1`, `p2`, `p3` |
| `rotation` | deg | θ_t − θ_t(relaxed frame), wrapped to [−180, 180) |
| `lattice_slip` | deg | θ_t wrapped to [−90, 90), with the three slip directions overlaid (line width ∝ \|γ_a\|) |
| `gnd` | m⁻² | `rhoG` × 10¹² (log scale) |
| `rho_mob` | m⁻² | `rho_mob` × 10¹² (log scale) |
| `eta_drag` | Pa·s | `eta_drag` (log scale) |
| `grain_theta`, `grain_territory`, `grain_misorientation` | deg | grain segmentation per frame (see Grains) |

The colour limits are set from 12 frames spread over the run and stored in `limits.json`; `manifest.json` lists the frames already drawn. Log-scale fields are clipped at their lower limit. `deliver_<tag>/` holds, for each field, `<field>.gif`, `<field>_contact_sheet.png` (12 frames) and `<field>_pngs.zip`, plus `<tag>_data.zip` (the run folder without `.out`).

### Stress statistics: `deliver_<tag>/tau/`

`tau_stats_<tag>.csv` has one row per snapshot and per element selection:

- `all`: every element;
- `bulk`: element-mean φ ≥ 0.5;
- `core`: element-mean φ < 0.5.

Percentiles and means are area-weighted. A weighted quantile interpolates the cumulative area at the element midpoints. τ_max = max_a \|τ_a\| per element.

| Column | Definition |
|---|---|
| `t_ns`, `snapshot`, `selection` | time (from the `.pvd`), file tag, selection |
| `n_elem`, `area_um2` | number of elements and their area |
| `taumax_max`, `_p99`, `_p95`, `_p90`, `_median`, `_mean` | statistics of τ_max (GPa) |
| `tau<a>_absmax`, `_absp99`, `_absmean`, `_min`, `_max` | statistics of \|τ_a\| and the signed extremes of τ_a (GPa) |
| `tauc_min`, `_median`, `_mean`, `_max` | statistics of τ_c (GPa) |
| `ratio_median`, `ratio_p99` | statistics of τ_max/τ_c |
| `area_frac_tau_ge_tauc` | area fraction with τ_max ≥ τ_c |
| `area_frac_taumax_ge_<x>` | area fraction with τ_max ≥ x GPa (x = 0.1, 0.2, 0.5 by default) |
| `T_max`, `T_p99`, `T_median`, `T_mean` | element-mean temperature (K) |
| `pressure_max`, `_p99`, `_median`, `_p01`, `_min` | p = −(σ₁₁ + σ₂₂ + σ₃₃)/3 (GPa, compression positive) |
| `beta_absmax_deg`, `_absp99_deg`, `_absmedian_deg` | Cosserat mismatch β = φ_e − θ̄, wrapped to [−180, 180): G = F F_p⁻¹ with F from the P1 gradient of the vertex u, φ_e = atan2(G₂₁ − G₁₂, G₁₁ + G₂₂), θ̄ = element mean of `dtheta_deg` |
| `dtheta_absmax_deg`, `dtheta_absp99_deg` | \|θ̄\| (deg) |

`tau_summary.csv` has one row per run, computed from the `all` rows:

| Column | Definition |
|---|---|
| `run`, `speed_m_s`, `angle_deg`, `n_snapshots` | run tag, plate speed and angle parsed from the tag, number of snapshots |
| `peak_<q>`, `t_of_peak_<q>` | maximum over time and its time, for q = taumax_max, taumax_p99, taumax_median, tauc_median, T_max, beta_absmax_deg, pressure_p99; `peak_beta_absp99_deg`; `peak_dtheta_absmax_deg_vtu` and its time |
| `t_of_largest_mean_T_rise`, `largest_mean_T_rise_K_per_ns` | the largest increase rate of T_mean between consecutive snapshots, at the middle of that interval |
| `t_of_largest_T_max_rise`, `largest_T_max_rise_K_per_ns` | the same for T_max |
| `window_mean_<q>`, `window` | mean over the snapshots in the loading window (TAUWIN0–TAUWIN1, default 0.5–30 ns), for q = taumax_max, taumax_p99, taumax_median, tauc_median |
| `log_peak_max_abs_dtheta_deg`, `log_t_of_peak_max_abs_dtheta_deg` | peak of `max_abs_dtheta_deg` from `<tag>_log.csv` and its time |

The figures are `tau_vs_time.png` (all panels) and each panel alone: `tau_stress`, `tau_slip_systems`, `tau_temperature`, `tau_mismatch`.

### Grains: `deliver_<tag>/grains/`

**Segmentation** (`grain_count.py`, applied to `theta_lattice_deg` of every frame):

1. Plateau seeds are connected triangles with \|∇θ_t\| < GRAD (2°/µm) and area ≥ MINAREA (0.1 µm²).
2. Every other triangle joins the nearest seed by graph distance along the mesh.
3. Neighbouring plateaus whose mean orientations differ by less than MERGE (2°) are merged.

A grain's **plateau** is its seed triangles; its **territory** is everything assigned to it, and the territories tile the domain. A grain's orientation θ_g is the area-weighted mean of θ_t over its plateau.

Misorientation between neighbours i and j:
- raw: \|θ_j − θ_i\|;
- sym: reduced by the cubic period SYM = 90°, \|((θ_j − θ_i + 45) mod 90) − 45\|.

The bins are <2°, [2, 5) very low angle (VL), [5, 15) low angle (L) and ≥15° high angle (H).

`grain_history.csv`: one row per frame.

| Column | Definition |
|---|---|
| `t_ns`, `step`, `file` | frame time, step, VTU file |
| `n_grains` | number of grains |
| `n_ge1um2` | grains with plateau area ≥ 1 µm² |
| `n_boundaries` | neighbouring grain pairs (territories sharing mesh edges) |
| `boundary_length_um` | total length of the shared edges |
| `plateau_area_frac` | plateau area / domain area |
| `mean_grain_area_um2` | mean plateau area |
| `pct_raw_<bin>`, `pct_sym_<bin>` | % of boundary length in each misorientation bin (raw, sym) |
| `mean_raw_deg`, `mean_sym_deg` | length-weighted mean misorientation |

`grain_table_all.csv`: one row per grain and frame.

| Column | Definition |
|---|---|
| `t_ns`, `step`, `file`, `id` | frame and grain id |
| `plateau_area_um2`, `territory_area_um2` | areas |
| `cx_um`, `cy_um` | area-weighted plateau centroid |
| `theta_deg`, `theta_std_deg` | θ_g and the area-weighted standard deviation of θ_t over the plateau |
| `theta_sym_deg` | θ_g reduced to [−45, 45) |
| `orig_label` | the initial grain (`Label`) covering most of the plateau |
| `n_neighbours` | number of neighbouring grains |

`boundary_table_all.csv`: one row per neighbour pair and frame.

| Column | Definition |
|---|---|
| `grain_i`, `grain_j`, `theta_i_deg`, `theta_j_deg` | the pair and their orientations |
| `misorient_deg`, `misorient_sym_deg` | raw and reduced misorientation |
| `length_um` | shared edge length |
| `width_um` | median over the shared edges of (distance to the nearest plateau triangle of i + of j): the width of the transition |
| `mean_phi_at_edge` | mean nodal φ on the shared edges |
| `mid_x_um`, `mid_y_um` | midpoint of the shared edges |

Other files in the folder:
- `grains_first_*` and `grains_last_*`: the same tables for the relaxed frame and the last frame (`_grains.csv`, `_boundaries.csv`), plus `_misorient_hist.csv` (length, % and count per bin, raw and sym) and figures.
- `state.json`: the segmentation parameters and the frames already analysed.

`grains_first_misorient_hist.csv`, `grains_last_misorient_hist.csv`: one row per misorientation bin.

| Column | Definition |
|---|---|
| `bin_deg` | `<2`, `[2,5)`, `[5,15)`, `>=15` |
| `raw_length_um`, `raw_percent`, `raw_n` | boundary length, % of the total length, number of boundaries in the bin (raw misorientation) |
| `sym_length_um`, `sym_percent`, `sym_n` | the same for the reduced misorientation |

### Grain summary: `grain_summary[<suffix>]/summary_table.csv`

One row per run. The equivalent diameter of a grain is d = √(4A/π), with A its territory area (`--area territory`, the default) or its plateau area. d is the mean over grains, and the subscript 0 refers to the first (relaxed) frame.

| Column | Definition |
|---|---|
| `tag`, `angle`, `up_ms` | run |
| `t_end` | time of the last frame |
| `N0`, `N_end`, `N_peak`, `N_ratio_end` | grain count at the start, at the end, its maximum, N_end/N0 |
| `d0_um`, `d_end_um`, `d_std_end_um` | mean d at the start and end, standard deviation over grains at the end |
| `refine_end`, `refine_max`, `t_refine_max` | refinement factor d0/d at the end, its maximum and the time of that maximum |
| `VL_end`, `L_end`, `H_end`, `lt2_end` | % of boundary length in [2, 5), [5, 15), ≥15 and <2° at the end (sym by default, `--misor raw` for raw) |

`end_grid.png` shows the end values over angle × speed; with `--tavg N` they are averaged over the last N ns instead of the last frame.
