# Derivative / smoothing model used for the derived animations

Mesh: P1 triangles. Solver output: nodal P1 fields (u, v, theta_dot, T, phi) and element-constant P0 fields (Cauchy stress, Fp, rhoG...).

1. Element gradient of a P1 field q_h = sum_i q_i N_i:  grad q_h|_e = (1/2A_e) sum_i q_i (b_i, c_i),
   b_i = y_j - y_k, c_i = x_k - x_j (cyclic). Constant per triangle, exact for the interpolant.
   Used for div v = d_x v_x + d_y v_y and curl v = d_x v_y - d_y v_x.
2. Nodal recovery of a P0 field (patch average = lumped L2 projection P0 -> P1):
   Q_i = sum_{e in patch(i)} A_e q_e / sum_{e in patch(i)} A_e.
3. Screened-Poisson (Helmholtz) smoothing, optional second stage, P1 FEM with natural BC:
   (M + l^2 K) Q~ = M Q, l = 0.04 um (default; argv[3] in animate.py).
   Equivalent to a Gaussian-like low-pass filter with transfer 1/(1+l^2 k^2); damps wavelengths below ~2 pi l.
4. GND from curl Fp (2D, in-plane slip, out-of-plane dislocation lines):
   (curl Fp)_{i3} = d_1 Fp_{i2} - d_2 Fp_{i1}; Fp components recovered by steps 2-3, gradient by step 1;
   rho_G = |curl Fp| / b, b = 2.556e-4 um, reported in m^-2. This is the Nye tensor magnitude / b (a lower bound on GND content; the
   scalar magnitude lumps edge and screw-like in-plane parts).
   The solver's own rhoG = phi^2 |grad theta - grad theta_R| / b is plotted alongside for comparison.
Known inaccuracies (accepted): smoothing blurs boundaries over ~l, constant-per-element gradients of a P1 field are discontinuous
across edges, and Fp is exp-map integrated so curl Fp includes mesh-scale noise in low-GND areas (visible as speckle at ~1e9-1e10 m^-2).
Wave quantities: pressure p = -(s11+s22)/2 and shear stress (s12+s21)/2 in MPa; div v, curl v in 1/ns; theta_dot in rad/ns; T-300 in K;
|u| in nm, shown on the mesh displaced by 400x u.
