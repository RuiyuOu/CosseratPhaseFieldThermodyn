---
title: "Stress and slip follow the lattice"
subtitle: "Cosserat–KWC crystal plasticity at finite deformation (production 7.5): statement, proof and numerical verification"
geometry: margin=1in
fontsize: 11pt
---

# 1. Statement

**Claim.** With the elasticity tensor and the slip systems rotated once by the initial orientation and never updated, the model's response refers to a single lattice at every point and time. That lattice has the current orientation $R_t = RR_0$. Three things follow:

1. The stress is the crystal's own elastic response, evaluated at the elastic stretch measured in the lattice frame and rotated by $R_t$.
2. The resolved shear stresses act on, and the slip occurs along, the slip systems of the current lattice.
3. The lattice carrying (1) and (2) is the same lattice whose orientation enters the KWC curvature.

All of this holds whatever made the lattice turn: rigid rotation, slip, the grain-boundary spin, or any combination, in any order.

# 2. Setting

**Kinematics.**
$$F=R\,\bar U F_p,\qquad \bar U=R^TFF_p^{-1},\qquad R_t=RR_0,\qquad R(X,0)=I,\quad F_p(X,0)=I .$$
Here $R$ is the relative Cosserat (lattice) rotation and $R_0(X)$ the initial orientation. Write $F_e:=FF_p^{-1}=R_eU_e$ (polar decomposition).

**Crystal data.** The data are rotated once, at setup, and never updated:
$$\mathbb C=R_0\star\mathbb C^{c},\qquad s_\alpha=R_0s^c_\alpha,\qquad n_\alpha=R_0n^c_\alpha,\qquad (Q\star\mathbb C)_{ijkl}=Q_{ia}Q_{jb}Q_{kc}Q_{ld}\,\mathbb C_{abcd}.$$

**Plastic flow.** Slip plus the grain-boundary spin:
$$L_p=\dot F_pF_p^{-1}=\sum_\alpha\dot\gamma_\alpha\,s_\alpha\otimes n_\alpha+\dot w^*\,\widehat W .$$

**Energy.**
$$\psi=\psi_{el}(\bar U)+2\lambda\cdot a+2\mu_c|a|^2+\psi_K(K,\phi,\dots),\qquad a=\operatorname{axl}(\operatorname{skw}\bar U),$$
where the elastic part is written in the Green strain $\bar E=\tfrac12(\bar U^T\bar U-I)$ (plus a volumetric term in $J=\det\bar U$). Its derivative is $\Pi=\partial\psi/\partial\bar U$.

**Stresses.**
$$P=R\,\Pi\,F_p^{-T},\qquad \tau=PF^T,\qquad M=\bar U^T\Pi,\qquad \tau_\alpha=M:(s_\alpha\otimes n_\alpha).$$

**Crystal-frame quantities.** For any reference-frame tensor $A$, write $A^c:=R_0^TAR_0$. Let $\psi^c$ be the same energy with $\mathbb C^c$ in place of $\mathbb C$, and $\Pi^c=\partial\psi^c/\partial\bar U^c$.

# 3. Three lemmas

**Lemma 1 (rotated material law).** $\psi(\bar U)=\psi^c(\bar U^c)$ and $\Pi=R_0\,\Pi^cR_0^T$.

*Proof.* Under $\bar U\mapsto R_0^T\bar UR_0$:

- $J$ is unchanged, and $\bar E\mapsto R_0^T\bar ER_0$.
- $(R_0\star\mathbb C^c):\bigl(R_0AR_0^T\bigr)=R_0(\mathbb C^c:A)R_0^T$, so $\tfrac12\bar E:\mathbb C:\bar E=\tfrac12\bar E^c:\mathbb C^c:\bar E^c$.
- $|a|$ is unchanged, since $\operatorname{axl}(Q^TWQ)=Q^T\operatorname{axl}W$.

The relation for $\Pi$ is the chain rule for the linear map $\bar U\mapsto R_0^T\bar UR_0$. $\square$

**Lemma 2 (no elastic torque).** $\bar U^T\bar U=F_e^TF_e$ for every $R$. Hence $\psi_{el}$ is independent of $R$, and the lattice rotation is set only by the constraint and curvature terms.

**Lemma 3 (the constraint selects the polar rotation).** Suppose $a=0$ and $\bar U\succ0$. Then $F_e=R\bar U$ is a polar decomposition of $F_e$, and by uniqueness
$$R=R_e=\operatorname{polar}(FF_p^{-1}),\qquad \bar U=U_e .$$
The branch $\bar U\succ0$ is selected by continuity from $R(0)=I=\operatorname{polar}(F_e(0))$. (See the constraint derivation, constraint_finite.pdf.) $\square$

# 4. Theorem: the stress follows the lattice

**Theorem 1.** For any $F$ and $F_p$,
$$\tau=R_t\,\mathcal T^c(\bar U^c)\,R_t^T,\qquad \mathcal T^c(U):=\Pi^c(U)\,U^T. \tag{1}$$
On the constraint ($a=0$), $\bar U^c=U_e^c$, where
$$U_e^c=R_t^TF_eR_0=R_0^TU_eR_0 \tag{2}$$
is the elastic stretch measured in the lattice frame. Equation (1) then reads: **the stress is the crystal's own response to its own elastic stretch, pushed forward by the current lattice orientation.**

*Proof.*
1. $\tau=PF^T=R\Pi F_p^{-T}F_p^T\bar U^TR^T=R\,\Pi\,\bar U^TR^T$.
2. Lemma 1 gives $\Pi=R_0\Pi^cR_0^T$ and $\bar U=R_0\bar U^cR_0^T$, so
$$\tau=RR_0\,\Pi^c\bar U^{cT}R_0^TR^T=R_t\,\mathcal T^c(\bar U^c)\,R_t^T .$$
3. On the constraint, Lemma 3 gives $\bar U=U_e$ and $R=R_e$. Hence $\bar U^c=R_0^TU_eR_0=R_0^TR_e^TF_eR_0=R_t^TF_eR_0$. $\square$

**Corollary 1 (lattice covariance).** Consider two states with the same lattice-frame stretch $U_e^c$ and lattice orientations $R_t$ and $R_t'$. Their stresses satisfy
$$\tau'=Q\,\tau\,Q^T,\qquad Q=R_t'R_t^T,$$
independently of $F$, $F_p$ and their histories.

**Corollary 2 (small elastic strain).** With $U_e^c=I+\varepsilon^c$ and the spatial strain $\varepsilon=R_t\varepsilon^cR_t^T$,
$$\tau\approx(R_t\star\mathbb C^c):\varepsilon .$$
The spatial stiffness is the crystal stiffness rotated by the current lattice orientation.

**Remark (the elastic stress follows $F_e$ even without the constraint).** For a Green-strain law $\Pi_{el}=\bar U\,S$ with $S=S(\bar U^T\bar U)$, so
$$\tau_{el}=R\bar US\bar U^TR^T=F_e\,S(F_e^TF_e)\,F_e^T .$$
This depends only on $F_e$, so the elastic stress always follows the lattice carried by $F_e$. What the constraint adds is that the Cosserat rotation $R$, which enters the curvature $K_i=\operatorname{axl}(R_t^T\partial_iR_t)$ and hence the grain-boundary energy, is that same lattice.

**Remark (skew stress).** The multiplier and the penalty add $R\,(2\lambda+4\mu_c a)\,(\partial a/\partial\bar U)\,\bar U^TR^T$ to $\tau$; in 2D, $\partial a/\partial\bar U=\tfrac12(e_2\otimes e_1-e_1\otimes e_2)$. This is the reaction to the curvature torque. It is expressed in the same lattice frame, and it vanishes when there is no curvature torque.

# 5. Theorem: slip follows the lattice

**Theorem 2.** On the constraint, with $s_\alpha$ and $n_\alpha$ fixed in the intermediate configuration:

1. **Driving force.** The resolved shear stress is the Kirchhoff stress resolved on the current lattice slip system:
$$\tau_\alpha=M:(s_\alpha\otimes n_\alpha)=\tau:\bigl(\tilde s_\alpha\otimes\tilde n_\alpha\bigr),\qquad
\tilde s_\alpha:=F_es_\alpha=R_tU_e^cs^c_\alpha,\quad \tilde n_\alpha:=F_e^{-T}n_\alpha=R_tU_e^{c\,-T}n^c_\alpha .\tag{3}$$
Equivalently, $\tau_\alpha=\Pi^c:(U_e^cs^c_\alpha\otimes n^c_\alpha)$, which is the crystal-frame expression.
2. **Kinematics.** Slip shears the material along the current slip direction on the current slip plane:
$$F_eL_pF_e^{-1}=\sum_\alpha\dot\gamma_\alpha\,\tilde s_\alpha\otimes\tilde n_\alpha+\dot w^*\,F_e\widehat WF_e^{-1}. \tag{4}$$

*Proof.*

(1) We have $M=\bar U^T\Pi$. Since $F_e=R\bar U$ and $\tau=R\Pi\bar U^TR^T$,
$$F_e^T\tau F_e^{-T}=\bar U^TR^TR\Pi\bar U^TR^TR\bar U^{-T}=\bar U^T\Pi=M .$$
Hence
$$M:(s\otimes n)=s\cdot Mn=s\cdot F_e^T\tau F_e^{-T}n=(F_es)\cdot\tau\,(F_e^{-T}n)=\tau:(F_es\otimes F_e^{-T}n).$$

The lattice-frame forms follow from $F_e=R_tU_e^cR_0^T$ (Theorem 1, Eq. 2) and $s_\alpha=R_0s^c_\alpha$, $n_\alpha=R_0n^c_\alpha$. The crystal-frame expression follows from $M=\bar U^T\Pi$ and Lemma 1.

(2) $F_e(s\otimes n)F_e^{-1}=(F_es)\otimes(F_e^{-T}n)$. $\square$

**Corollary 3.** The slip directions and normals must not be rotated in time. Doing so in addition to the $F_p$ update would resolve the stress on a slip system rotated twice. The same holds for $\mathbb C$. (See elasticity_invariance.pdf, Proposition 4, for the equivalent "rotating-lattice" bookkeeping, in which $F_p$ is updated with slip only.)

# 6. Why the lattice orientation is $R_t$: the plastic rotation

Theorem 2 of constraint_finite.pdf gives, exactly,
$$R=\operatorname{polar}(FF_p^{-1})=R_F\,R_{\mathrm{mix}}\,R_p^T,\qquad R_{\mathrm{mix}}=\operatorname{polar}(U_FU_p^{-1}),$$
where $F=R_FU_F$, $F_p=R_pU_p$, and $R_{\mathrm{mix}}=I$ for coaxial stretches. In 2D this is $\theta=\theta_F-\theta_p+\theta_{\mathrm{mix}}$. For a pure plastic rotation with $F=I$, the lattice turns by $R_p^T$ while the material stays put (Figure 1).

![Pure plastic rotation with $F=I$. $F_p=R(\theta_p)$ rotates the material, and the lattice keeps its orientation in the intermediate configuration. $F_e=FF_p^{-1}=R(-\theta_p)$ returns the material and carries the lattice. The constraint gives $R=\operatorname{polar}(FF_p^{-1})=R(-\theta_p)$. Red: material axes. Blue: lattice.](lattice_frame_figs/polar_rotation_figure_ms.png){width=100%}

The orientation $R_t(X)$ is spatial: it says where the crystal axes at material point $X$ point now. The fields themselves are labelled by the reference position $X$.

# 7. Numerical verification

**Setup.**

- **Material.** Cu ($C_{11}=168.4$, $C_{12}=121.4$, $C_{44}=75.4$ GPa), in the solver's lattice frame $x=[110]$, $y=[001]$, $z=[\bar110]$. The setup reproduces the code's $C_{xx}=220.3$, $C_{yy}=168.4$, $C_{xy}=121.4$ and $C_{ss}=75.4$ GPa.
- **Energy.** Plane-strain St Venant–Kirchhoff, $\psi=\tfrac12E:\mathbb C:E+2\mu_ca^2$ with $E=\tfrac12(\bar U^T\bar U-I)$, $\mu_c=754$ GPa.
- **Stress.** Production formulas: exact constraint $R=\operatorname{polar}(F_e)$, multiplier from the torque balance, and $P=R(\Pi+2\lambda\,\partial a/\partial\bar U)F_p^{-T}$. Derivatives by complex step (exact to round-off).
- **Load set.** 35 stretch states (uniaxial, biaxial, compression, pure shear, 20 random large stretches, $J\approx0.75$–$1.3$), each at 12 lattice angles $\alpha$ from $-100^\circ$ to $170^\circ$. The stress changes by up to 138% of its own size as $\alpha$ varies, so agreement is not trivial.
- **Reference.** A crystal that starts at orientation $\alpha$: $\mathbb C=R(\alpha)\star\mathbb C^c$, slip systems rotated by $\alpha$, $F=V$, $F_p=I$.

**Tests** ($\mathbb C$ and the slip systems fixed at orientation $0$, lattice turned by $\alpha$ during the history):

- **T1, rigid rotation:** $F=VR(\alpha)$, $F_p=I$.
- **T2, grain-boundary spin:** $F=V$, $F_p=R(-\alpha)$.
- **T3, slip plus spin:** $F_p$ = three slips ($\gamma=0.5,-0.2,0.3$ on systems at $0^\circ$, $60^\circ$, $120^\circ$) followed by $R(-\alpha)$, with $F=VR(\alpha)F_p$.
- **T4, time-stepped spin:** 200 steps, $F_p\leftarrow\exp(\Delta w^*\widehat W)F_p$, $w^*\to-73^\circ$, while $F$ goes from $I$ to $V$.
- **W, wrong:** T2 with $\mathbb C$ also rotated by $\alpha$.

| Quantity | T1 | T2 | T3 | T4 | W |
|---|---|---|---|---|---|
| Kirchhoff stress $\tau$ (max rel. error) | $1.4\times10^{-14}$ | $1.9\times10^{-14}$ | $1.6\times10^{-14}$ | $2.4\times10^{-13}$ | **1.70** |
| W compared with a crystal at $2\alpha$ | | | | | $2.3\times10^{-14}$ |

Slip systems, over T1–T3 and the full load set:

| Quantity | max error |
|---|---|
| Resolved shear stress $\tau_\alpha$ (relative) | $3\times10^{-13}$ |
| Current slip direction $F_es_\alpha$ | $7\times10^{-16}$ |
| Current slip-plane normal $F_e^{-T}n_\alpha$ | $7\times10^{-16}$ |
| $F_eL_pF_e^{-1}=\sum\dot\gamma_\alpha\,\tilde s_\alpha\otimes\tilde n_\alpha$ | $6\times10^{-17}$ |

The penalty was also tested with an applied couple, so that $a\neq0$ and the penalty carries stress ($|a|$ up to $2.3\times10^{-2}$). T1–T3 still agree to $\le6\times10^{-14}$; W does not.

![Stress against lattice rotation for $V=\operatorname{diag}(1.10,0.95)$. Lines: reference crystal at $\alpha$. Markers: tests with $\mathbb C$ fixed (T1–T4) and the wrong variant (W, which falls on the curve for $2\alpha$). Last panel: maximum error over the full load set.](lattice_frame_figs/c_fixed_all5.png){width=100%}

# 8. Scope

- **Exact versus approximate.** The statements are exact when $a=0$. In production the P0 multiplier makes the element mean of $a$ zero exactly; within an element, the $2\mu_ca^2$ penalty keeps $|\theta-\theta_e|$ of order $(\text{curvature torque})/(\mu_c(\operatorname{tr}U_e)^2)$.
- **Diffuse grain boundaries.** $\mathbb C$, $s_\alpha$ and $n_\alpha$ use each grain's $R_0$, taken from the mesh region. The curvature uses the smoothed $\theta_0$. Inside the diffuse boundary the two initial orientations differ.
- **The $\pi$ branch.** A lattice rotated by $\pi$ about an eigenvector of $U_e$ gives an identical stress and identical $\tau_\alpha$. "The lattice" is made unique by continuity from $R(0)=I$, and in 2D by the symmetry of the square lattice.
- **What the tests cover.** They are pointwise constitutive checks: no finite elements, no curvature, no time integration of the balance laws.

**Tests in the repository.** `checks/test_stress_follows_lattice.py` and `checks/test_slip_follows_lattice.py` (both laws, pass/fail, figures); see `checks/LATTICE_FRAME_TESTS.md`.
