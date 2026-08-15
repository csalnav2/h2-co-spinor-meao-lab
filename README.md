# EAO Spinor Lab: Core Equation Spine

## Scope

The project is organized mathematically as

```math
\text{molecular envelope}
\rightarrow
\text{spinor field}
\rightarrow
\text{structured-light evolution}
\rightarrow
\text{bond and information diagnostics}
\rightarrow
\text{Clifford descriptors}
\rightarrow
\text{PCA/ridge-DMD}
\rightarrow
S^2/S^3\text{ topology}.
```

Sections 1 through 14 describe equations implemented in the current reduced-order framework. Section 15 contains proposed future braid experiments and is not yet an implemented entanglement measurement.

> **Reconstruction note**
>
> The source text had stripped display-math delimiters, missing equation labels, malformed matrix row separators, and several incomplete expressions. This version restores GitHub-compatible Markdown and labels equations from their surrounding section context. Sections 11.1 and 11.2 remain explicitly marked as partial because the uploaded source omitted part of those state definitions.

---

## 1. Grid, inner product, and normalization

### 1.1 Spatial grid

For a square domain of side length $L$ containing $N\times N$ grid points,

```math
\Delta x=\frac{L}{N},
```

and

```math
x_i=-\frac{L}{2}+i\Delta x,
\qquad
i=0,\ldots,N-1.
```

The grid is periodic because the finite-difference operators use wrapped neighbors.

### 1.2 Discrete inner product

For two scalar fields $a$ and $b$,

```math
\langle a,b\rangle_h
=
\sum_{i,j}
a_{ij}^{*}b_{ij}
(\Delta x)^2.
```

### 1.3 Scalar-orbital normalization

A scalar orbital is normalized by

```math
\phi
\longrightarrow
\frac{\phi}
{
\sqrt{
\sum_{i,j}|\phi_{ij}|^2(\Delta x)^2+\varepsilon
}
}.
```

Thus,

```math
\langle\phi,\phi\rangle_h\approx 1,
```

up to the regularizer $\varepsilon$ and floating-point error.

### 1.4 Two-component spinor

The simulated field is

```math
\Psi(x,y,t)
=
\begin{pmatrix}
\psi_{\uparrow}(x,y,t)\\
\psi_{\downarrow}(x,y,t)
\end{pmatrix}
=
\begin{pmatrix}
u(x,y,t)\\
d(x,y,t)
\end{pmatrix}.
```

Its density is

```math
\rho(x,y,t)=|u|^2+|d|^2.
```

The integrated field mass is

```math
M(t)
=
\sum_{i,j}
\rho_{ij}(t)
(\Delta x)^2.
```

### 1.5 Mean-density normalization

The current full model sets the initial target mass to

```math
M_{\mathrm{target}}
=
\rho_{\mathrm{target}}A,
\qquad
A=L^2.
```

The rescaling is

```math
\Psi
\longrightarrow
\Psi
\sqrt{
\frac{
\rho_{\mathrm{target}}A
}{
M+\varepsilon
}
}.
```

Equivalently,

```math
\frac{1}{A}
\sum_{i,j}
\rho_{ij}(t)(\Delta x)^2
\approx
\rho_{\mathrm{target}}.
```

This is a mean-field amplitude convention, not a one-electron probability normalization. The conservative launch disables repeated normalization after every time step.

These definitions are implemented directly in the grid, inner-product, scalar-normalization, and spinor-normalization functions.

---

## 2. Analytic H₂/CO molecular scaffold

### 2.1 Atomic-center locations

For bond coordinate $R$,

```math
x_A=-\frac{R}{2},
\qquad
x_B=+\frac{R}{2}.
```

### 2.2 Softened radial coordinates

For atom $A$,

```math
r_A
=
\sqrt{
(x-x_A)^2+(y-y_A)^2+s_{\mathrm{AO}}^2
}.
```

Likewise,

```math
r_B
=
\sqrt{
(x-x_B)^2+(y-y_B)^2+s_{\mathrm{AO}}^2
}.
```

The implementation uses

```math
s_{\mathrm{AO}}=0.08.
```

### 2.3 Slater-like atomic envelopes

```math
\phi_A
=
\mathcal N_A e^{-\zeta_A r_A},
```

```math
\phi_B
=
\mathcal N_B e^{-\zeta_B r_B},
```

where $\mathcal N_A$ and $\mathcal N_B$ are numerical normalization factors.

### 2.4 Atomic-orbital overlap

```math
S
=
\operatorname{Re}
\langle\phi_A,\phi_B\rangle_h.
```

### 2.5 Bonding channel

```math
\phi_b^{\mathrm{raw}}
=
c_A\phi_A+c_B\phi_B,
```

```math
\phi_b
=
\frac{
\phi_b^{\mathrm{raw}}
}{
\sqrt{
\langle\phi_b^{\mathrm{raw}},
\phi_b^{\mathrm{raw}}\rangle_h
+\varepsilon
}
}.
```

### 2.6 Antibonding channel

```math
\phi_{ab}^{\mathrm{raw}}
=
c_B\phi_A-c_A\phi_B,
```

```math
\phi_{ab}
=
\frac{
\phi_{ab}^{\mathrm{raw}}
}{
\sqrt{
\langle\phi_{ab}^{\mathrm{raw}},
\phi_{ab}^{\mathrm{raw}}\rangle_h
+\varepsilon
}
}.
```

### 2.7 Reduced molecular potential

Define softened potential distances

```math
r_A^{(V)}
=
\sqrt{
(x-x_A)^2+y^2+s_V^2
},
```

```math
r_B^{(V)}
=
\sqrt{
(x-x_B)^2+y^2+s_V^2
},
```

with

```math
s_V=0.18.
```

The reduced two-center potential is

```math
V_{\mathrm{mol}}
=
-V_0
\left[
\frac{Z_A}{Z_A+Z_B}\frac{1}{r_A^{(V)}}
+
\frac{Z_B}{Z_A+Z_B}\frac{1}{r_B^{(V)}}
\right]
+
0.045\,p\,x,
```

where $p$ is the heteronuclear polarity parameter. For H₂, $p=0$; for CO, the nonzero polarity term breaks left-right symmetry.

These molecular envelopes, overlaps, channels, and reduced potentials are implemented as analytic two-center approximations.

---

## 3. Initial spinor texture

Define the initial spin imbalance

```math
\delta_z=0.10p,
```

and relative phase

```math
\theta_{\mathrm{rel}}(x,y)
=
0.35\tanh y+0.15px.
```

The initial components are

```math
u(x,y,0)
=
\sqrt{0.5+\delta_z}\,
\phi_b(x,y)
e^{0.12i\sin y},
```

```math
d(x,y,0)
=
\sqrt{0.5-\delta_z}\,
\phi_b(x,y)
e^{i\theta_{\mathrm{rel}}(x,y)}.
```

The complete initial spinor is

```math
\Psi(x,y,0)
=
\begin{pmatrix}
u(x,y,0)\\
d(x,y,0)
\end{pmatrix},
```

followed by vortex imprinting and mean-density normalization.

---

## 4. Reduced bond and vibrational dynamics

### 4.1 Morse potential

```math
V_M(R)
=
D_e
\left[
1-e^{-a(R-R_e)}
\right]^2
-D_e.
```

At equilibrium,

```math
V_M(R_e)=-D_e,
```

while

```math
\lim_{R\rightarrow\infty}V_M(R)=0.
```

### 4.2 Morse force

```math
F_M(R)
=
-\frac{dV_M}{dR},
```

or explicitly,

```math
F_M(R)
=
-2D_ea
e^{-a(R-R_e)}
\left[
1-e^{-a(R-R_e)}
\right].
```

### 4.3 Approximate Morse vibrational levels

The uploaded source retained only the anharmonic correction term. Written as the usual reduced Morse-level expression,

```math
E_v
\approx
\hbar\omega_e\left(v+\frac12\right)
-
\frac{
\left[
\hbar\omega_e\left(v+\frac12\right)
\right]^2
}{
4D_e
}.
```

> Verify this line against the implementation if the code stores only the anharmonic correction rather than the full approximate level.

### 4.4 Reduced vibrational energy

```math
E_{\mathrm{vib}}
=
V_M(R)
+
\frac12\mu\dot R^2\kappa_R,
```

where the code uses the reduced scaling

```math
\kappa_R=0.025.
```

### 4.5 Density-weighted optical pressure

```math
F_{\mathrm{opt}}
=
\operatorname{mean}
\left[
V_{\mathrm{laser}}
\frac{x}{L/2+\varepsilon}
\rho
\right].
```

### 4.6 Total reduced bond force

The source retained only the damping fragment. The surrounding definitions imply the reduced force form

```math
F_R
=
F_M+F_{\mathrm{opt}}-\gamma_R\dot R.
```

The corresponding acceleration is

```math
a_R
=
\frac{F_R}{\max(\mu,0.2)}.
```

### 4.7 Velocity-Verlet update

```math
\dot R_{n+1/2}
=
\dot R_n+\frac{\Delta t}{2}a_n,
```

```math
R_{n+1}
=
R_n+\Delta t\,\dot R_{n+1/2},
```

```math
\dot R_{n+1}
=
\dot R_{n+1/2}
+
\frac{\Delta t}{2}a_{n+1}.
```

The implementation clips $R$ to a specified interval and applies a reduced rebound factor when a boundary is reached.

The Morse potential, force, levels, bond observables, and reduced velocity-Verlet coupling are implemented explicitly.

---

## 5. Structured light and vortices

### 5.1 Hermite-Gaussian mode

```math
E_{mn}^{\mathrm{HG}}(x,y,t)
=
H_m
\left(
\frac{\sqrt2x}{w}
\right)
H_n
\left(
\frac{\sqrt2y}{w}
\right)
e^{-(x^2+y^2)/w^2}
e^{-i(\omega t+\phi_0)}.
```

This is a transverse waist-plane profile rather than a complete propagating Maxwell solution.

### 5.2 Laguerre-Gaussian mode

Let

```math
r=\sqrt{x^2+y^2},
\qquad
\theta=\operatorname{atan2}(y,x),
```

and

```math
\rho_L=\frac{2r^2}{w^2}.
```

Then

```math
E_{p\ell}^{\mathrm{LG}}(x,y,t)
=
\left(
\frac{\sqrt2r}{w}
\right)^{|\ell|}
L_p^{|\ell|}(\rho_L)
e^{-r^2/w^2}
e^{i(\ell\theta-\omega t+\phi_0)}.
```

The phase factor

```math
e^{i\ell\theta}
```

provides the programmed optical-vortex winding.

### 5.3 Optical coupling channels

For mixed HG/LG operation, the reduced scalar potential is

```math
V_{\mathrm{light}}
=
A_Lc_s
\left[
\operatorname{Re}(E^{\mathrm{HG}})
+
\operatorname{Re}(E^{\mathrm{LG}})
\right].
```

The reduced spin-mixing drive is

```math
s_{\mathrm{light}}
=
A_Lc_{\mathrm{spin}}
\left[
\operatorname{Re}(E^{\mathrm{HG}})
+
\operatorname{Im}(E^{\mathrm{LG}})
\right].
```

These are phenomenological coupling channels rather than transition-dipole or Maxwell-derived interaction Hamiltonians.

### 5.4 Static vortex coordinates

For vortex $j$,

```math
r_j^2
=
(x-x_j)^2+(y-y_j)^2,
```

```math
\theta_j
=
\operatorname{atan2}(y-y_j,x-x_j).
```

The accumulated phase is

```math
\Phi_v
=
\sum_jq_j\theta_j.
```

The depleted amplitude envelope is

```math
A_v
=
\prod_j
\operatorname{clip}
\left[
1-d_j
e^{-r_j^2/(2\xi_j^2)},
0.04,
1
\right].
```

The full implementation imprints

```math
u
\longrightarrow
A_v e^{i\Phi_v}u,
```

```math
d
\longrightarrow
A_v e^{-0.55i\Phi_v}d.
```

The second component therefore receives a project-specific phase-defect texture rather than a conventional integer-winding vortex.

### 5.5 Moving-vortex drives

Define

```math
G_j
=
e^{-r_j^2/(2\xi_j^2)}.
```

The scalar drive is

```math
V_{\mathrm{vortex}}
=
\kappa_v
\sum_jd_jG_j,
```

and the spin drive is

```math
s_{\mathrm{vortex}}
=
\kappa_v
\sum_j
d_jq_jG_j\sin\theta_j.
```

The current source explicitly labels these moving-vortex terms as phenomenological perturbations.

---

## 6. Nonlinear spinor-field evolution

### 6.1 Total potential and spin drive

```math
V_{\mathrm{tot}}
=
V_{\mathrm{mol}}
+
V_{\mathrm{light}}
+
V_{\mathrm{vortex}}
+
V_{\mathrm{vertex}},
```

```math
\Omega
=
\Omega_0
+
s_{\mathrm{light}}
+
s_{\mathrm{vortex}}
+
s_{\mathrm{vertex}}.
```

Here $\Omega_0$ is the parameter named `zeeman` in the code. It multiplies $\sigma_x$, so it acts as a transverse mixing term rather than a $\sigma_z$ energy splitting.

### 6.2 Compact Hamiltonian

The uploaded source retained the spin-orbit term but not the complete scalar/nonlinear operator. A syntax-safe decomposition is

```math
H[\Psi,t]
=
H_0[\Psi,t]
+
\Omega\sigma_x
+
i\lambda
\left(
\sigma_xD_x+\sigma_yD_y
\right),
```

where $H_0[\Psi,t]$ denotes the implemented scalar kinetic, potential, and nonlinear contributions.

### 6.3 Density-relaxation term

```math
\Gamma(\rho)
=
\gamma
\frac{
\rho_{\mathrm{target}}-\rho
}{
1+\rho
}.
```

### 6.4 Evolution equation

```math
\frac{\partial\Psi}{\partial t}
=
-iH[\Psi,t]\Psi
+
\Gamma(\rho)\Psi.
```

Equivalently,

```math
i\frac{\partial\Psi}{\partial t}
=
H[\Psi,t]\Psi
+
i\Gamma(\rho)\Psi.
```

### 6.5 Component form

Writing the scalar part of $H_0$ acting on each component as $H_s$, the upper component has the source-preserved derivative structure

```math
i\frac{\partial u}{\partial t}
=
H_su
+
\Omega d
+
i\lambda D_xd
+
\lambda D_yd
+
i\Gamma u.
```

The lower component is recorded in the source as

```math
i\frac{\partial d}{\partial t}
=
H_sd
+
\Omega u
+
i\lambda D_xu
+
\lambda D_yu
+
i\Gamma d.
```

> **Implementation check:** verify the sign of the lower-component $D_yu$ term against the actual $\sigma_y$ convention used in code.

The component equations follow directly from the implemented scalar, spin-mixing, derivative-coupling, and gain terms.

---

## 7. Finite differences and Strang propagation

### 7.1 Periodic Laplacian

```math
(\nabla_h^2f)_{ij}
=
\frac{
f_{i+1,j}
+
f_{i-1,j}
+
f_{i,j+1}
+
f_{i,j-1}
-
4f_{ij}
}{
(\Delta x)^2
}.
```

### 7.2 Centered derivatives

```math
(D_xf)_{ij}
=
\frac{
f_{i+1,j}-f_{i-1,j}
}{
2\Delta x
},
```

```math
(D_yf)_{ij}
=
\frac{
f_{i,j+1}-f_{i,j-1}
}{
2\Delta x
}.
```

### 7.3 Exact finite-difference Fourier symbols

```math
\xi_x
=
\frac{\sin(k_x\Delta x)}{\Delta x},
```

```math
\xi_y
=
\frac{\sin(k_y\Delta x)}{\Delta x},
```

```math
\epsilon_h(\mathbf k)
=
\frac{4\kappa}{(\Delta x)^2}
\left[
\sin^2
\left(
\frac{k_x\Delta x}{2}
\right)
+
\sin^2
\left(
\frac{k_y\Delta x}{2}
\right)
\right].
```

These are the symbols of the existing finite-difference operators, not a replacement continuous $k^2$ Laplacian.

### 7.4 Kinetic and SOC Fourier operator

```math
H_K(\mathbf k)
=
\epsilon_h(\mathbf k)I
+
\lambda
\left(
\xi_x\sigma_x+\xi_y\sigma_y
\right).
```

Define

```math
r_{\mathrm{soc}}
=
\lambda
\sqrt{
\xi_x^2+\xi_y^2
}.
```

The exact Fourier-space substep over duration $\tau$ is

```math
U_K(\tau,\mathbf k)
=
e^{-i\tau\epsilon_h(\mathbf k)}
\left[
\cos(\tau r_{\mathrm{soc}})I
-
i
\frac{
\sin(\tau r_{\mathrm{soc}})
}{
r_{\mathrm{soc}}
}
\lambda
\left(
\xi_x\sigma_x+\xi_y\sigma_y
\right)
\right].
```

At $r_{\mathrm{soc}}=0$, the ratio is evaluated by its continuous limit.

### 7.5 Exact local substep

Using $V_{\mathrm{eff}}$ for the local scalar/nonlinear phase term,

```math
U_{\mathrm{local}}(\Delta t)
=
e^{-i\Delta t V_{\mathrm{eff}}}
\left[
\cos(\Delta t\Omega)I
-
i\sin(\Delta t\Omega)\sigma_x
\right].
```

### 7.6 Second-order Strang composition

```math
U(\Delta t)
\approx
U_K
\left(
\frac{\Delta t}{2}
\right)
U_{\mathrm{local}}(\Delta t)
U_K
\left(
\frac{\Delta t}{2}
\right).
```

With zero damping and no repeated mass projection, each subflow is unitary up to floating-point error.

### 7.7 LG phase-kick rate

```math
\theta_{\mathrm{kick}}
=
\delta t\,
k_{\phi}
A_L
E_{\mathrm{LG}}^{\mathrm{real}}.
```

The components are updated by

```math
u
\longrightarrow
e^{i\theta_{\mathrm{kick}}}u,
```

```math
d
\longrightarrow
e^{-i\theta_{\mathrm{kick}}/2}d.
```

---

## 8. Density, bonding, information, and geometric diagnostics

### 8.1 Density and collective phase

```math
\rho
=
|u|^2+|d|^2,
```

```math
\phi
=
\arg(u+d).
```

### 8.2 Phase-modulated density

```math
\rho_{\mathrm{mod}}
=
\rho
\left[
1+0.28\cos\phi
\right].
```

This is a nonnegative phase-modulated density display.

### 8.3 Signed phase quadrature

```math
Q_{\rho}
=
\rho\cos\phi.
```

Negative $Q_{\rho}$ means phase inversion. It does not mean negative electron density, and it is not a Wigner function.

### 8.4 Bonding projection

```math
n_b
=
\sum_{s\in\{\uparrow,\downarrow\}}
\left|
\langle\phi_b,\psi_s\rangle_h
\right|^2.
```

### 8.5 Antibonding projection

```math
n_{ab}
=
\sum_{s\in\{\uparrow,\downarrow\}}
\left|
\langle\phi_{ab},\psi_s\rangle_h
\right|^2.
```

### 8.6 Transition score

```math
T
=
\frac{
n_{ab}
}{
n_b+n_{ab}+\varepsilon
}.
```

### 8.7 Effective bond-order proxy

```math
B
=
\frac{
n_b-n_{ab}
}{
n_b+n_{ab}+\varepsilon
}.
```

### 8.8 Two-site sharing proxy

Define

```math
n_A
=
\sum_s
\left|
\langle\phi_A,\psi_s\rangle_h
\right|^2,
```

```math
n_B
=
\sum_s
\left|
\langle\phi_B,\psi_s\rangle_h
\right|^2,
```

and

```math
p_A
=
\frac{n_A}{n_A+n_B+\varepsilon},
\qquad
p_B
=
\frac{n_B}{n_A+n_B+\varepsilon}.
```

The sharing proxy is

```math
C_{\mathrm{share}}
=
\operatorname{clip}
\left[
4p_Ap_B\max(S,0),
0,
1
\right].
```

These are projection and sharing diagnostics, not ab-initio bond orders or orbital entanglement measurements.

---

## 9. Bloch geometry and topology

### 9.1 Local Bloch vector

```math
n_x
=
\frac{
2\operatorname{Re}(u^{*}d)
}{
\rho+\varepsilon
},
```

```math
n_y
=
\frac{
2\operatorname{Im}(u^{*}d)
}{
\rho+\varepsilon
},
```

```math
n_z
=
\frac{
|u|^2-|d|^2
}{
\rho+\varepsilon
}.
```

Collectively,

```math
\mathbf n
=
\frac{
\Psi^{\dagger}
\boldsymbol{\sigma}
\Psi
}{
\Psi^{\dagger}\Psi+\varepsilon
}.
```

### 9.2 Spin-texture curvature diagnostic

```math
K_{\mathrm{spin}}
=
|\nabla\mathbf n|^2,
```

or discretely,

```math
K_{\mathrm{spin}}
=
\sum_{a=x,y,z}
\left[
(D_xn_a)^2+(D_yn_a)^2
\right].
```

This is spin-texture gradient energy, not spacetime curvature.

### 9.3 Topological-charge density

```math
q(x,y,t)
=
\frac{1}{4\pi}
\mathbf n
\cdot
\left(
\partial_x\mathbf n
\times
\partial_y\mathbf n
\right).
```

### 9.4 Integrated charge

```math
Q(t)
=
\int_{\Omega}
q(x,y,t)\,dA.
```

On the masked grid,

```math
Q_h(t)
=
\sum_{i,j}
m_{ij}
q_{ij}
(\Delta x)^2.
```

Because low-density points and derivative neighbors are masked, the computed value need not be exactly integer.

---

## 10. QFI, entropy, memory, OAM, and winding

### 10.1 Spin-$z$ expectation

```math
\langle\sigma_z\rangle
=
\frac{
\int
\left(
|u|^2-|d|^2
\right)dA
}{
\int\rho\,dA
}.
```

### 10.2 QFI diagnostic

Since

```math
\sigma_z^2=I,
```

the implemented pure-state-style diagnostic is

```math
F_Q^{(z)}
=
4
\left[
1-\langle\sigma_z\rangle^2
\right].
```

### 10.3 Spatial entropy diagnostic

Define

```math
p_{ij}
=
\frac{
\rho_{ij}
}{
\sum_{k,l}\rho_{kl}+\varepsilon
}.
```

Then

```math
S_{\mathrm{spatial}}
=
-\sum_{i,j}
p_{ij}\ln(p_{ij}+\varepsilon).
```

The QFI and entropy formulas are implemented as relative diagnostics.

### 10.4 BLP-inspired density distinguishability

For two nearby spinor trajectories $a$ and $b$,

```math
D_{\rho}(t)
=
\frac{
\frac12
\int
|\rho_a-\rho_b|\,dA
}{
\frac12
\int
(\rho_a+\rho_b)\,dA
+
\varepsilon
}.
```

### 10.5 Memory current

```math
J_{\mathrm{mem}}(t_n)
=
\frac{
D_{\rho}(t_n)-D_{\rho}(t_{n-1})
}{
\Delta t
}.
```

A positive value indicates revival of density distinguishability in this reduced diagnostic. It is BLP-inspired, but it is not the canonical trace distance between optimized reduced density matrices.

### 10.6 OAM expectation

```math
\langle L_z\rangle
=
\frac{
\sum_s
\operatorname{Im}
\left[
\int
\psi_s^{*}
\left(
x\partial_y-y\partial_x
\right)
\psi_s\,dA
\right]
}{
\int\rho\,dA
}.
```

This quantity is expressed in reduced units.

### 10.7 Boundary phase winding

```math
\nu_{\partial\Omega}
=
\frac{1}{2\pi}
\oint_{\partial\Omega}
d\phi.
```

The code evaluates the wrapped phase differences around the periodic grid boundary.

---

## 11. Toy MEAO/QIT reduced-state diagnostics

The local Fock-space labels are

```math
\left\{
0,\uparrow,\downarrow,\uparrow\downarrow
\right\}
```

for each of two orbital subsystems.

These formulas describe the live reduced proxy only. They are not a validated MEAO orbital optimization and are not time-resolved ab-initio entanglement values.

### 11.1 Covalent reference state

The uploaded source omitted the opening coefficient and at least part of this state definition. The retained fragment is

```math
|\psi_{\mathrm{cov}}\rangle
\propto
|\downarrow,\uparrow\rangle
+
|\uparrow\downarrow,0\rangle.
```

> Replace this partial expression with the exact implemented state before treating it as a normative definition.

### 11.2 Antibonding-contaminated reference

The retained source fragment is

```math
|\psi_{\mathrm{ab}}\rangle
\propto
|\uparrow,\downarrow\rangle
+
|\downarrow,\uparrow\rangle
+
|\uparrow\downarrow,0\rangle.
```

> Replace the proportionality with the exact normalization and verify whether additional basis terms appear in the implementation.

### 11.3 Ionic reference

For molecular polarity $p$,

```math
|\psi_{\mathrm{ionic}}\rangle
=
\mathcal N_I
\left[
\sqrt{\frac{1+p}{2}}
|0,\uparrow\downarrow\rangle
+
\sqrt{\frac{1-p}{2}}
|\uparrow\downarrow,0\rangle
\right].
```

### 11.4 Ionic mixing coefficient

```math
\iota
=
\operatorname{clip}
\left[
0.12|p|
+
0.22T
+
0.12\max(0,1-B),
0,
0.90
\right].
```

### 11.5 Toy two-orbital state

Using

```math
[x]_+=\max(0,x),
```

the proxy state is

```math
|\psi_{\mathrm{toy}}\rangle
=
\mathcal N
\left[
\sqrt{
[1-\iota-0.15T]_+
}
|\psi_{\mathrm{cov}}\rangle
+
\sqrt{\iota}
|\psi_{\mathrm{ionic}}\rangle
+
\sqrt{0.15T}
|\psi_{\mathrm{ab}}\rangle
\right].
```

### 11.6 Reduced density operator

```math
\rho_{AB}
=
|\psi_{\mathrm{toy}}\rangle
\langle\psi_{\mathrm{toy}}|.
```

### 11.7 Partial traces

```math
\rho_A
=
\operatorname{Tr}_B(\rho_{AB}),
```

```math
\rho_B
=
\operatorname{Tr}_A(\rho_{AB}).
```

### 11.8 Von Neumann entropy

```math
S(\rho)
=
-\operatorname{Tr}
\left[
\rho\ln(\rho+\varepsilon)
\right].
```

### 11.9 Mutual information

The source retained only the final $S(\rho_{AB})$ fragment. The complete standard reduced-state expression is

```math
I(A:B)
=
S(\rho_A)
+
S(\rho_B)
-
S(\rho_{AB}).
```

The normalized proxy is

```math
I_{\mathrm{norm}}
=
\frac{
I(A:B)
}{
2\ln4
}.
```

### 11.10 Logarithmic negativity

```math
E_N
=
\max
\left[
0,
\log_2
\left(
\left\|
\rho_{AB}^{T_B}
\right\|_1+\varepsilon
\right)
\right].
```

### 11.11 Four-level I-concurrence-style diagnostic

```math
C_I
=
\frac{
\sqrt{
2
\left[
1-\operatorname{Tr}(\rho_A^2)
\right]
}
}{
\sqrt{
2(1-1/4)
}
}.
```

This is not Wootters' two-qubit concurrence.

### 11.12 MEAO coherence proxy

```math
F_{\mathrm{MEAO}}^{\mathrm{proxy}}
=
\left|
\rho_{(0,\uparrow\downarrow),
(\uparrow\downarrow,0)}
\right|^2
+
\left|
\rho_{(\uparrow,\downarrow),
(\downarrow,\uparrow)}
\right|^2.
```

The normalized form is

```math
F_{\mathrm{norm}}
=
\operatorname{clip}
\left[
\frac{
F_{\mathrm{MEAO}}^{\mathrm{proxy}}
}{
0.125
},
0,
1.5
\right].
```

### 11.13 Fractional bond-order proxy

```math
m_{\mathrm{frac}}
=
m_0
\left[
0.50\max(0,B)
+
0.30I_{\mathrm{norm}}
+
0.20\min(F_{\mathrm{norm}},1)
\right]
(1-0.18T),
```

followed by scenario-dependent clipping.

The complete state construction and the resulting mutual-information, negativity, concurrence, coherence, and fractional-bond formulas are implemented in the toy MEAO/QIT analyzer.

---

## 12. Clifford-algebra encoding

### 12.1 Clifford relation

The real Clifford generators satisfy

```math
e_ie_j+e_je_i
=
2g_{ij}.
```

For $\mathrm{Cl}(3,0)$,

```math
g
=
\operatorname{diag}(1,1,1).
```

For $\mathrm{Cl}(1,3)$,

```math
g
=
\operatorname{diag}(1,-1,-1,-1).
```

The compact engine represents $n=p+q$ generators with

```math
\dim\mathrm{Cl}(p,q)=2^n.
```

It implements genuine geometric products and reversion.

### 12.2 Density-weighted global spinor

Define

```math
w_{ij}
=
\frac{
\rho_{ij}
}{
\sum_{k,l}\rho_{kl}+\varepsilon
}.
```

Then

```math
\alpha_{\mathrm{raw}}
=
\sum_{i,j}
w_{ij}u_{ij},
```

```math
\beta_{\mathrm{raw}}
=
\sum_{i,j}
w_{ij}d_{ij}.
```

Normalize:

```math
\chi
=
\frac{
1
}{
\sqrt{
|\alpha_{\mathrm{raw}}|^2
+
|\beta_{\mathrm{raw}}|^2
+\varepsilon
}
}
\begin{pmatrix}
\alpha_{\mathrm{raw}}\\
\beta_{\mathrm{raw}}
\end{pmatrix}.
```

### 12.3 Complex coefficients

```math
\alpha=a+ib,
\qquad
\beta=c+id.
```

### 12.4 Even subalgebra

```math
\mathrm{Cl}^{+}(3,0)
=
\operatorname{span}
\left\{
1,e_{12},e_{23},e_{31}
\right\}.
```

This algebra is quaternion-like.

### 12.5 Spinor-to-Clifford map

```math
\Phi
=
a
+
b e_{12}
+
c e_{23}
+
d e_{31}.
```

### 12.6 Primitive idempotent

```math
f
=
\frac{1+e_3}{2},
\qquad
f^2=f.
```

### 12.7 Minimal-left-ideal representative

```math
\Phi_L
=
\Phi f.
```

### 12.8 Clifford density element

```math
\rho_c
=
\Phi_L
\widetilde{\Phi_L}.
```

For a grade-$k$ component,

```math
\widetilde{A_k}
=
(-1)^{k(k-1)/2}A_k.
```

### 12.9 Rotor-angle diagnostic

```math
\theta_R
=
2\,\operatorname{atan2}
\left(
|\langle\Phi\rangle_2|,
|\langle\Phi\rangle_0|+\varepsilon
\right).
```

### 12.10 Purity-like Clifford diagnostic

Let

```math
s_c=\langle\rho_c\rangle_0,
```

```math
v_c=|\langle\rho_c\rangle_1|.
```

Then

```math
P_c
=
\frac{
s_c^2+v_c^2
}{
|\rho_c|^2+\varepsilon
}.
```

The source implements this encoding as a representation and descriptor system, not as the PDE time-evolution operator.

---

## 13. Feature chart, PCA, and ridge-DMD

At each stored time, the project forms a diagnostic vector

```math
\mathbf x_t
=
\begin{pmatrix}
\text{memory current}\\
\text{QFI}\\
\text{curvature}\\
\text{OAM}\\
\text{phase winding}\\
\text{bond order}\\
\text{MEAO/QIT proxies}\\
\text{bond coordinate}\\
\text{Clifford descriptors}\\
\text{laser controls}\\
\text{local patch descriptors}
\end{pmatrix}.
```

The actual feature chart includes both global observables and compact spatial Clifford summaries.

### 13.1 Training-only standardization

```math
\boldsymbol{\mu}
=
\frac{1}{N_{\mathrm{train}}}
\sum_{t\in\mathrm{train}}
\mathbf x_t,
```

```math
\mathbf s
=
\operatorname{std}_{t\in\mathrm{train}}
(\mathbf x_t),
```

```math
\widehat{\mathbf x}_t
=
\frac{
\mathbf x_t-\boldsymbol{\mu}
}{
\mathbf s
}.
```

The division is componentwise.

### 13.2 Training-block SVD

For the standardized training matrix $\widehat X_{\mathrm{train}}$,

```math
\widehat X_{\mathrm{train}}
=
U\Sigma V^{T}.
```

Let $W$ contain the first $k=3$ rows of $V^{T}$.

### 13.3 PCA latent coordinates

```math
\mathbf z_t
=
\widehat{\mathbf x}_tW^{T}.
```

These are learned diagnostic coordinates. They are not gamma-matrix axes and are not an intrinsic manifold reconstruction.

### 13.4 Ridge-DMD transition matrix

Construct

```math
X
=
\begin{pmatrix}
\mathbf z_0\\
\mathbf z_1\\
\vdots\\
\mathbf z_{N_{\mathrm{train}}-2}
\end{pmatrix},
```

and

```math
Y
=
\begin{pmatrix}
\mathbf z_1\\
\mathbf z_2\\
\vdots\\
\mathbf z_{N_{\mathrm{train}}-1}
\end{pmatrix}.
```

Then

```math
A_{\mathrm{DMD}}
=
\left(
X^{T}X+\lambda_RI
\right)^{-1}
X^{T}Y.
```

### 13.5 One-step prediction

```math
\widehat{\mathbf z}_{t+1}
=
\mathbf z_tA_{\mathrm{DMD}}.
```

### 13.6 One-step error

```math
e_{t+1}
=
\left\|
\mathbf z_{t+1}
-
\widehat{\mathbf z}_{t+1}
\right\|_2.
```

### 13.7 Recursive rollout

```math
\widehat{\mathbf z}_{t+h}
=
\mathbf z_t
A_{\mathrm{DMD}}^h.
```

### 13.8 Skill against persistence

```math
\operatorname{Skill}
=
1-
\frac{
\operatorname{RMSE}_{\mathrm{DMD}}
}{
\operatorname{RMSE}_{\mathrm{persistence}}
}.
```

The scaler, PCA basis, and DMD operator are fitted only on an early contiguous block and evaluated on later held-out states.

---

## 14. Hopf geometry and $S^3$ state-space fibers

### 14.1 Normalized complex spinor on $S^3$

Let

```math
z
=
\begin{pmatrix}
z_1\\
z_2
\end{pmatrix},
```

with

```math
|z_1|^2+|z_2|^2=1.
```

The four real coordinates are

```math
q
=
\begin{pmatrix}
\operatorname{Re}z_1\\
\operatorname{Im}z_1\\
\operatorname{Re}z_2\\
\operatorname{Im}z_2
\end{pmatrix}
\in S^3\subset\mathbb R^4.
```

### 14.2 Hopf map

For a possibly unnormalized pair,

```math
\mathbf n
=
\frac{1}{
|z_1|^2+|z_2|^2
}
\begin{pmatrix}
2\operatorname{Re}(z_1^{*}z_2)\\
2\operatorname{Im}(z_1^{*}z_2)\\
|z_1|^2-|z_2|^2
\end{pmatrix}.
```

For a normalized spinor, the denominator is one.

### 14.3 Local chart phase

The source also records the chart-local phase

```math
\chi_{\mathrm{local}}
=
\frac12\arg(ud).
```

This is not a globally defined fiber coordinate at component zeros.

### 14.4 Canonical lift from $S^2$

For

```math
\mathbf n
=
\begin{pmatrix}
\sin\theta\cos\varphi\\
\sin\theta\sin\varphi\\
\cos\theta
\end{pmatrix},
```

one canonical lift is

```math
z_1
=
\cos\frac{\theta}{2},
```

```math
z_2
=
e^{i\varphi}
\sin\frac{\theta}{2}.
```

### 14.5 Complete $U(1)$ fiber

Every common phase maps to the same Bloch point:

```math
(z_1,z_2)
\longrightarrow
e^{i\alpha}(z_1,z_2).
```

The complete fiber is

```math
\mathcal F_z
=
\left\{
e^{i\alpha}(z_1,z_2)
:
0\leq\alpha<2\pi
\right\}.
```

In real $\mathbb R^4$ coordinates,

```math
q(\alpha)
=
\mathbf u\cos\alpha+\mathbf v\sin\alpha.
```

### 14.6 Clifford torus coordinates

```math
z_1
=
\cos\eta\,e^{i\xi_1},
```

```math
z_2
=
\sin\eta\,e^{i\xi_2}.
```

Therefore,

```math
q
=
\begin{pmatrix}
\cos\eta\cos\xi_1\\
\cos\eta\sin\xi_1\\
\sin\eta\cos\xi_2\\
\sin\eta\sin\xi_2
\end{pmatrix}.
```

The corresponding Bloch latitude satisfies

```math
n_z=\cos(2\eta).
```

### 14.7 Stereographic projection

After a fixed $\mathbb R^4$ rotation $A$, let

```math
q'=Aq.
```

The raw $\mathbb R^3$ stereographic coordinate is

```math
\mathbf p
=
\frac{
(q'_1,q'_2,q'_3)
}{
1-q'_4
}.
```

### 14.8 Finite display compaction

The finite display coordinate is

```math
\mathbf p_{\mathrm{display}}
=
\frac{
R_c\mathbf p
}{
R_c+|\mathbf p|
}.
```

These are complete $U(1)$ state-space fibers projected for visualization. They are not physical-space three-dimensional Hopfion solitons.

---

## 15. Future braid-resolved information experiments

This section is proposed future work. It is not a current measured entanglement observable.

### 15.1 Braid word

For resolved temporal crossings,

```math
\mathcal B
=
\sigma_{i_1}^{\epsilon_1}
\sigma_{i_2}^{\epsilon_2}
\cdots
\sigma_{i_m}^{\epsilon_m},
```

where

```math
\epsilon_k\in\{-1,+1\}
```

records the orientation of each crossing.

### 15.2 Information change over lag $\tau$

```math
\Delta_{\tau}I(t)
=
I(t+\tau)-I(t).
```

### 15.3 Braid-feature prediction model

```math
\Delta_{\tau}I(t)
=
\beta_0
+
\sum_{j=1}^{p}
\beta_jT_j(t)
+
\varepsilon(t+\tau).
```

Here $T_j(t)$ may include:

- crossing rate;
- braid-word complexity;
- braid entropy;
- phase holonomy;
- writhe-like descriptors;
- motif persistence.

### 15.4 Controlled prediction model

To separate braiding from shared external drivers,

```math
\Delta_{\tau}I(t)
=
\beta_0
+
\boldsymbol{\beta}^{T}\mathbf T(t)
+
\boldsymbol{\gamma}^{T}\mathbf C(t)
+
\varepsilon(t+\tau),
```

where $\mathbf C(t)$ includes controls such as:

- molecule identity;
- current $I(t)$;
- laser amplitude;
- LG index $\ell$;
- vortex radius and depth;
- bond coordinate $R(t)$;
- curvature;
- dissipative or conservative regime.

The key test is whether

```math
\text{controls + braid features}
```

predict held-out information changes better than

```math
\text{controls only}.
```

### 15.5 Future normal-mode extension

For a nonlinear $N$-atom molecule,

```math
q_k(t),
\qquad
k=1,\ldots,3N-6.
```

For a linear molecule,

```math
q_k(t),
\qquad
k=1,\ldots,3N-5.
```

This would replace the present single bond-stretch coordinate $R(t)$ with a vector of normal-mode coordinates.

---

## Scientific interpretation boundary

The implemented hierarchy is

```math
\text{PDE state}
\rightarrow
\text{derived diagnostics}
\rightarrow
\text{Clifford/PCA representation}
\rightarrow
\text{Hopf visualization}.
```

It does not currently establish

```math
\text{molecular entanglement}.
```

Any future braid relationship must be reported separately as

```math
\text{correlation},
```

```math
\text{held-out prediction},
```

or

```math
\text{causal evidence}.
```
