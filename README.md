MEAO Spinor Lab: Core Equation Spine
Scope

The project is organized mathematically as

[
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
]

Sections 1 through 9 describe equations implemented in the current reduced-order framework. Section 10 contains proposed future braid experiments and is not yet an implemented entanglement measurement.

1. Grid, inner product, and normalization
1.1 Spatial grid

For a square domain of side length (L) containing (N\times N) grid points,

[
\Delta x=\frac{L}{N},
]

[
x_i=-\frac{L}{2}+i\Delta x,
\qquad
i=0,\ldots,N-1.
]

The grid is periodic because the finite-difference operators use wrapped neighbors.

1.2 Discrete inner product

For two scalar fields (a) and (b),

\sum_{i,j}
a_{ij}^{*}b_{ij}
(\Delta x)^2.
]

1.3 Scalar-orbital normalization

A scalar orbital is normalized by

[
\phi
\longrightarrow
\frac{\phi}
{
\sqrt{
\sum_{i,j}|\phi_{ij}|^2(\Delta x)^2+\varepsilon
}
}.
]

Thus,




]

1.4 Two-component spinor

The simulated field is

\begin{pmatrix}
\psi_{\uparrow}(x,y,t)\
\psi_{\downarrow}(x,y,t)
\end{pmatrix}.
]

Its density is

|u|^2+|d|^2.
]

The integrated field mass is

\sum_{i,j}
\rho_{ij}(t)
(\Delta x)^2.
]

1.5 Mean-density normalization

The current full model sets the initial target mass to

\rho_{\mathrm{target}}A,
\qquad
A=L^2.
]

The rescaling is

[
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
]

Equivalently,

\rho_{\mathrm{target}}.
]

This is a mean-field amplitude convention, not a one-electron probability normalization. The conservative launch disables repeated normalization after every time step.

These definitions are implemented directly in the grid, inner-product, scalar-normalization, and spinor-normalization functions.

2. Analytic H₂/CO molecular scaffold
2.1 Atomic-center locations

For bond coordinate (R),

[
x_A=-\frac{R}{2},
\qquad
x_B=+\frac{R}{2}.
]

2.2 Softened radial coordinates

For atom (A),

\sqrt{
(x-x_A)^2+(y-y_A)^2+s_{\mathrm{AO}}^2
}.
]

Likewise,

\sqrt{
(x-x_B)^2+(y-y_B)^2+s_{\mathrm{AO}}^2
}.
]

The implementation uses

[
s_{\mathrm{AO}}=0.08.
]

2.3 Slater-like atomic envelopes

\mathcal N_A e^{-\zeta_A r_A},
]

\mathcal N_B e^{-\zeta_B r_B},
]

where (\mathcal N_A) and (\mathcal N_B) are numerical normalization factors.

2.4 Atomic-orbital overlap

\mathrm{Re}
\langle\phi_A,\phi_B\rangle_h.
]

2.5 Bonding channel

c_A\phi_A+c_B\phi_B,
]

\frac{
\phi_b^{\mathrm{raw}}
}{
\sqrt{
\langle\phi_b^{\mathrm{raw}},
\phi_b^{\mathrm{raw}}\rangle_h
+\varepsilon
}
}.
]

2.6 Antibonding channel

c_B\phi_A-c_A\phi_B,
]

\frac{
\phi_{ab}^{\mathrm{raw}}
}{
\sqrt{
\langle\phi_{ab}^{\mathrm{raw}},
\phi_{ab}^{\mathrm{raw}}\rangle_h
+\varepsilon
}
}.
]

2.7 Reduced molecular potential

Define softened potential distances

\sqrt{
(x-x_A)^2+y^2+s_V^2
},
]

\sqrt{
(x-x_B)^2+y^2+s_V^2
},
]

with

[
s_V=0.18.
]

The reduced two-center potential is

-V_0
\left[
\frac{Z_A}{Z_A+Z_B}\frac{1}{r_A^{(V)}}
+
\frac{Z_B}{Z_A+Z_B}\frac{1}{r_B^{(V)}}
\right]
+
0.045,p,x,
]

where (p) is the heteronuclear polarity parameter. For H₂, (p=0); for CO, the nonzero polarity term breaks left-right symmetry.

These molecular envelopes, overlaps, channels, and reduced potentials are implemented as analytic two-center approximations.

3. Initial spinor texture

Define the initial spin imbalance

[
\delta_z=0.10p,
]

and relative phase

0.35\tanh y+0.15px.
]

The initial components are

\sqrt{0.5+\delta_z},
\phi_b(x,y)
e^{0.12i\sin y},
]

\sqrt{0.5-\delta_z},
\phi_b(x,y)
e^{i\theta_{\mathrm{rel}}(x,y)}.
]

The complete initial spinor is

\begin{pmatrix}
u(x,y,0)\
d(x,y,0)
\end{pmatrix},
]

followed by vortex imprinting and mean-density normalization.

4. Reduced bond and vibrational dynamics
4.1 Morse potential

D_e
\left[
1-e^{-a(R-R_e)}
\right]^2
-D_e.
]

At equilibrium,

[
V_M(R_e)=-D_e,
]

while

[
\lim_{R\rightarrow\infty}V_M(R)=0.
]

4.2 Morse force

-\frac{dV_M}{dR},
]

or explicitly,

-2D_ea
e^{-a(R-R_e)}
\left[
1-e^{-a(R-R_e)}
\right].
]

4.3 Approximate Morse vibrational levels

\frac{
\left[
\hbar\omega_e(v+\frac12)
\right]^2
}{
4D_e
}.
]

4.4 Reduced vibrational energy

V_M(R)
+
\frac12\mu\dot R^2\kappa_R,
]

where the code uses the reduced scaling

[
\kappa_R=0.025.
]

4.5 Density-weighted optical pressure

\mathrm{mean}
\left[
V_{\mathrm{laser}}
\frac{x}{L/2+\varepsilon}
\rho
\right].
]

4.6 Total reduced bond force

\gamma_R\dot R.
]

The corresponding acceleration is

\frac{F_R}{\max(\mu,0.2)}.
]

4.7 Velocity-Verlet update

\dot R_n+\frac{\Delta t}{2}a_n,
]

R_n+\Delta t,\dot R_{n+1/2},
]

\dot R_{n+1/2}
+
\frac{\Delta t}{2}a_{n+1}.
]

The implementation clips (R) to a specified interval and applies a reduced rebound factor when a boundary is reached.

The Morse potential, force, levels, bond observables, and reduced velocity-Verlet coupling are implemented explicitly.
5. Structured light and vortices
5.1 Hermite-Gaussian mode

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
]

This is a transverse waist-plane profile rather than a complete propagating Maxwell solution.

5.2 Laguerre-Gaussian mode

Let

[
r=\sqrt{x^2+y^2},
\qquad
\theta=\mathrm{atan2}(y,x),
]

and

[
\rho_L=\frac{2r^2}{w^2}.
]

Then

\left(
\frac{\sqrt2r}{w}
\right)^{|\ell|}
L_p^{|\ell|}(\rho_L)
e^{-r^2/w^2}
e^{i(\ell\theta-\omega t+\phi_0)}.
]

The phase factor

[
e^{i\ell\theta}
]

provides the programmed optical-vortex winding.

5.3 Optical coupling channels

For mixed HG/LG operation, the reduced scalar potential is

A_Lc_s
\left[
\mathrm{Re}(E^{\mathrm{HG}})
+
\mathrm{Re}(E^{\mathrm{LG}})
\right].
]

The reduced spin-mixing drive is

A_Lc_{\mathrm{spin}}
\left[
\mathrm{Re}(E^{\mathrm{HG}})
+
\mathrm{Im}(E^{\mathrm{LG}})
\right].
]

These are phenomenological coupling channels rather than transition-dipole or Maxwell-derived interaction Hamiltonians.

5.4 Static vortex coordinates

For vortex (j),

(x-x_j)^2+(y-y_j)^2,
]

\mathrm{atan2}(y-y_j,x-x_j).
]

The accumulated phase is

\sum_jq_j\theta_j.
]

The depleted amplitude envelope is

\prod_j
\mathrm{clip}
\left[
1-d_j
e^{-r_j^2/(2\xi_j^2)},
0.04,
1
\right].
]

The full implementation imprints

[
u
\longrightarrow
A_v e^{i\Phi_v}u,
]

[
d
\longrightarrow
A_v e^{-0.55i\Phi_v}d.
]

The second component therefore receives a project-specific phase-defect texture rather than a conventional integer-winding vortex.

5.5 Moving-vortex drives

Define

e^{-r_j^2/(2\xi_j^2)}.
]

The scalar drive is

\kappa_v
\sum_jd_jG_j,
]

and the spin drive is

\kappa_v
\sum_j
d_jq_jG_j\sin\theta_j.
]

The current source explicitly labels these moving-vortex terms as phenomenological perturbations.

6. Nonlinear spinor-field evolution
6.1 Total potential and spin drive

V_{\mathrm{mol}}
+
V_{\mathrm{light}}
+
V_{\mathrm{vortex}}
+
V_{\mathrm{vertex}},
]

\Omega_0
+
s_{\mathrm{light}}
+
s_{\mathrm{vortex}}
+
s_{\mathrm{vertex}}.
]

Here (\Omega_0) is the parameter named zeeman in the code. It multiplies (\sigma_x), so it acts as a transverse mixing term rather than a (\sigma_z) energy splitting.

6.2 Compact Hamiltonian

i\lambda
\left(
\sigma_xD_x+\sigma_yD_y
\right).
]

6.3 Density-relaxation term

\gamma
\frac{
\rho_{\mathrm{target}}-\rho
}{
1+\rho
}.
]

6.4 Evolution equation

-iH[\Psi,t]\Psi
+
\Gamma(\rho)\Psi.
]

Equivalently,

H[\Psi,t]\Psi
+
i\Gamma(\rho)\Psi.
]

6.5 Component form

For the upper component,

i\lambda D_xd

\lambda D_yd
+
i\Gamma u.
]

For the lower component,

i\lambda D_xu
+
\lambda D_yu
+
i\Gamma d.
]

The component equations follow directly from the implemented scalar, spin-mixing, derivative-coupling, and gain terms.

7. Finite differences and Strang propagation
7.1 Periodic Laplacian

4f_{ij}
}{
(\Delta x)^2
}.
]

7.2 Centered derivatives

\frac{
f_{i+1,j}-f_{i-1,j}
}{
2\Delta x
},
]

\frac{
f_{i,j+1}-f_{i,j-1}
}{
2\Delta x
}.
]

7.3 Exact finite-difference Fourier symbols

\frac{\sin(k_x\Delta x)}{\Delta x},
]

\frac{\sin(k_y\Delta x)}{\Delta x},
]

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
]

These are the symbols of the existing finite-difference operators, not a replacement continuous (k^2) Laplacian.

7.4 Kinetic and SOC Fourier operator

\epsilon_h(\mathbf k)I
+
\lambda
\left(
\xi_x\sigma_x+\xi_y\sigma_y
\right).
]

Define

\lambda
\sqrt{
\xi_x^2+\xi_y^2
}.
]

The exact half-step is

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
]

7.5 Exact local substep

i\sin(\Delta t\Omega)\sigma_x
\right].
]

7.6 Second-order Strang composition

[
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
]

With zero damping and no repeated mass projection, each subflow is unitary up to floating-point error.

7.7 LG phase-kick rate

\delta t,
k_{\phi}
A_L
E_{\mathrm{LG}}^{\mathrm{real}}.
]

The components are updated by

[
u
\longrightarrow
e^{i\theta_{\mathrm{kick}}}u,
]

[
d
\longrightarrow
e^{-i\theta_{\mathrm{kick}}/2}d.
]

8. Density, bonding, information, and geometric diagnostics
8.1 Density and collective phase

|u|^2+|d|^2,
]

\arg(u+d).
]

8.2 Phase-modulated density

\rho
\left[
1+0.28\cos\phi
\right].
]

This is a nonnegative phase-modulated density display.

8.3 Signed phase quadrature

\rho\cos\phi.
]

Negative (Q_{\rho}) means phase inversion. It does not mean negative electron density, and it is not a Wigner function.

8.4 Bonding projection

\sum_{s\in{\uparrow,\downarrow}}
\left|
\langle\phi_b,\psi_s\rangle_h
\right|^2.
]

8.5 Antibonding projection

\sum_{s\in{\uparrow,\downarrow}}
\left|
\langle\phi_{ab},\psi_s\rangle_h
\right|^2.
]

8.6 Transition score

\frac{
n_{ab}
}{
n_b+n_{ab}+\varepsilon
}.
]

8.7 Effective bond-order proxy

\frac{
n_b-n_{ab}
}{
n_b+n_{ab}+\varepsilon
}.
]

8.8 Two-site sharing proxy

Define

\sum_s
\left|
\langle\phi_A,\psi_s\rangle_h
\right|^2,
]

\sum_s
\left|
\langle\phi_B,\psi_s\rangle_h
\right|^2,
]

\frac{n_B}{n_A+n_B+\varepsilon}.
]

The sharing proxy is

\mathrm{clip}
\left[
4p_Ap_B\max(S,0),
0,
1
\right].
]

These are projection and sharing diagnostics, not ab-initio bond orders or orbital entanglement measurements.

9. Bloch geometry and topology
9.1 Local Bloch vector

\frac{
2\mathrm{Re}(u^{*}d)
}{
\rho+\varepsilon
},
]

\frac{
2\mathrm{Im}(u^{*}d)
}{
\rho+\varepsilon
},
]

\frac{
|u|^2-|d|^2
}{
\rho+\varepsilon
}.
]

Collectively,

\frac{
\Psi^{\dagger}
\boldsymbol{\sigma}
\Psi
}{
\Psi^{\dagger}\Psi+\varepsilon
}.
]

9.2 Spin-texture curvature diagnostic

|\nabla\mathbf n|^2,
]

or discretely,

\sum_{a=x,y,z}
\left[
(D_xn_a)^2+(D_yn_a)^2
\right].
]

This is spin-texture gradient energy, not spacetime curvature.

9.3 Topological-charge density

\frac{1}{4\pi}
\mathbf n
\cdot
\left(
\partial_x\mathbf n
\times
\partial_y\mathbf n
\right).
]

9.4 Integrated charge

\int_{\Omega}
q(x,y,t),dA.
]

On the masked grid,

\sum_{i,j}
m_{ij}
q_{ij}
(\Delta x)^2.
]

Because low-density points and derivative neighbors are masked, the computed value need not be exactly integer.

10. QFI, entropy, memory, OAM, and winding
10.1 Spin-(z) expectation

\frac{
\int
\left(
|u|^2-|d|^2
\right)dA
}{
\int\rho,dA
}.
]

10.2 QFI diagnostic

Since

[
\sigma_z^2=I,
]

the implemented pure-state-style diagnostic is

4
\left[
1-\langle\sigma_z\rangle^2
\right].
]

10.3 Spatial entropy diagnostic

Define

\frac{
\rho_{ij}
}{
\sum_{k,l}\rho_{kl}+\varepsilon
}.
]

Then

-\sum_{i,j}
p_{ij}\ln(p_{ij}+\varepsilon).
]

The QFI and entropy formulas are implemented as relative diagnostics.

10.4 BLP-inspired density distinguishability

For two nearby spinor trajectories (a) and (b),

\frac{
\frac12
\int
|\rho_a-\rho_b|,dA
}{
\frac12
\int
(\rho_a+\rho_b),dA
+
\varepsilon
}.
]

10.5 Memory current

\frac{
D_{\rho}(t_n)-D_{\rho}(t_{n-1})
}{
\Delta t
}.
]

A positive value indicates revival of density distinguishability in this reduced diagnostic. It is BLP-inspired, but it is not the canonical trace distance between optimized reduced density matrices.

10.6 OAM expectation

\frac{
\sum_s
\mathrm{Im}
\left[
\int
\psi_s^{*}
\left(
x\partial_y-y\partial_x
\right)
\psi_s,dA
\right]
}{
\int\rho,dA
}.
]

This quantity is expressed in reduced units.

10.7 Boundary phase winding

\frac{1}{2\pi}
\oint_{\partial\Omega}
d\phi.
]

The code evaluates the wrapped phase differences around the periodic grid boundary.

11. Toy MEAO/QIT reduced-state diagnostics

The local Fock-space labels are

[
{0,\uparrow,\downarrow,\uparrow\downarrow}
]

for each of two orbital subsystems.

These formulas describe the live reduced proxy only. They are not a validated MEAO orbital optimization and are not time-resolved ab-initio entanglement values.

11.1 Covalent reference state

|\downarrow,\uparrow\rangle
+
|\uparrow\downarrow,0\rangle
\right).
]

11.2 Antibonding-contaminated reference

|\uparrow,\downarrow\rangle
+
|\downarrow,\uparrow\rangle
+
|\uparrow\downarrow,0\rangle
\right).
]

11.3 Ionic reference

For molecular polarity (p),

\mathcal N_I
\left[
\sqrt{\frac{1+p}{2}}
|0,\uparrow\downarrow\rangle
+
\sqrt{\frac{1-p}{2}}
|\uparrow\downarrow,0\rangle
\right].
]

11.4 Ionic mixing coefficient

\mathrm{clip}
\left[
0.12|p|
+
0.22T
+
0.12\max(0,1-B),
0,
0.90
\right].
]

11.5 Toy two-orbital state

Using

[
[x]_+=\max(0,x),
]

the proxy state is

\mathcal N
\left[
\sqrt{
[1-\iota-0.15T]+
}
|\psi{\mathrm{cov}}\rangle
+
\sqrt{\iota}
|\psi_{\mathrm{ionic}}\rangle
+
\sqrt{0.15T}
|\psi_{\mathrm{ab}}\rangle
\right].
]

11.6 Reduced density operator

|\psi_{\mathrm{toy}}\rangle
\langle\psi_{\mathrm{toy}}|.
]

11.7 Partial traces

\mathrm{Tr}B(\rho{AB}),
]

\mathrm{Tr}A(\rho{AB}).
]

11.8 Von Neumann entropy

-\mathrm{Tr}
\left[
\rho\ln(\rho+\varepsilon)
\right].
]

11.9 Mutual information

S(\rho_{AB}).
]

The normalized proxy is

\frac{
I(A)
}{
2\ln4
}.
]

11.10 Logarithmic negativity

\max
\left[
0,
\log_2
\left(
|\rho_{AB}^{T_B}|_1+\varepsilon
\right)
\right].
]

11.11 Four-level I-concurrence-style diagnostic

\frac{
\sqrt{
2
\left[
1-\mathrm{Tr}(\rho_A^2)
\right]
}
}{
\sqrt{
2(1-1/4)
}
}.
]

This is not Wootters' two-qubit concurrence.

11.12 MEAO coherence proxy

\left|
\rho_{(0,\uparrow\downarrow),
(\uparrow\downarrow,0)}
\right|^2
+
\left|
\rho_{(\uparrow,\downarrow),
(\downarrow,\uparrow)}
\right|^2.
]

The normalized form is

\mathrm{clip}
\left[
\frac{
F_{\mathrm{MEAO}}^{\mathrm{proxy}}
}{
0.125
},
0,
1.5
\right].
]

11.13 Fractional bond-order proxy

m_0
\left[
0.50\max(0,B)
+
0.30I_{\mathrm{norm}}
+
0.20\min(F_{\mathrm{norm}},1)
\right]
(1-0.18T),
]

followed by scenario-dependent clipping.

The complete state construction and the resulting mutual-information, negativity, concurrence, coherence, and fractional-bond formulas are implemented in the toy MEAO/QIT analyzer.

12. Clifford-algebra encoding
12.1 Clifford relation

The real Clifford generators satisfy

2g_{ij}.
]

For (\mathrm{Cl}(3,0)),

\mathrm{diag}(1,1,1).
]

For (\mathrm{Cl}(1,3)),

\mathrm{diag}(1,-1,-1,-1).
]

The compact engine represents (n=p+q) generators with

[
\dim\mathrm{Cl}(p,q)=2^n.
]

It implements genuine geometric products and reversion.

12.2 Density-weighted global spinor

Define

\frac{
\rho_{ij}
}{
\sum_{k,l}\rho_{kl}+\varepsilon
}.
]

Then

\sum_{i,j}
w_{ij}u_{ij},
]

\sum_{i,j}
w_{ij}d_{ij}.
]

Normalize:

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
\alpha_{\mathrm{raw}}\
\beta_{\mathrm{raw}}
\end{pmatrix}.
]

12.3 Complex coefficients

[
\alpha=a+ib,
\qquad
\beta=c+id.
]

12.4 Even subalgebra

\mathrm{span}
{1,e_{12},e_{23},e_{31}}.
]

This algebra is quaternion-like.

12.5 Spinor-to-Clifford map

a
+
b e_{12}
+
c e_{23}
+
d e_{31}.
]

12.6 Primitive idempotent

\frac{1+e_3}{2},
\qquad
f^2=f.
]

12.7 Minimal-left-ideal representative

\Phi f.
]

12.8 Clifford density element

\Phi_L
\widetilde{\Phi_L}.
]

For a grade-(k) component,

(-1)^{k(k-1)/2}A_k.
]

12.9 Rotor-angle diagnostic

2,\mathrm{atan2}
\left(
|\langle\Phi\rangle_2|,
|\langle\Phi\rangle_0|+\varepsilon
\right).
]

12.10 Purity-like Clifford diagnostic

Let

[
s_c=\langle\rho_c\rangle_0,
]

[
v_c=|\langle\rho_c\rangle_1|.
]

Then

\frac{
s_c^2+v_c^2
}{
|\rho_c|^2+\varepsilon
}.
]

The source implements this encoding as a representation and descriptor system, not as the PDE time-evolution operator.

13. Feature chart, PCA, and ridge-DMD

At each stored time, the project forms a diagnostic vector

\begin{pmatrix}
\text{memory current}\
\text{QFI}\
\text{curvature}\
\text{OAM}\
\text{phase winding}\
\text{bond order}\
\text{MEAO/QIT proxies}\
\text{bond coordinate}\
\text{Clifford descriptors}\
\text{laser controls}\
\text{local patch descriptors}
\end{pmatrix}.
]

The actual feature chart includes both global observables and compact spatial Clifford summaries.

13.1 Training-only standardization

\frac{1}{N_{\mathrm{train}}}
\sum_{t\in\mathrm{train}}
\mathbf x_t,
]

\mathrm{std}_{t\in\mathrm{train}}
(\mathbf x_t),
]

\frac{
\mathbf x_t-\boldsymbol{\mu}
}{
\mathbf s
}.
]

13.2 Training-block SVD

U\Sigma V^{T}.
]

Let (W) contain the first (k=3) rows of (V^{T}).

13.3 PCA latent coordinates

\widehat{\mathbf x}_tW^{T}.
]

These are learned diagnostic coordinates. They are not gamma-matrix axes and are not an intrinsic manifold reconstruction.

13.4 Ridge-DMD transition matrix

Construct

\begin{pmatrix}
\mathbf z_0\
\mathbf z_1\
\vdots\
\mathbf z_{N_{\mathrm{train}}-2}
\end{pmatrix},
]

\begin{pmatrix}
\mathbf z_1\
\mathbf z_2\
\vdots\
\mathbf z_{N_{\mathrm{train}}-1}
\end{pmatrix}.
]

Then

\left(
X^{T}X+\lambda_RI
\right)^{-1}
X^{T}Y.
]

13.5 One-step prediction

\mathbf z_tA_{\mathrm{DMD}}.
]

13.6 One-step error

\widehat{\mathbf z}_{t+1}
|_2.
]

13.7 Recursive rollout

\mathbf z_t
A_{\mathrm{DMD}}^h.
]

13.8 Skill against persistence

1-
\frac{
\mathrm{RMSE}{\mathrm{DMD}}
}{
\mathrm{RMSE}{\mathrm{persistence}}
}.
]

The scaler, PCA basis, and DMD operator are fitted only on an early contiguous block and evaluated on later held-out states.

14. Hopf geometry and (S^3) state-space fibers
14.1 Normalized complex spinor on (S^3)

Let

\begin{pmatrix}
z_1\
z_2
\end{pmatrix},
]

with

[
|z_1|^2+|z_2|^2=1.
]

The four real coordinates are

\begin{pmatrix}
\mathrm{Re},z_1\
\mathrm{Im},z_1\
\mathrm{Re},z_2\
\mathrm{Im},z_2
\end{pmatrix}
\in S^3\subset R^4.
]

14.2 Hopf map

For a possibly unnormalized pair,

\frac{1}{
|z_1|^2+|z_2|^2
}
\begin{pmatrix}
2\mathrm{Re}(z_1^{}z_2)\
2\mathrm{Im}(z_1^{}z_2)\
|z_1|^2-|z_2|^2
\end{pmatrix}.
]

For a normalized spinor, the denominator is one.

14.3 Local chart phase

The source also records the chart-local phase

\frac12\arg(ud).
]

This is not a globally defined fiber coordinate at component zeros.

14.4 Canonical lift from (S^2)

For

\begin{pmatrix}
\sin\theta\cos\varphi\
\sin\theta\sin\varphi\
\cos\theta
\end{pmatrix},
]

one canonical lift is

\cos\frac{\theta}{2},
]

e^{i\varphi}
\sin\frac{\theta}{2}.
]

14.5 Complete (U(1)) fiber

Every common phase maps to the same Bloch point:

[
(z_1,z_2)
\longrightarrow
e^{i\alpha}(z_1,z_2).
]

The complete fiber is

{
e^{i\alpha}(z_1,z_2)
:
0\leq\alpha<2\pi
}.
]

In real (R^4) coordinates,

u\cos\alpha+v\sin\alpha.
]

14.6 Clifford torus coordinates

\cos\eta,e^{i\xi_1},
]

\sin\eta,e^{i\xi_2}.
]

Therefore,

\begin{pmatrix}
\cos\eta\cos\xi_1\
\cos\eta\sin\xi_1\
\sin\eta\cos\xi_2\
\sin\eta\sin\xi_2
\end{pmatrix}.
]

The corresponding Bloch latitude satisfies

[
n_z=\cos(2\eta).
]

14.7 Stereographic projection

After a fixed (R^4) rotation (A), let

[
q'=Aq.
]

The raw (R^3) stereographic coordinate is

\frac{
(q'_1,q'_2,q'_3)
}{
1-q'_4
}.
]

14.8 Finite display compaction

The finite display coordinate is

\frac{
R_c\mathbf p
}{
R_c+|\mathbf p|
}.
]

These are complete (U(1)) state-space fibers projected for visualization. They are not physical-space three-dimensional Hopfion solitons.

15. Future braid-resolved information experiments

This section is proposed future work. It is not a current measured entanglement observable.

15.1 Braid word

For resolved temporal crossings,

\sigma_{i_1}^{\epsilon_1}
\sigma_{i_2}^{\epsilon_2}
\cdots
\sigma_{i_m}^{\epsilon_m},
]

where

[
\epsilon_k\in{-1,+1}
]

records the orientation of each crossing.

15.2 Information change over lag (\tau)

I(t+\tau)-I(t).
]

15.3 Braid-feature prediction model

\beta_0
+
\sum_{j=1}^{p}
\beta_jT_j(t)
+
\varepsilon(t+\tau).
]

Here (T_j(t)) may include:

crossing rate;
braid-word complexity;
braid entropy;
phase holonomy;
writhe-like descriptors;
motif persistence.
15.4 Controlled prediction model

To separate braiding from shared external drivers,

\beta_0
+
\boldsymbol{\beta}^{T}\mathbf T(t)
+
\boldsymbol{\gamma}^{T}\mathbf C(t)
+
\varepsilon(t+\tau),
]

where (\mathbf C(t)) includes controls such as:

molecule identity;
current (I(t));
laser amplitude;
LG index (\ell);
vortex radius and depth;
bond coordinate (R(t));
curvature;
dissipative or conservative regime.

The key test is whether

[
\text{controls + braid features}
]

predict held-out information changes better than

[
\text{controls only}.
]

15.5 Future normal-mode extension

For a nonlinear (N)-atom molecule,

[
q_k(t),
\qquad
k=1,\ldots,3N-6.
]

For a linear molecule,

[
q_k(t),
\qquad
k=1,\ldots,3N-5.
]

This would replace the present single bond-stretch coordinate (R(t)) with a vector of normal-mode coordinates.

Scientific interpretation boundary

The implemented hierarchy is

[
\text{PDE state}
\rightarrow
\text{derived diagnostics}
\rightarrow
\text{Clifford/PCA representation}
\rightarrow
\text{Hopf visualization}.
]

It does not currently establish

\text{molecular entanglement}.
]

Any future braid relationship must be reported separately as:

[
\text{correlation},
]

[
\text{held-out prediction},
]

or

[
\text{causal evidence}.
]

