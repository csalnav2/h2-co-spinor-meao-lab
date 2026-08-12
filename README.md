


# MEAO Spinor Lab

## A Reduced-Order Spinor-Field Framework for H₂ and CO: MEAO-Inspired Bonding, Information Flow, and Hopf-Fiber Diagnostics

[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Apache--2.0-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-research%20prototype-orange.svg)](#scientific-scope-and-claim-boundaries)

**MEAO Spinor Lab** is a public research and collaboration scaffold for a private, physics-informed simulation of **molecularly shaped two-component spinor fields** in analytic H₂ and CO landscapes.

The project asks a focused question:

> How do structured-light and topological perturbations reorganize a molecularly shaped spinor field, and which reduced bonding, information, geometric, and memory diagnostics co-evolve during that response?

The complete private framework combines analytic diatomic orbital shapes, nonlinear spinor dynamics, programmable Gaussian-family light, vortices, MEAO-inspired reduced-state diagnostics, Clifford-algebra descriptors, held-out latent forecasting, and Hopf-fiber state-space visualization.

This public repository provides:

- the scientific framing and equations;
- explicit claim boundaries;
- a runnable abbreviated Python concept demo;
- deterministic CSV, JSON, and image outputs;
- tests and continuous-integration configuration;
- a pathway for research collaboration without publishing the complete private solver.

> **Scope in one sentence:** this is a reduced-order hypothesis-generation and diagnostic framework, not a replacement for ab-initio electronic-structure theory.

---

## Table of contents

- [Why this project exists](#why-this-project-exists)
- [Conceptual pipeline](#conceptual-pipeline)
- [Core mathematical model](#core-mathematical-model)
- [Structured-light and vortex controls](#structured-light-and-vortex-controls)
- [Diagnostic layers](#diagnostic-layers)
- [Clifford encoding and latent dynamics](#clifford-encoding-and-latent-dynamics)
- [Hopf fibers and future braid diagnostics](#hopf-fibers-and-future-braid-diagnostics)
- [What the public demo implements](#what-the-public-demo-implements)
- [Quick start](#quick-start)
- [Generated outputs](#generated-outputs)
- [Scientific scope and claim boundaries](#scientific-scope-and-claim-boundaries)
- [Validation philosophy](#validation-philosophy)
- [Repository structure](#repository-structure)
- [Roadmap](#roadmap)
- [Collaboration](#collaboration)
- [Citation and references](#citation-and-references)
- [License](#license)

---

## Why this project exists

Electronic-structure methods such as Hartree-Fock, post-Hartree-Fock active-space calculations, density-functional theory, and explicitly correlated dynamics provide the quantitative foundation for molecular physics. They are indispensable, but full time-dependent calculations under broad control sweeps can become expensive, especially when the goal is exploratory visualization, controller design, or rapid comparison of many perturbation regimes.

MEAO Spinor Lab explores a complementary layer:

1. represent H₂ and CO with analytic, molecule-shaped two-center envelopes;
2. initialize a two-component Pauli spinor field;
3. perturb it with continuous Gaussian, Hermite-Gaussian, and Laguerre-Gaussian fields and programmable vortices;
4. evolve a nonlinear reduced-order wave model;
5. monitor reduced bonding, information, memory, spin, topology, Clifford, and latent diagnostics;
6. use separately labelled quantum-chemistry calculations as calibration or validation anchors when appropriate.

The objective is not to replace quantum chemistry. It is to build an **instrument panel for hypotheses**: a fast environment in which unexpected correlations, control responses, numerical pathologies, or candidate topological signatures can be identified before committing to more expensive calculations or experiments.

---

## Conceptual pipeline

```text
Analytic H₂ / CO geometry
          │
          ▼
Two-center atomic envelopes ──► bonding / antibonding channels
          │
          ▼
Two-component spinor field Ψ = (ψ↑, ψ↓)ᵀ
          │
          ├── Gaussian / HG / LG structured light
          ├── optical OAM index ℓ
          ├── integer or programmable vortices
          ├── molecular potential and bond coordinate
          └── optional conservative / dissipative terms
          │
          ▼
Reduced nonlinear spinor propagation
          │
          ├── density, phase, signed phase quadrature
          ├── bond-order and transition proxies
          ├── MEAO-inspired information proxies
          ├── BLP-inspired density distinguishability
          ├── QFI, OAM, Bloch texture, curvature, topology
          ├── Clifford minimal-left-ideal descriptors
          └── local patch / low-frequency summaries
          │
          ▼
Train-only PCA ──► three-coordinate diagnostic latent trajectory
          │
          ▼
Ridge-DMD prediction on held-out future states
          │
          ▼
S² Bloch views and selected S³ Hopf U(1) fibers
```

Three layers must remain conceptually separate:

| Layer | Role | What it is not |
|---|---|---|
| Spinor PDE | Evolves the reduced field | An ab-initio time-dependent many-electron calculation |
| Clifford encoding | Re-expresses spinor geometry and produces descriptors | The time propagator or gamma-matrix dynamics |
| PCA/ridge-DMD | Compresses and forecasts diagnostic trajectories | An intrinsic manifold reconstruction or proof of a symmetry |

---

## Core mathematical model

### 1. Two-component spinor field

The reduced state is

$$
\Psi(x,y,t)=
\begin{pmatrix}
\psi_\uparrow(x,y,t)\\
\psi_\downarrow(x,y,t)
\end{pmatrix},
$$

with total field density

$$
\rho(x,y,t)=|\psi_\uparrow|^2+|\psi_\downarrow|^2.
$$

Depending on the selected normalization convention, $\rho$ can be interpreted either as a unit-normalized effective probability field or as the intensity of a nonlinear mean-field order parameter. The public demo exposes both conventions explicitly:

$$
\int_\Omega \rho\,dA=1
\qquad\mathrm{or}\qquad
\frac{1}{A}\int_\Omega\rho\,dA=\rho_{\mathrm{target}}.
$$

### 2. Reduced nonlinear evolution

The complete private solver uses a controlled Pauli/Gross-Pitaevskii-like model with kinetic, molecular, nonlinear, spin-mixing, spin-orbit-like, optical, vortex, and optional dissipative terms. A compact schematic form is

$$
i\,\partial_t\Psi=
[
-\kappa\nabla^2
+V_{\mathrm{mol}}
+V_{\mathrm{light}}
+g\rho
+\Omega(x,y,t)\sigma_x
+H_{\mathrm{SOC}}
]\Psi
+i\Gamma(\rho)\Psi.
$$

The public mini-demo deliberately omits the private spin-orbit kernel, controller, production parameterization, and checkpoint/rendering system. Its abbreviated conservative propagator is

$$
i\,\partial_t\Psi=
[
-\kappa\nabla^2
+V(x,y,t)
+g\rho
+\Omega(x,y,t)\sigma_x
]\Psi,
$$

advanced with a unitary Strang composition:

$$
U(\Delta t)\approx
U_K\!(\frac{\Delta t}{2})
U_{\mathrm{local}}(\Delta t)
U_K\!(\frac{\Delta t}{2}).
$$

### 3. Molecular scaffold

The public demo constructs softened Slater-like two-center envelopes,

$$
\phi_A\propto e^{-\zeta_A r_A},
\qquad
\phi_B\propto e^{-\zeta_B r_B},
$$

then forms bonding and antibonding channels,

$$
\phi_b\propto c_A\phi_A+c_B\phi_B,
$$

$$
\phi_{ab}\propto c_B\phi_A-c_A\phi_B.
$$

For the public concept code, the antibonding channel is Gram-Schmidt orthogonalized against the bonding channel so that the projection diagnostics are unambiguous.

A reduced Morse coordinate supplies an illustrative bond plant,

$$
V_M(R)=D_e(1-e^{-a(R-R_e)})^2-D_e.
$$

These orbitals and parameters are intentionally generic. They are **molecular shape models**, not basis-set electronic-structure solutions.

---

## Structured-light and vortex controls

### Gaussian and Hermite-Gaussian modes

The public demo includes Gaussian and HG envelopes of the form

$$
E_{mn}^{\mathrm{HG}}(x,y,t)
\propto
H_m\!(\frac{\sqrt2x}{w})
H_n\!(\frac{\sqrt2y}{w})
e^{-(x^2+y^2)/w^2}
e^{-i\omega t}.
$$

### Laguerre-Gaussian modes and orbital angular momentum

For radial index $p=0$, the abbreviated LG field is

$$
E_{0\ell}^{\mathrm{LG}}(r,\theta,t)
\propto
(\frac{\sqrt2r}{w})^{|\ell|}
e^{-r^2/w^2}
e^{i(\ell\theta-\omega t)}.
$$

The azimuthal phase factor $e^{i\ell\theta}$ provides the optical orbital-angular-momentum control parameter $\ell$.

### Vortices

The public code imprints an integer spin vortex,

$$
\psi_\uparrow\rightarrow A(r)e^{iq\theta}\psi_\uparrow,
\qquad
\psi_\downarrow\rightarrow A(r)e^{-iq\theta}\psi_\downarrow,
$$

with a smooth depleted core $A(r)$. The complete private project explores a broader programmable vortex parameter space, including population, core size, depth, geometry, orbit radius, angular velocity, and drive strength.

### Current versus planned optical channels

| Capability | Status |
|---|---|
| Continuous Gaussian/HG/LG fields | Implemented in the research framework and public mini-demo |
| Mixed HG/LG control | Implemented |
| Optical OAM index $\ell$ | Implemented |
| Explicit pulse envelopes | Planned extension |
| Explicit one-photon transition model | Planned extension |
| Explicit two-photon interaction model | Planned extension |

A genuine two-photon model will require an explicit second-order interaction structure, pulse bandwidth and timing, intermediate-state or effective-transition physics, and calibrated couplings. It should not be represented merely by squaring a field amplitude.

---

## Diagnostic layers

### Density, phase, and signed quadrature

The field density is

$$
\rho=|\psi_\uparrow|^2+|\psi_\downarrow|^2.
$$

The aggregate phase used in the display is

$$
\phi=\arg(\psi_\uparrow+\psi_\downarrow),
$$

and the signed phase quadrature is

$$
Q_\rho=\rho\cos\phi.
$$

Negative $Q_\rho$ indicates phase inversion. It is **not negative electron density**.

### Bonding and transition proxies

Projecting the field onto orthonormalized bonding and antibonding channels gives

$$
n_b=\sum_{s\in\{\uparrow,\downarrow\}}
|\langle\phi_b|\psi_s\rangle|^2,
$$

$$
n_{ab}=\sum_{s\in\{\uparrow,\downarrow\}}
|\langle\phi_{ab}|\psi_s\rangle|^2.
$$

The public demo then reports

$$
B_{\mathrm{proxy}}=
\frac{n_b-n_{ab}}{n_b+n_{ab}+\varepsilon},
$$

and

$$
T_{\mathrm{proxy}}=
\frac{n_{ab}}{n_b+n_{ab}+\varepsilon}.
$$

These are reduced control diagnostics, not experimentally calibrated bond orders or electronic transition probabilities.

### MEAO-inspired information diagnostics

**MEAO** means **Maximally Entangled Atomic Orbitals**. In the private live PDE pathway, mutual information, concurrence, logarithmic negativity, and related values are explicitly labelled **MEAO/QIT-inspired reduced-state proxies** unless they are derived from a separately identified quantum-chemistry calculation.

The public mini-demo uses intentionally simple site-sharing proxies to illustrate the data pathway. They are not orbital reduced-density-matrix entanglement measures.

A production quantum-chemistry pathway is conceptually

```text
geometry
  → basis and molecular orbitals
  → localized orbitals / IAOs / Meta-Löwdin orbitals
  → correlated active-space wavefunction
  → 1RDM / 2RDM and one-/two-orbital reduced states
  → validated MEAO optimization and information graph
```

### BLP-inspired memory current

The public demo evolves two nearby spinor preparations and computes a normalized density-field distance,

$$
D_\rho(t)=
\frac{1}{2\mathcal N}
\int|\rho_a(x,y,t)-\rho_b(x,y,t)|\,dA,
$$

followed by

$$
\sigma_\rho(t)=\frac{dD_\rho}{dt}.
$$

Positive intervals can indicate a revival of **density distinguishability** in this reduced model. This is BLP-inspired, but it is not the canonical Breuer-Laine-Piilo optimization over reduced density operators.

### Bloch texture and QFI

The local Bloch vector is

$$
\mathbf n(x,y,t)=
\frac{\Psi^\dagger\boldsymbol\sigma\Psi}
{\Psi^\dagger\Psi+\varepsilon}.
$$

For a normalized local spinor, $\mathbf n$ lies on $S^2$. The public demo also reports a spin-$z$ QFI-like quantity for the aggregate state,

$$
F_Q[\sigma_z]=4(\langle\sigma_z^2\rangle-\langle\sigma_z\rangle^2).
$$

### Curvature and topological charge

The spin-texture gradient diagnostic is

$$
C(x,y,t)=|\nabla\mathbf n|^2.
$$

This is a geometric spin-texture diagnostic. It is not spacetime curvature.

The density-masked topological-charge form used by the research framework is

$$
q(x,y,t)=
\frac{1}{4\pi}
\mathbf n\cdot
(\partial_x\mathbf n\times\partial_y\mathbf n),
$$

with integrated charge

$$
Q=\int q(x,y,t)\,dA.
$$

### OAM expectation

The abbreviated field OAM diagnostic is

$$
\langle L_z\rangle
\propto
\frac{
\sum_s\mathrm{Im}
(\int\psi_s^*(x\partial_y-y\partial_x)\psi_s\,dA)
}{\int\rho\,dA}.
$$

It is expressed in reduced units and should not be read as a calibrated photon or electron angular momentum measurement.

---

## Clifford encoding and latent dynamics

### Pauli spinor to the even subalgebra of $\mathrm{Cl}(3,0)$

For a normalized global spinor

$$
\psi=
\begin{pmatrix}
\alpha\\
\beta
\end{pmatrix},
\qquad
\alpha=a+ib,
\qquad
\beta=c+id,
$$

the Clifford encoder uses

$$
\Phi=a+b e_{12}+c e_{23}+d e_{31}.
$$

The even subalgebra

$$
\mathrm{Cl}^+(3,0)=
\mathrm{span}\{1,e_{12},e_{23},e_{31}\}
$$

is quaternion-like. A minimal-left-ideal representative is formed with

$$
f=\frac{1+e_3}{2},
\qquad
\Phi_L=\Phi f,
$$

and the Clifford density element is

$$
\rho_c=\Phi_L\widetilde{\Phi_L},
$$

For any multivector $M$, $\widetilde{M}$ denotes reversion.

The public demo includes a compact, generic $\mathrm{Cl}(3,0)$ geometric-product engine and extracts scalar, vector, bivector, pseudoscalar, rotor-angle, and purity-like descriptors.

### What this does not mean

The Clifford encoder is **not** the time-evolution operator. The code does not propagate the field by sequentially multiplying a latent vector by $\gamma^0,\gamma^1,\gamma^2,\gamma^3$.

The actual order is

$$
\Psi_t
\longrightarrow
\Psi_{t+\Delta t}
\longrightarrow
x_t
\longrightarrow
z_t.
$$

The three arrows denote, in order, the PDE update, diagnostic extraction, and PCA projection.

### PCA compression

A diagnostic chart is assembled,

$$
x_t=
\begin{pmatrix}
B_{\mathrm{proxy}}(t)\\
I_{\mathrm{proxy}}(t)\\
F_Q(t)\\
C(t)\\
\rho_{c,0}(t)\\
\vdots
\end{pmatrix}
\in\mathbb R^p.
$$

After train-only standardization, PCA gives

$$
z_t=W\frac{x_t-\mu}{s},
\qquad
z_t\in\mathbb R^3.
$$

PC1, PC2, and PC3 are learned mixtures of observables. They are not gamma matrices, physical coordinate axes, or guaranteed rotational generators.

### Ridge-DMD

A regularized linear map is fitted only on the early contiguous block,

$$
z_{t+1}\approx z_tA,
$$

with

$$
A=
(X^TX+\lambda I)^{-1}X^TY.
$$

Later states are held out for validation. The DMD matrix may contain rotation, shear, contraction, expansion, or mixed behavior. It is not constrained to be a quaternion rotor or Lorentz transformation.

---

## Hopf fibers and future braid diagnostics

### Hopf state-space visualization

A normalized two-component complex spinor

$$
(z_1,z_2),
\qquad
|z_1|^2+|z_2|^2=1,
$$

is a point on $S^3$. The Hopf map sends it to $S^2$:

$$
H(z_1,z_2)=
\begin{pmatrix}
2\mathrm{Re}(z_1^*z_2)\\
2\mathrm{Im}(z_1^*z_2)\\
|z_1|^2-|z_2|^2
\end{pmatrix}.
$$

Every common phase

$$
(z_1,z_2)\rightarrow e^{i\alpha}(z_1,z_2)
$$

maps to the same point on $S^2$. Its complete preimage is a $U(1)$ fiber,

$$
\mathcal F_{\mathbf n}=
\{e^{i\alpha}(z_1,z_2):0\le\alpha<2\pi\}.
$$

The complete private dashboard tracks 236 stable featured $S^2$ display identities and lifts 14 selected representatives to complete Hopf fibers. The public demo constructs one complete fiber from its final global spinor and stereographically projects it from $S^3$ into $\mathbb R^3$.

These fibers are **state-space geometry**. They are not automatically physical-space Hopfion solitons.

### Future direction: braid-resolved diagnostics

The future hypothesis is not that visible braids already measure entanglement. It is:

> Can temporal braid descriptors derived from consistently tracked fiber markers provide robust, supplementary correlates of independently computed molecular information dynamics?

Candidate descriptors include

- crossing-event rate;
- braid words and conjugacy-invariant summaries;
- braid entropy or complexity;
- relative phase holonomy;
- persistence of crossing motifs;
- projection-averaged writhe-like statistics;
- lagged predictability of changes in bond order or information metrics.

A future test could model

$$
\Delta I(t+\tau)
=
\beta_0+
\sum_j\beta_jT_j(t)+
\varepsilon(t),
$$

where $T_j(t)$ are braid features and $I$ is an independently evaluated information quantity.

Credible evidence will require:

1. gauge robustness;
2. projection robustness;
3. fiber-selection bootstrapping;
4. randomized-label and time-shuffle null models;
5. held-out perturbation regimes;
6. independent RDM-based quantum-chemistry validation;
7. explicit separation of correlation, prediction, and causal claims.

A suitable future-work heading is:

> **Toward Braid-Resolved Topological Correlates of Molecular Information Dynamics**

---

## What the public demo implements

The repository provides two public code paths:

- `demo/quick_gist.py`, a compact roughly 300-line architectural sketch for readers who want the idea quickly;
- `demo/mini_spinor_lab.py`, a more defensive public concept implementation with explicit normalization modes, a compact geometric-product engine, richer exports, metadata, and tests.

Neither file reproduces the complete private implementation.

### Included

- H₂-like and CO-like analytic two-center envelopes;
- orthonormalized bonding and antibonding channels;
- a two-component complex spinor field;
- unit-norm or mean-density normalization modes;
- Gaussian, HG, LG, mixed, or disabled optical fields;
- integer spin-vortex imprinting;
- a simplified unitary Strang split-step propagator;
- a reduced Morse bond coordinate;
- density, phase, signed quadrature, and field plots;
- bond-order and transition proxies;
- deliberately labelled MEAO mutual-information and concurrence proxies;
- BLP-inspired density distinguishability and memory current;
- QFI, OAM, Bloch curvature, and topological charge;
- a compact $\mathrm{Cl}(3,0)$ geometric-product engine;
- minimal-left-ideal and Clifford-density descriptors;
- train-only PCA and ridge-DMD;
- one complete Hopf $U(1)$ fiber in $S^3$ with an $\mathbb R^3$ projection;
- deterministic CSV, JSON, PNG, and test outputs.

### Intentionally omitted

- the complete private finite-difference-symbol kinetic/SOC kernel;
- the production H₂/CO parameter set;
- the full MPC controller and action search;
- private checkpoint, cache, provenance, and dashboard orchestration;
- private rendering and selected-fiber ranking logic;
- full local-patch and spectral feature tensors;
- generated research data;
- credentials and network services;
- calibrated one-photon or two-photon pulse physics;
- a validated MEAO orbital-optimization backend;
- a physical-space three-dimensional Hopfion model.

The public values are generic reduced units. They are not a disguised copy of the private parameterization.

---

## Quick start

### 1. Create an environment

```bash
python -m venv .venv
```

Activate it:

```bash
# Linux / macOS
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

### 2. Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### 3. Run the compact architectural gist

```bash
python demo/quick_gist.py \
  --molecule H2 \
  --light lg \
  --ell 1 \
  --grid 48 \
  --steps 100 \
  --out demo_output/quick_gist
```

This produces a single compact figure, scalar diagnostics, Hopf-fiber coordinates, and metadata.

### 4. Run the more complete public H₂-like LG/OAM concept demo

```bash
python demo/mini_spinor_lab.py \
  --molecule H2 \
  --light lg \
  --lg-l 1 \
  --vortex \
  --grid 64 \
  --steps 180 \
  --out demo_output/h2_lg
```

### 5. Run a CO-like mixed HG/LG experiment

```bash
python demo/mini_spinor_lab.py \
  --molecule CO \
  --light mixed \
  --hg 2 1 \
  --lg-l 2 \
  --grid 64 \
  --steps 180 \
  --out demo_output/co_mixed
```

### 6. Compare normalization conventions

Unit-normalized effective state:

```bash
python demo/mini_spinor_lab.py \
  --molecule H2 \
  --normalization-mode unit \
  --out demo_output/h2_unit
```

Mean-field amplitude convention:

```bash
python demo/mini_spinor_lab.py \
  --molecule H2 \
  --normalization-mode mean-density \
  --mean-density 0.05 \
  --out demo_output/h2_mean_field
```

### 7. Enable a vortex or disable the light

```bash
python demo/mini_spinor_lab.py \
  --molecule CO \
  --light lg \
  --vortex \
  --out demo_output/co_vortex

python demo/mini_spinor_lab.py \
  --molecule CO \
  --light off \
  --out demo_output/co_no_light
```

### 8. Run tests

```bash
python -m unittest discover -s tests -v
```

---

## Generated outputs

The compact `quick_gist.py` path writes:

| File | Contents |
|---|---|
| `diagnostics.csv` | Abbreviated scalar timeline |
| `hopf_fiber.csv` | One complete state-space fiber |
| `quick_gist.png` | Density, topology, diagnostics, PCA, and Hopf summary |
| `metadata.json` | DMD matrix, feature names, and claim boundary |

The richer `mini_spinor_lab.py` path writes:

| File | Contents |
|---|---|
| `diagnostics.csv` | Complete scalar time series, including latent coordinates and DMD error |
| `hopf_fiber.csv` | Raw $S^3\subset\mathbb R^4$ fiber coordinates and projected $\mathbb R^3$ coordinates |
| `field_and_diagnostics.png` | Final density, phase, signed quadrature, curvature, topology, and 1D diagnostics |
| `latent_and_hopf.png` | PCA trajectory and Hopf-fiber projection |
| `metadata.json` | Configuration, normalization audit, PCA loadings, DMD matrix, and claim boundaries |

A compact gist render and a richer public concept render are included below:

![MEAO Spinor Lab quick gist](docs/assets/quick_gist.png)

![MEAO Spinor Lab public concept demo](docs/assets/concept_demo_h2.png)

Both images are generated by public educational code. Neither is output from the complete private research solver.

---

## Scientific scope and claim boundaries

This section is not legal padding. It is the guardrail that keeps an exploratory model from wearing a lab coat three sizes too large.

### 1. Reduced-order, not ab-initio

The analytic molecular envelopes, nonlinear spinor field, and live MEAO/QIT quantities are phenomenological. Quantitative chemistry should be validated with PySCF, Psi4, ORCA, Gaussian, or another appropriate electronic-structure workflow.

### 2. MEAO-inspired proxies are not time-resolved MEAO measurements

The live PDE-side mutual-information, concurrence, log-negativity, and $F_{\mathrm{MEAO}}$-style channels are qualitative reduced-state proxies. A genuine MEAO analysis requires validated orbital reduced states and an orbital-rotation optimization.

### 3. The BLP channel is a density-field proxy

The project monitors density distinguishability and its derivative. It does not yet implement the canonical optimized trace-distance measure over reduced density operators.

### 4. Signed quadrature is not signed density

$$
Q_\rho=\rho\cos\phi
$$

is gauge-sensitive phase information. Negative values do not imply negative probability or electron density.

### 5. Display samples are not automatically physical modes

The full dashboard's 236 featured $S^2$ identities are stable display samples selected from the spin texture. They should not be described as 236 independently established physical eigenmodes without a separate modal analysis.

### 6. Hopf fibers are state-space fibers

The selected $U(1)$ loops are mathematically valid preimages in normalized spinor state space. Their linked appearance is not, by itself, an entanglement measurement, a braid invariant, or a physical-space Hopfion.

### 7. Clifford latent coordinates are not gamma-matrix dynamics

The Clifford algebra encodes state geometry. PCA produces learned statistical axes, and ridge-DMD fits a linear transition model. The latent ribbon is not a Dirac propagator, Lorentz transformation, or intrinsic manifold reconstruction.

### 8. Symmetry claims require explicit tests

The Bloch vector and several topology observables are invariant under a common local spinor phase. Other phase-sensitive quantities are not. The present framework does not claim a complete gauge theory, derived Lagrangian, Noether current, or demonstrated Lorentz covariance.

### 9. Optical pulse physics remains future work

The current optical controls are continuous Gaussian-family fields. Explicit pulsed one-photon and two-photon interaction models remain roadmap items.

---

## Validation philosophy

The private research framework and public demo follow the same general philosophy even though their numerical kernels differ.

### Numerical checks

- finiteness of every state and diagnostic;
- conservation of mass for conservative subflows;
- time-step refinement and second-order convergence studies;
- grid and box-size sensitivity;
- orthonormality of public bonding/antibonding channels;
- analytic sign checks for topological charge;
- separation of raw, derived, and render data;
- deterministic seeds and reproducible exports.

### Latent-model checks

- scaler, PCA basis, and DMD matrix fit on the early training block only;
- later states retained for held-out evaluation;
- comparison against persistence baselines;
- explicit export of PCA loadings and DMD matrix;
- no claim that a visually smooth ribbon is an intrinsic manifold.

### Topological checks

- each Hopf fiber remains on $S^3$;
- every point on a complete fiber maps to the same $S^2$ base;
- projection changes do not alter the underlying $S^3$ coordinates;
- future braid metrics must survive gauge, projection, selection, and null-model controls.

### Chemistry checks

- reduced diagnostics compared with independent quantum-chemistry anchors;
- proxy terminology retained until calibrated;
- molecule-shape differences distinguished from numerical normalization artifacts;
- no electron-count claims inferred from one effective spinor norm.

---

## Repository structure

```text
meao-spinor-lab/
├── README.md
├── LICENSE
├── NOTICE
├── CITATION.cff
├── requirements.txt
├── pyproject.toml
├── references.bib
├── COLLABORATION.md
├── CORE_CODE_ACCESS.md
├── CONTRIBUTING.md
├── ROADMAP.md
├── SECURITY.md
├── CODE_OF_CONDUCT.md
├── demo/
│   ├── quick_gist.py
│   ├── mini_spinor_lab.py
│   └── synthetic_diagnostics.py
├── tests/
│   ├── test_quick_gist.py
│   ├── test_mini_spinor_lab.py
│   └── test_synthetic_diagnostics.py
├── docs/
│   ├── equations.md
│   ├── references.md
│   ├── scientific_scope.md
│   └── assets/
│       ├── quick_gist.png
│       └── concept_demo_h2.png
├── social/
│   └── ...
└── .github/
    ├── ISSUE_TEMPLATE/
    ├── PULL_REQUEST_TEMPLATE.md
    └── workflows/ci.yml
```

---

## Roadmap

### Near term

- [ ] Add side-by-side H₂/CO public comparison scripts.
- [ ] Add unit-versus-mean-density box-padding tests.
- [ ] Add timestep and grid-convergence reports.
- [ ] Add richer null models for the memory channel.
- [ ] Add explicit PCA-loading and DMD-eigenvalue reports.
- [ ] Add a public fiber-marker tracking prototype.

### Medium term

- [ ] Add calibrated pulse envelopes.
- [ ] Add explicit one-photon and effective two-photon interaction models.
- [ ] Replace site-sharing proxies with reduced states derived from a quantum-chemistry backend.
- [ ] Implement validated MEAO orbital optimization over a selected active space.
- [ ] Expand the bond coordinate to normal modes $q_k(t)$.
- [ ] Compare conservative and explicitly dissipative open-system models.

### Long term

- [ ] Test braid-resolved correlates of independently computed information dynamics.
- [ ] Extend from 2D state-space visualization to a true 3D spatial spinor field.
- [ ] Evaluate physical-space preimages and converged Hopf invariants.
- [ ] Integrate experimental structured-light and spectroscopy calibration data.
- [ ] Build molecule-general interfaces for larger active spaces and reaction paths.

See [ROADMAP.md](ROADMAP.md) for the collaboration-facing roadmap.

---

## Collaboration

The complete research implementation is not distributed in this public repository. Focused private access may be considered for well-scoped collaborations.

Relevant areas include:

- quantum chemistry and orbital reduced states;
- nonlinear and open-system dynamics;
- structured-light and optical OAM experiments;
- MEAO validation;
- geometric and Clifford algebra;
- topological data analysis and braid theory;
- numerical PDE verification;
- reduced-order modeling and Koopman/DMD methods;
- scientific visualization and provenance engineering.

A useful collaboration proposal should identify:

1. the precise scientific question;
2. the validation target or falsifiable hypothesis;
3. the contributor's expertise and expected work product;
4. required code, data, compute, or experimental access;
5. reproducibility and publication expectations;
6. confidentiality, attribution, and authorship expectations;
7. a milestone that can be completed without releasing the entire private core.

Read:

- [COLLABORATION.md](COLLABORATION.md)
- [CORE_CODE_ACCESS.md](CORE_CODE_ACCESS.md)
- [CONTRIBUTING.md](CONTRIBUTING.md)

Then open the repository's **Collaboration Proposal** issue template.

---

## Citation and references

Use [CITATION.cff](CITATION.cff) to cite this public software scaffold.

The project draws on several scientific families rather than one monolithic formalism:

- maximally entangled atomic orbitals and orbital information analysis;
- Breuer-Laine-Piilo information-backflow diagnostics;
- PySCF and active-space reduced-density matrices;
- Hermite-Gaussian and Laguerre-Gaussian structured light;
- optical orbital angular momentum;
- geometric/Clifford algebra and spinors;
- Bloch and Hopf geometry;
- dynamic mode decomposition and reduced-order forecasting.

See:

- [docs/references.md](docs/references.md)
- [references.bib](references.bib)
- [docs/equations.md](docs/equations.md)
- [docs/scientific_scope.md](docs/scientific_scope.md)

When publishing results obtained with the private framework, cite the actual quantum-chemistry, numerical, information-theoretic, and topological methods used in that specific study. Do not cite the repository title as a substitute for those methods.

---

## Security and data hygiene

The public demo:

- performs no network requests;
- reads no API keys;
- does not launch a server or tunnel;
- writes only to the requested output directory;
- uses deterministic local computation;
- contains no private checkpoints or generated research datasets.

Never commit:

- `OPENAI_API_KEY`, `NGROK_AUTHTOKEN`, or other credentials;
- private solver files;
- checkpoint roots;
- collaborator-restricted data;
- unpublished experimental data;
- outputs whose release has not been approved.

See [SECURITY.md](SECURITY.md).

---

## License

Files actually published in this repository are licensed under the **Apache License 2.0**. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

The license applies only to materials included in this repository. It does not grant rights to unpublished source code, private datasets, private model outputs, credentials, research notebooks, or other materials that are not distributed here.

---

## Project identity

**Repository slug**

```text
meao-spinor-lab
```

**Display name**

```text
MEAO Spinor Lab
```

**Full title**

> **A Reduced-Order Spinor-Field Framework for H₂ and CO: MEAO-Inspired Bonding, Information Flow, and Hopf-Fiber Diagnostics**

**Suggested GitHub description**

```text
Reduced-order H₂/CO spinor-field dynamics with structured light, MEAO-inspired information diagnostics, Clifford embeddings, latent forecasting, and Hopf-fiber visualization.
```

**Suggested topics**

```text
spinor-dynamics
quantum-chemistry
meao
structured-light
optical-oam
hopf-fibration
clifford-algebra
non-markovianity
pyscf
reduced-order-model
dynamic-mode-decomposition
scientific-visualization
```
