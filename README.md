

## Future experiment: braid-resolved prediction of information dynamics

A central future question is whether temporal braid structure contains information
about subsequent changes in independently evaluated molecular-information
diagnostics.

Define the change over a time lag $\tau$ as

$$
\Delta_{\tau} I(t)
=
I(t+\tau)-I(t).
$$

A first predictive model is

$$
\Delta_{\tau} I(t)
=
\beta_0
+
\sum_{j=1}^{p}\beta_j T_j(t)
+
\varepsilon(t+\tau).
$$

Here:

- $I(t)$ is an independently evaluated information quantity;
- $T_j(t)$ is the $j$th braid-derived feature measured at time $t$;
- $\tau$ is the prediction lag;
- $\beta_0$ is the fitted intercept;
- $\beta_j$ is the fitted contribution of braid feature $T_j$;
- $\varepsilon(t+\tau)$ represents unexplained variation.

Possible information quantities include mutual information, concurrence,
logarithmic negativity, or independently calculated reduced-density-matrix
diagnostics. Candidate braid features may include crossing statistics, braid-word
complexity, phase holonomy, motif persistence, or other explicitly defined
topological descriptors.

This equation represents a proposed future experiment. It is not a current claim
that braid formation measures entanglement.

### Evidence requirements

Credible evidence for a braid-information relationship will require:

1. **Gauge robustness:** the result must survive permitted common-phase changes.
2. **Projection robustness:** the conclusion must not depend on one stereographic
   projection, camera orientation, or viewing axis.
3. **Fiber-selection bootstrapping:** the analysis must be repeated across many
   subsets of selected fibers.
4. **Null-model testing:** randomized fiber labels, phase scrambling, and temporal
   shuffling must not reproduce the observed result.
5. **Held-out perturbation regimes:** predictive performance must generalize to
   laser, vortex, or dissipative regimes excluded during model fitting.
6. **Independent RDM validation:** braid features must be compared with separately
   generated reduced-density-matrix quantities rather than only with diagnostics
   constructed from the same reduced spinor trajectory.
7. **Claim separation:** correlation, out-of-sample prediction, and causal influence
   must be reported as distinct levels of evidence.

A suitable future-work heading is:

> **Toward Braid-Resolved Topological Correlates of Molecular Information Dynamics**

---

## What the public demo implements

The repository provides two public code paths:

- `demo/quick_gist.py`, a compact architectural sketch for readers who want the
  central idea quickly;
- `demo/mini_spinor_lab.py`, a more defensive public concept implementation with
  explicit normalization modes, a compact geometric-product engine, richer
  exports, metadata, and tests.

Neither file reproduces the complete private implementation.

### Included

- H₂-like and CO-like analytic two-center envelopes;
- orthonormalized bonding and antibonding channels;
- a two-component complex spinor field;
- unit-norm and mean-density normalization modes;
- Gaussian, Hermite-Gaussian, Laguerre-Gaussian, mixed, and disabled optical modes;
- integer-winding spin-vortex imprinting;
- a simplified unitary Strang split-step propagator;
- a reduced Morse bond coordinate;
- density, phase, signed-quadrature, and field plots;
- bond-order and transition proxies;
- explicitly labeled MEAO-inspired mutual-information and concurrence proxies;
- BLP-inspired density distinguishability and memory current;
- QFI, OAM, Bloch curvature, and topological-charge diagnostics;
- a compact $\mathrm{Cl}(3,0)$ geometric-product engine;
- minimal-left-ideal and Clifford-density descriptors;
- train-only PCA and ridge-DMD;
- one complete Hopf $U(1)$ fiber in $S^3$, with a stereographic projection into
  $R^3$;
- deterministic CSV, JSON, PNG, and test outputs.

### Intentionally omitted

- the complete private finite-difference-symbol kinetic and spin-orbit-coupling
  kernel;
- the production H₂/CO parameter set;
- the complete MPC controller and action search;
- private checkpoint, cache, provenance, and dashboard orchestration;
- private rendering and selected-fiber ranking logic;
- complete local-patch and spectral feature tensors;
- generated research data;
- credentials and network services;
- calibrated one-photon or two-photon pulse physics;
- a validated MEAO orbital-optimization backend;
- a physical-space three-dimensional Hopfion model.

The public values use generic reduced units. They are not a disguised copy of the
private parameterization.
