# Programmable Molecular Spinor-Field Lab

A public-facing research scaffold for a private, physics-informed simulation of **molecularly shaped two-component spinor mean fields** in analytic H₂ and CO landscapes.

The private research code explores how structured-light and vortex controls co-evolve with reduced diagnostics for bonding, information flow, topology, and latent dynamics. The public repository intentionally exposes the scientific framing, equations, interfaces, a synthetic demonstration, and collaboration process without publishing the full solver, proprietary orchestration, unpublished parameter studies, or generated research data.

> **Scope in one sentence:** this is a reduced-order diagnostic and visualization framework, not a replacement for ab-initio electronic-structure theory.

## Research question

How do programmable optical and topological perturbations reorganize a molecularly shaped spinor field, and which reduced diagnostics move together during that evolution?

The private framework evolves

$$
\Psi(x,y,t)=
\begin{pmatrix}
\psi_\uparrow(x,y,t)\\
\psi_\downarrow(x,y,t)
\end{pmatrix},
\qquad
\rho=|\psi_\uparrow|^2+|\psi_\downarrow|^2,
$$

with a controlled Pauli/Gross-Pitaevskii-like nonlinear wave model. Analytic H₂ and CO orbitals provide molecular shape and polarity, while structured-light and programmable vortices provide perturbations.

## Current private implementation

- Analytic 2D H₂ and CO molecular envelopes and a Morse bond coordinate.
- Continuous Hermite-Gaussian and Laguerre-Gaussian structured-light controls, including orbital-angular-momentum index $\ell$.
- Vortex charge, core radius, depletion depth, lattice, orbit radius, angular velocity, and drive strength.
- Conservative Strang split-step evolution and explicit optional dissipative terms.
- Density, phase, and signed phase quadrature

  $$Q_\rho=\rho\cos\phi,$$

  where negative values indicate phase inversion, not negative probability density.
- Bond-order and transition proxies, plus explicitly labelled MEAO/QIT-inspired reduced-state diagnostics such as mutual information, concurrence, and log-negativity proxies.
- A BLP-inspired density-distinguishability current, not the canonical density-matrix BLP measure.
- Bloch texture, curvature, and density-masked topological-charge diagnostics:

  $$
  \mathbf n=\frac{\Psi^\dagger\boldsymbol\sigma\Psi}{\Psi^\dagger\Psi+\varepsilon},
  \qquad
  q=\frac{1}{4\pi}\mathbf n\cdot
  (\partial_x\mathbf n\times\partial_y\mathbf n).
  $$
- Clifford embeddings in $\mathrm{Cl}(3,0)$ or $\mathrm{Cl}(1,3)$, local Clifford descriptors, PCA, and held-out ridge-DMD diagnostics.
- A 236-identity S² display sampling and 14 selected representatives lifted to complete Hopf $U(1)$ preimage circles in S³, shown by stereographic projection.
- Optional PySCF RHF plus active-space CI anchors stored separately from the reduced PDE trajectory.

## Important claim boundaries

1. **MEAO means Maximally Entangled Atomic Orbitals.** The live PDE-side quantities are MEAO-inspired proxies unless a separately labelled quantum-chemistry anchor is used.
2. **The 236 beads are stable display samples, not 236 independently established physical field modes.**
3. **The 14 Hopf fibers are exact geometric preimages of selected S² display states.** Their visual linking or braid-like appearance is not, by itself, a measured molecular entanglement invariant.
4. **The Clifford latent ribbon is a derived PCA/ridge-DMD coordinate view.** It is not an intrinsic manifold reconstruction and is not evidence of Lorentz covariance.
5. **Common-phase $U(1)$ invariance holds for Bloch, mass, and topology observables.** Some Clifford phase features remain explicitly gauge-sensitive. The implementation does not claim a complete gauge theory or a demonstrated Noether current.
6. **Single-photon and two-photon pulsed interaction channels are roadmap items.** The current private build implements continuous Gaussian-family HG/LG controls.

See [Scientific Scope](docs/scientific_scope.md), [Equations](docs/equations.md), and [References](docs/references.md).

## Public demo

The demo generates deterministic **synthetic signals only**. It illustrates the public diagnostic vocabulary without reproducing the private PDE, controller, molecular parameters, selection policies, or rendering pipeline.

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python demo/synthetic_diagnostics.py --out demo_output
python -m unittest discover -s tests -v
```

Expected outputs:

- `demo_output/synthetic_diagnostics.csv`
- `demo_output/synthetic_diagnostics.png`
- `demo_output/metadata.json`

## Collaboration and core-code access

The complete research implementation is not distributed in this public repository. Researchers interested in validation, quantum chemistry, nonlinear dynamics, geometric algebra, topological visualization, or experimental calibration can open the **Collaboration Proposal** issue template.

Review criteria include:

- a concrete research question or validation target;
- relevant expertise and expected contribution;
- a reproducibility and publication plan;
- data, compute, confidentiality, and authorship expectations;
- agreement on which materials may be shared privately.

Read [COLLABORATION.md](COLLABORATION.md) and [CORE_CODE_ACCESS.md](CORE_CODE_ACCESS.md) before opening a proposal.

## Licensing boundary

Files actually published in this repository are licensed under the Apache License 2.0. The license does **not** grant rights to unpublished source code, private datasets, model outputs, credentials, or research materials that are not included here. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

## Citation

Use [CITATION.cff](CITATION.cff) for the public scaffold. Cite the underlying scientific methods separately using [docs/references.md](docs/references.md).
