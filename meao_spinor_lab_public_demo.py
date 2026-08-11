#!/usr/bin/env python3
"""MEAO Spinor Lab public concept demo.

This runnable example exposes the *architecture* of the research project without
reproducing the private solver, calibrated parameter studies, controller,
selection policies, dashboards, or unpublished datasets.

The demo intentionally uses a compact, educational model:

1. Build analytic two-center H2- or CO-like orbital envelopes.
2. Initialize a two-component Pauli spinor from an orthonormalized bonding mode.
3. Apply continuous Gaussian/Hermite-Gaussian/Laguerre-Gaussian perturbations.
4. Evolve the field with a simplified unitary Strang split-step method.
5. Extract bonding, MEAO-inspired, memory, Bloch, topology, Clifford, PCA, and
   ridge-DMD diagnostics.
6. Lift the final global spinor to a complete Hopf U(1) fiber in S3 and show a
   stereographic R3 projection.

Scientific boundary
-------------------
This is a reduced-order concept demonstration in dimensionless units. It is not
ab-initio electronic-structure theory, a time-resolved MEAO calculation, a
canonical BLP non-Markovianity measure, or a physical-space Hopfion solver.

Python 3.11+, NumPy, and Matplotlib are the only dependencies.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EPS = 1.0e-12


@dataclass(frozen=True)
class MoleculeSpec:
    """Illustrative, dimensionless two-center parameters.

    These values are deliberately generic and are not the private project's
    calibrated or production parameter set.
    """

    name: str
    atom_a: str
    atom_b: str
    zeta_a: float
    zeta_b: float
    coeff_a: float
    coeff_b: float
    equilibrium_separation: float
    polarity: float
    well_depth: float
    morse_a: float
    reduced_mass: float


MOLECULES: Mapping[str, MoleculeSpec] = {
    "H2": MoleculeSpec(
        name="H2",
        atom_a="H",
        atom_b="H",
        zeta_a=1.00,
        zeta_b=1.00,
        coeff_a=1.00,
        coeff_b=1.00,
        equilibrium_separation=1.45,
        polarity=0.00,
        well_depth=1.00,
        morse_a=1.15,
        reduced_mass=1.00,
    ),
    "CO": MoleculeSpec(
        name="CO",
        atom_a="C",
        atom_b="O",
        zeta_a=1.10,
        zeta_b=1.55,
        coeff_a=0.78,
        coeff_b=1.00,
        equilibrium_separation=1.85,
        polarity=0.22,
        well_depth=1.55,
        morse_a=1.30,
        reduced_mass=4.00,
    ),
}


@dataclass(frozen=True)
class DemoConfig:
    molecule: str = "H2"
    light: str = "lg"  # off | gaussian | hg | lg | mixed
    grid: int = 64
    length: float = 10.0
    dt: float = 0.004
    steps: int = 180
    seed: int = 13424

    # Simplified dimensionless PDE coefficients.
    kinetic: float = 0.075
    nonlinearity: float = 0.28
    transverse_spin_mix: float = 0.075

    # Structured-light controls.
    light_amplitude: float = 0.28
    light_waist: float = 2.8
    light_omega: float = 3.2
    lg_l: int = 1
    hg_m: int = 1
    hg_n: int = 0

    # Integer spin-vortex imprint.
    vortex: bool = False
    vortex_charge: int = 1
    vortex_core: float = 0.50

    # Normalization is explicit because molecular probability amplitudes and
    # mean-field order parameters use different conventions.
    normalization_mode: str = "unit"  # unit | mean-density
    mean_density: float = 0.05

    # Tiny reference perturbation for the density-distinguishability diagnostic.
    reference_phase_perturbation: float = 0.018

    # Reduced bond-coordinate plant. This is illustrative, not spectroscopy.
    bond_damping: float = 0.10
    bond_information_force: float = 0.020
    bond_light_force: float = 0.008

    # Latent model.
    train_fraction: float = 0.65
    ridge: float = 1.0e-4


@dataclass
class SimulationResult:
    config: DemoConfig
    molecule: MoleculeSpec
    x: np.ndarray
    y: np.ndarray
    dx: float
    diagnostics: dict[str, np.ndarray]
    chart_names: list[str]
    pca_components: np.ndarray
    explained_variance_ratio: np.ndarray
    dmd_matrix: np.ndarray
    final_state: np.ndarray
    final_reference_state: np.ndarray
    final_fields: dict[str, np.ndarray]
    final_bonding: np.ndarray
    final_antibonding: np.ndarray
    hopf_q4: np.ndarray
    hopf_r3: np.ndarray
    hopf_base: np.ndarray


def _validate_config(cfg: DemoConfig) -> None:
    if cfg.molecule.upper() not in MOLECULES:
        raise ValueError(f"molecule must be one of {sorted(MOLECULES)}")
    if cfg.light not in {"off", "gaussian", "hg", "lg", "mixed"}:
        raise ValueError("light must be off, gaussian, hg, lg, or mixed")
    if cfg.grid < 24:
        raise ValueError("grid must be at least 24")
    if cfg.length <= 0.0 or cfg.dt <= 0.0 or cfg.steps < 8:
        raise ValueError("length and dt must be positive; steps must be at least 8")
    if cfg.light_waist <= 0.0 or cfg.vortex_core <= 0.0:
        raise ValueError("light_waist and vortex_core must be positive")
    if cfg.normalization_mode not in {"unit", "mean-density"}:
        raise ValueError("normalization_mode must be unit or mean-density")
    if cfg.mean_density <= 0.0:
        raise ValueError("mean_density must be positive")
    if not 0.50 <= cfg.train_fraction <= 0.85:
        raise ValueError("train_fraction must be between 0.50 and 0.85")
    if cfg.ridge <= 0.0:
        raise ValueError("ridge must be positive")


def make_grid(n: int, length: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return a periodic square grid in dimensionless coordinates."""

    x = np.linspace(-0.5 * length, 0.5 * length, n, endpoint=False, dtype=np.float64)
    dx = float(length / n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    return x, X, Y, dx


def inner2d(a: np.ndarray, b: np.ndarray, dx: float) -> complex:
    return complex(np.sum(np.conj(a) * b, dtype=np.complex128) * dx * dx)


def normalize_scalar(phi: np.ndarray, dx: float) -> np.ndarray:
    norm = math.sqrt(max(float(np.sum(np.abs(phi) ** 2) * dx * dx), EPS))
    return np.asarray(phi, dtype=np.complex128) / norm


def spinor_mass(psi: np.ndarray, dx: float) -> float:
    rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
    return float(np.sum(rho, dtype=np.float64) * dx * dx)


def target_mass(cfg: DemoConfig) -> float:
    if cfg.normalization_mode == "unit":
        return 1.0
    return float(cfg.mean_density * cfg.length * cfg.length)


def normalize_spinor(psi: np.ndarray, dx: float, desired_mass: float) -> np.ndarray:
    mass = spinor_mass(psi, dx)
    if not np.isfinite(mass) or mass <= EPS:
        raise FloatingPointError("cannot normalize an empty or non-finite spinor")
    return np.asarray(psi, dtype=np.complex128) * math.sqrt(desired_mass / mass)


def atomic_orbital(
    X: np.ndarray,
    Y: np.ndarray,
    center_x: float,
    zeta: float,
    *,
    softening: float = 0.12,
) -> np.ndarray:
    radius = np.sqrt((X - center_x) ** 2 + Y**2 + softening**2)
    return np.exp(-zeta * radius)


def molecular_basis(
    X: np.ndarray,
    Y: np.ndarray,
    dx: float,
    spec: MoleculeSpec,
    separation: float,
) -> dict[str, Any]:
    """Build orthonormalized two-center channels and a soft molecular potential."""

    x_a = -0.5 * separation
    x_b = +0.5 * separation
    phi_a = normalize_scalar(atomic_orbital(X, Y, x_a, spec.zeta_a), dx)
    phi_b = normalize_scalar(atomic_orbital(X, Y, x_b, spec.zeta_b), dx)
    overlap = inner2d(phi_a, phi_b, dx)

    bonding = normalize_scalar(spec.coeff_a * phi_a + spec.coeff_b * phi_b, dx)
    antibonding_raw = spec.coeff_b * phi_a - spec.coeff_a * phi_b
    # Gram-Schmidt makes the public projection channels unambiguous even for CO.
    antibonding_raw = antibonding_raw - inner2d(bonding, antibonding_raw, dx) * bonding
    antibonding = normalize_scalar(antibonding_raw, dx)

    r_a = np.sqrt((X - x_a) ** 2 + Y**2 + 0.30**2)
    r_b = np.sqrt((X - x_b) ** 2 + Y**2 + 0.30**2)
    charge_a = 0.5 * (1.0 - spec.polarity)
    charge_b = 0.5 * (1.0 + spec.polarity)
    potential = -0.85 * (charge_a / r_a + charge_b / r_b)
    potential += 0.025 * spec.polarity * X

    return {
        "phi_a": phi_a,
        "phi_b": phi_b,
        "bonding": bonding,
        "antibonding": antibonding,
        "overlap": overlap,
        "potential": potential.astype(np.float64),
        "centers": (x_a, x_b),
    }


def morse_potential(separation: float, spec: MoleculeSpec) -> float:
    u = math.exp(-spec.morse_a * (separation - spec.equilibrium_separation))
    return float(spec.well_depth * (1.0 - u) ** 2 - spec.well_depth)


def morse_force(separation: float, spec: MoleculeSpec) -> float:
    u = math.exp(-spec.morse_a * (separation - spec.equilibrium_separation))
    d_v = 2.0 * spec.well_depth * spec.morse_a * u * (1.0 - u)
    return float(-d_v)


def imprint_integer_vortex(
    psi: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    charge: int,
    core: float,
) -> np.ndarray:
    radius = np.sqrt(X**2 + Y**2)
    angle = np.arctan2(Y, X)
    depletion = np.tanh(radius / max(core, EPS)) ** abs(int(charge))
    out = np.array(psi, dtype=np.complex128, copy=True)
    out[0] *= depletion * np.exp(1j * charge * angle)
    out[1] *= depletion * np.exp(-1j * charge * angle)
    return out


def initialize_spinor(
    X: np.ndarray,
    Y: np.ndarray,
    dx: float,
    basis: Mapping[str, Any],
    spec: MoleculeSpec,
    cfg: DemoConfig,
) -> tuple[np.ndarray, np.ndarray]:
    bonding = np.asarray(basis["bonding"], dtype=np.complex128)
    bias = 0.04 * spec.polarity
    rng = np.random.default_rng(cfg.seed)
    relative_phase = (
        0.26 * np.tanh(Y)
        + 0.08 * spec.polarity * X
        + 0.002 * rng.standard_normal(X.shape)
    )
    up = math.sqrt(0.5 + bias) * bonding * np.exp(0.08j * np.sin(Y))
    down = math.sqrt(0.5 - bias) * bonding * np.exp(1j * relative_phase)
    psi = np.stack([up, down], axis=0)
    if cfg.vortex:
        psi = imprint_integer_vortex(
            psi, X, Y, cfg.vortex_charge, cfg.vortex_core
        )
    psi = normalize_spinor(psi, dx, target_mass(cfg))

    # A nearby reference preparation for a BLP-inspired density-distance channel.
    perturb = cfg.reference_phase_perturbation * np.sin(1.3 * X - 0.7 * Y)
    reference = np.array(psi, copy=True)
    reference[0] *= np.exp(1j * perturb)
    reference[1] *= np.exp(-1j * perturb)
    reference = normalize_spinor(reference, dx, target_mass(cfg))
    return psi, reference


def hermite_polynomial(order: int, x: np.ndarray) -> np.ndarray:
    if order < 0:
        raise ValueError("Hermite order must be non-negative")
    if order == 0:
        return np.ones_like(x)
    if order == 1:
        return 2.0 * x
    h_nm2 = np.ones_like(x)
    h_nm1 = 2.0 * x
    for n in range(2, order + 1):
        h_n = 2.0 * x * h_nm1 - 2.0 * (n - 1) * h_nm2
        h_nm2, h_nm1 = h_nm1, h_n
    return h_nm1


def structured_light(
    X: np.ndarray,
    Y: np.ndarray,
    time_value: float,
    cfg: DemoConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return scalar drive, spin-mixing drive, and complex optical envelope."""

    if cfg.light == "off" or cfg.light_amplitude == 0.0:
        zeros = np.zeros_like(X, dtype=np.float64)
        return zeros, zeros, zeros.astype(np.complex128)

    waist = max(cfg.light_waist, EPS)
    radius2 = X**2 + Y**2
    gaussian = np.exp(-radius2 / waist**2)
    phase_t = cfg.light_omega * time_value

    xi = math.sqrt(2.0) * X / waist
    eta = math.sqrt(2.0) * Y / waist
    hg_envelope = (
        hermite_polynomial(cfg.hg_m, xi)
        * hermite_polynomial(cfg.hg_n, eta)
        * gaussian
    )
    hg = hg_envelope * np.exp(-1j * phase_t)

    radius = np.sqrt(radius2 + EPS)
    theta = np.arctan2(Y, X)
    lg_envelope = (math.sqrt(2.0) * radius / waist) ** abs(cfg.lg_l) * gaussian
    lg = lg_envelope * np.exp(1j * (cfg.lg_l * theta - phase_t))

    if cfg.light == "gaussian":
        optical = gaussian * np.exp(-1j * phase_t)
    elif cfg.light == "hg":
        optical = hg
    elif cfg.light == "lg":
        optical = lg
    else:
        optical = (hg + lg) / math.sqrt(2.0)

    scale = float(np.max(np.abs(optical)))
    if scale > EPS:
        optical = optical / scale
    scalar_drive = cfg.light_amplitude * np.real(optical)
    spin_drive = 0.38 * cfg.light_amplitude * np.imag(optical)
    return scalar_drive.astype(np.float64), spin_drive.astype(np.float64), optical


def kinetic_phase(shape: tuple[int, int], dx: float, dt: float, kinetic: float) -> np.ndarray:
    """Continuum spectral kinetic phase used only by this abbreviated demo."""

    kx = 2.0 * np.pi * np.fft.fftfreq(shape[0], d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(shape[1], d=dx)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    return np.exp(-0.5j * dt * kinetic * (KX**2 + KY**2))


def split_step(
    psi: np.ndarray,
    scalar_potential: np.ndarray,
    spin_drive: np.ndarray,
    phase_half: np.ndarray,
    cfg: DemoConfig,
) -> np.ndarray:
    """One second-order, norm-preserving Strang step for the public mini-model."""

    spectrum = np.fft.fft2(psi, axes=(-2, -1), norm="ortho")
    state = np.fft.ifft2(
        spectrum * phase_half[None, :, :], axes=(-2, -1), norm="ortho"
    )

    rho = np.abs(state[0]) ** 2 + np.abs(state[1]) ** 2
    local_phase = np.exp(-1j * cfg.dt * (scalar_potential + cfg.nonlinearity * rho))
    mixing = cfg.transverse_spin_mix + spin_drive
    cosine = np.cos(cfg.dt * mixing)
    sine = np.sin(cfg.dt * mixing)
    up, down = state[0], state[1]
    state = local_phase[None, :, :] * np.stack(
        [
            cosine * up - 1j * sine * down,
            cosine * down - 1j * sine * up,
        ],
        axis=0,
    )

    spectrum = np.fft.fft2(state, axes=(-2, -1), norm="ortho")
    state = np.fft.ifft2(
        spectrum * phase_half[None, :, :], axes=(-2, -1), norm="ortho"
    )
    if not np.all(np.isfinite(state)):
        raise FloatingPointError("split-step produced non-finite values")
    return np.asarray(state, dtype=np.complex128)


def central_gradient(field: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray]:
    gx = (np.roll(field, -1, axis=-2) - np.roll(field, 1, axis=-2)) / (2.0 * dx)
    gy = (np.roll(field, -1, axis=-1) - np.roll(field, 1, axis=-1)) / (2.0 * dx)
    return gx, gy


def local_spin(psi: np.ndarray) -> np.ndarray:
    up, down = psi[0], psi[1]
    rho = np.abs(up) ** 2 + np.abs(down) ** 2
    denominator = rho + EPS
    sx = 2.0 * np.real(np.conj(up) * down) / denominator
    sy = 2.0 * np.imag(np.conj(up) * down) / denominator
    sz = (np.abs(up) ** 2 - np.abs(down) ** 2) / denominator
    return np.stack([sx, sy, sz], axis=0)


def curvature_field(spin: np.ndarray, dx: float) -> np.ndarray:
    curvature = np.zeros_like(spin[0], dtype=np.float64)
    for component in spin:
        gx, gy = central_gradient(component, dx)
        curvature += np.real(gx * np.conj(gx) + gy * np.conj(gy))
    return curvature


def topological_charge_density(spin: np.ndarray, dx: float) -> np.ndarray:
    norm = np.sqrt(np.sum(spin * spin, axis=0) + EPS)
    n = spin / norm[None, :, :]
    dnx = np.empty_like(n)
    dny = np.empty_like(n)
    for component in range(3):
        dnx[component], dny[component] = central_gradient(n[component], dx)
    cross = np.cross(np.moveaxis(dnx, 0, -1), np.moveaxis(dny, 0, -1))
    chirality = np.sum(np.moveaxis(n, 0, -1) * cross, axis=-1)
    return chirality / (4.0 * np.pi)


def oam_expectation(psi: np.ndarray, X: np.ndarray, Y: np.ndarray, dx: float) -> float:
    mass = spinor_mass(psi, dx)
    total = 0.0
    for component in psi:
        d_dx, d_dy = central_gradient(component, dx)
        generator = X * d_dy - Y * d_dx
        total += float(
            np.sum(np.imag(np.conj(component) * generator), dtype=np.float64)
            * dx
            * dx
        )
    return total / max(mass, EPS)


def density_distinguishability(psi_a: np.ndarray, psi_b: np.ndarray, dx: float) -> float:
    rho_a = np.abs(psi_a[0]) ** 2 + np.abs(psi_a[1]) ** 2
    rho_b = np.abs(psi_b[0]) ** 2 + np.abs(psi_b[1]) ** 2
    mass_scale = 0.5 * (float(np.sum(rho_a)) + float(np.sum(rho_b))) * dx * dx
    return float(0.5 * np.sum(np.abs(rho_a - rho_b)) * dx * dx / max(mass_scale, EPS))


def binary_entropy(probability: float) -> float:
    p = float(np.clip(probability, EPS, 1.0 - EPS))
    return float(-(p * math.log(p) + (1.0 - p) * math.log(1.0 - p)))


def bond_and_information_proxies(
    psi: np.ndarray,
    basis: Mapping[str, Any],
    dx: float,
) -> dict[str, float]:
    bonding = np.asarray(basis["bonding"])
    antibonding = np.asarray(basis["antibonding"])
    phi_a = np.asarray(basis["phi_a"])
    phi_b = np.asarray(basis["phi_b"])

    n_bonding = sum(abs(inner2d(bonding, component, dx)) ** 2 for component in psi)
    n_antibonding = sum(abs(inner2d(antibonding, component, dx)) ** 2 for component in psi)
    channel_total = n_bonding + n_antibonding + EPS
    transition = float(n_antibonding / channel_total)
    bond_order = float((n_bonding - n_antibonding) / channel_total)

    n_a = sum(abs(inner2d(phi_a, component, dx)) ** 2 for component in psi)
    n_b = sum(abs(inner2d(phi_b, component, dx)) ** 2 for component in psi)
    p_a = float(n_a / (n_a + n_b + EPS))
    p_b = 1.0 - p_a
    overlap_weight = float(np.clip(abs(complex(basis["overlap"])), 0.0, 1.0))

    # These are deliberately labelled proxies. They are not orbital RDM measures.
    mi_proxy = float(binary_entropy(p_a) / math.log(2.0) * overlap_weight)
    concurrence_proxy = float(
        np.clip(2.0 * math.sqrt(max(p_a * p_b, 0.0)) * overlap_weight, 0.0, 1.0)
    )
    return {
        "bond_order_proxy": bond_order,
        "transition_proxy": transition,
        "meao_mutual_information_proxy": mi_proxy,
        "meao_concurrence_proxy": concurrence_proxy,
        "site_probability_a": p_a,
        "site_probability_b": p_b,
    }


# -----------------------------------------------------------------------------
# Compact Cl(3,0) engine for the public concept demo
# -----------------------------------------------------------------------------


def _popcount(value: int) -> int:
    return int(value.bit_count())


def _blade_gp(a: int, b: int) -> tuple[float, int]:
    """Geometric product of canonical Cl(3,0) basis blades stored as bit masks."""

    sign = 1.0
    for i in range(3):
        if (a >> i) & 1 and _popcount(b & ((1 << i) - 1)) % 2:
            sign *= -1.0
    # All three basis vectors square to +1 in Cl(3,0), so no metric sign is added.
    return sign, a ^ b


def clifford_gp(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.zeros(8, dtype=np.float64)
    for i in np.flatnonzero(np.abs(a) > 0.0):
        for j in np.flatnonzero(np.abs(b) > 0.0):
            sign, mask = _blade_gp(int(i), int(j))
            out[mask] += sign * float(a[i]) * float(b[j])
    return out


def clifford_reverse(multivector: np.ndarray) -> np.ndarray:
    signs = np.array(
        [(-1.0) ** (_popcount(mask) * (_popcount(mask) - 1) // 2) for mask in range(8)],
        dtype=np.float64,
    )
    return np.asarray(multivector, dtype=np.float64) * signs


def grade_norm(multivector: np.ndarray, grade: int) -> float:
    indices = [mask for mask in range(8) if _popcount(mask) == grade]
    return float(np.linalg.norm(np.asarray(multivector)[indices]))


def global_spinor(psi: np.ndarray) -> tuple[complex, complex]:
    rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
    weights = rho / (float(np.sum(rho)) + EPS)
    alpha = complex(np.sum(weights * psi[0]))
    beta = complex(np.sum(weights * psi[1]))
    norm = math.sqrt(abs(alpha) ** 2 + abs(beta) ** 2)
    if norm <= 1.0e-10:
        alpha = complex(np.sum(psi[0]))
        beta = complex(np.sum(psi[1]))
        norm = math.sqrt(abs(alpha) ** 2 + abs(beta) ** 2)
    if norm <= 1.0e-10:
        return 1.0 + 0.0j, 0.0 + 0.0j
    return alpha / norm, beta / norm


def clifford_descriptors(psi: np.ndarray) -> dict[str, float]:
    """Encode a Pauli spinor in the even subalgebra and form a left-ideal density."""

    alpha, beta = global_spinor(psi)
    # Phi = a + b e12 + c e23 + d e31.
    # Canonical bit-mask basis stores e13, so e31 = -e13.
    phi = np.zeros(8, dtype=np.float64)
    phi[0] = alpha.real
    phi[3] = alpha.imag  # e12
    phi[6] = beta.real  # e23
    phi[5] = -beta.imag  # e31 = -e13

    idempotent = np.zeros(8, dtype=np.float64)
    idempotent[0] = 0.5
    idempotent[4] = 0.5  # (1 + e3)/2
    phi_left = clifford_gp(phi, idempotent)
    rho_c = clifford_gp(phi_left, clifford_reverse(phi_left))

    scalar = float(rho_c[0])
    vector_norm = grade_norm(rho_c, 1)
    bivector_norm = grade_norm(rho_c, 2)
    pseudoscalar_norm = grade_norm(rho_c, 3)
    total = float(np.linalg.norm(rho_c) + EPS)
    purity_like = float((scalar**2 + vector_norm**2) / total**2)
    rotor_angle = float(2.0 * math.atan2(grade_norm(phi, 2), abs(phi[0]) + EPS))
    return {
        "global_alpha_real": alpha.real,
        "global_alpha_imag": alpha.imag,
        "global_beta_real": beta.real,
        "global_beta_imag": beta.imag,
        "clifford_rho_scalar": scalar,
        "clifford_rho_vector_norm": vector_norm,
        "clifford_rho_bivector_norm": bivector_norm,
        "clifford_rho_pseudoscalar_norm": pseudoscalar_norm,
        "clifford_left_ideal_purity_like": purity_like,
        "clifford_rotor_angle": rotor_angle,
    }


def field_diagnostics(
    psi: np.ndarray,
    reference: np.ndarray,
    basis: Mapping[str, Any],
    X: np.ndarray,
    Y: np.ndarray,
    dx: float,
    separation: float,
    separation_velocity: float,
    spec: MoleculeSpec,
    light_scalar: np.ndarray,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    rho = np.abs(psi[0]) ** 2 + np.abs(psi[1]) ** 2
    phase = np.angle(psi[0] + psi[1])
    signed_quadrature = rho * np.cos(phase)
    spin = local_spin(psi)
    density_floor = max(1.0e-8, 1.0e-3 * float(np.max(rho)))
    density_mask = rho >= density_floor
    stencil_valid = np.array(density_mask, copy=True)
    for axis in (-2, -1):
        stencil_valid &= np.roll(density_mask, 1, axis=axis)
        stencil_valid &= np.roll(density_mask, -1, axis=axis)
    curvature = np.where(stencil_valid, curvature_field(spin, dx), 0.0)
    q_density = np.where(
        stencil_valid, topological_charge_density(spin, dx), 0.0
    )
    phase_masked = np.where(density_mask, phase, np.nan)

    scalar = bond_and_information_proxies(psi, basis, dx)
    scalar.update(clifford_descriptors(psi))
    mass = spinor_mass(psi, dx)
    mean_sz = float(np.sum((np.abs(psi[0]) ** 2 - np.abs(psi[1]) ** 2)) * dx * dx / max(mass, EPS))
    qfi_z = float(max(0.0, 4.0 * (1.0 - mean_sz**2)))
    scalar.update(
        {
            "mass": mass,
            "density_distinguishability_proxy": density_distinguishability(
                psi, reference, dx
            ),
            "qfi_sigma_z": qfi_z,
            "curvature_mean": float(np.mean(curvature)),
            "curvature_p95": float(np.percentile(curvature, 95.0)),
            "topological_charge": float(np.sum(q_density) * dx * dx),
            "oam_expectation": oam_expectation(psi, X, Y, dx),
            "bond_separation": float(separation),
            "bond_velocity": float(separation_velocity),
            "morse_potential": morse_potential(separation, spec),
            "light_rms": float(np.sqrt(np.mean(light_scalar**2))),
        }
    )
    arrays = {
        "density": rho,
        "phase": phase,
        "phase_masked": phase_masked,
        "signed_quadrature": signed_quadrature,
        "spin": spin,
        "curvature": curvature,
        "topological_charge_density": q_density,
    }
    return scalar, arrays


def _append_row(storage: dict[str, list[float]], row: Mapping[str, float]) -> None:
    for key, value in row.items():
        storage.setdefault(key, []).append(float(value))


def fit_pca_ridge_dmd(
    diagnostics: Mapping[str, np.ndarray],
    cfg: DemoConfig,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Train-only PCA followed by a regularized linear latent transition model."""

    preferred = [
        "bond_order_proxy",
        "transition_proxy",
        "meao_mutual_information_proxy",
        "meao_concurrence_proxy",
        "density_distinguishability_proxy",
        "memory_current_proxy",
        "qfi_sigma_z",
        "curvature_mean",
        "curvature_p95",
        "topological_charge",
        "oam_expectation",
        "bond_separation",
        "morse_potential",
        "light_rms",
        "clifford_rho_scalar",
        "clifford_rho_vector_norm",
        "clifford_rho_bivector_norm",
        "clifford_left_ideal_purity_like",
        "clifford_rotor_angle",
    ]
    names = [name for name in preferred if name in diagnostics]
    chart = np.column_stack([np.asarray(diagnostics[name], dtype=np.float64) for name in names])
    n = chart.shape[0]
    train_end = int(np.clip(math.floor(cfg.train_fraction * n), 6, n - 2))
    train = np.nan_to_num(chart[:train_end], nan=0.0, posinf=0.0, neginf=0.0)
    mean = np.mean(train, axis=0, keepdims=True)
    scale = np.std(train, axis=0, keepdims=True)
    scale = np.where(scale > 1.0e-9, scale, 1.0)
    standardized = np.nan_to_num((chart - mean) / scale)

    _, singular_values, vt = np.linalg.svd(standardized[:train_end], full_matrices=False)
    k = min(3, vt.shape[0])
    components = np.asarray(vt[:k], dtype=np.float64)
    for component in components:
        pivot = int(np.argmax(np.abs(component)))
        if component[pivot] < 0.0:
            component *= -1.0
    latent = standardized @ components.T
    if k < 3:
        components = np.pad(components, ((0, 3 - k), (0, 0)))
        latent = np.pad(latent, ((0, 0), (0, 3 - k)))

    x_train = latent[: train_end - 1]
    y_train = latent[1:train_end]
    regularizer = cfg.ridge * np.eye(3)
    dmd = np.linalg.solve(
        x_train.T @ x_train + regularizer,
        x_train.T @ y_train,
    )
    one_step = np.full_like(latent, np.nan)
    one_step[1:] = latent[:-1] @ dmd
    error = np.full(n, np.nan, dtype=np.float64)
    error[1:] = np.linalg.norm(latent[1:] - one_step[1:], axis=1)

    variance = singular_values**2
    explained = np.zeros(3, dtype=np.float64)
    if variance.size and float(np.sum(variance)) > EPS:
        explained[:k] = variance[:k] / float(np.sum(variance))
    return names, components, explained, dmd, latent, error


def hopf_map(z1: complex, z2: complex) -> np.ndarray:
    return np.asarray(
        [
            2.0 * np.real(np.conj(z1) * z2),
            2.0 * np.imag(np.conj(z1) * z2),
            abs(z1) ** 2 - abs(z2) ** 2,
        ],
        dtype=np.float64,
    )


def _plane_rotation(i: int, j: int, angle: float) -> np.ndarray:
    matrix = np.eye(4, dtype=np.float64)
    c = math.cos(angle)
    s = math.sin(angle)
    matrix[i, i] = c
    matrix[j, j] = c
    matrix[i, j] = -s
    matrix[j, i] = s
    return matrix


def hopf_fiber_from_spinor(
    z1: complex,
    z2: complex,
    samples: int = 240,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    norm = math.sqrt(abs(z1) ** 2 + abs(z2) ** 2)
    if norm <= EPS:
        raise ValueError("spinor must be nonzero")
    z1 /= norm
    z2 /= norm
    alpha = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    phase = np.exp(1j * alpha)
    fiber_z1 = phase * z1
    fiber_z2 = phase * z2
    q4 = np.column_stack(
        [fiber_z1.real, fiber_z1.imag, fiber_z2.real, fiber_z2.imag]
    )

    # One fixed orthogonal view avoids choosing the raw pole through the fiber.
    rotation = _plane_rotation(0, 3, 0.37) @ _plane_rotation(1, 2, -0.41)
    rotated = q4 @ rotation.T
    denominator = np.maximum(1.0 - rotated[:, 3], 1.0e-8)
    raw = rotated[:, :3] / denominator[:, None]
    raw_norm = np.linalg.norm(raw, axis=1, keepdims=True)
    compact_radius = 4.0
    projected = compact_radius * raw / (compact_radius + raw_norm)
    return q4, projected, hopf_map(z1, z2)


def run_simulation(cfg: DemoConfig) -> SimulationResult:
    _validate_config(cfg)
    spec = MOLECULES[cfg.molecule.upper()]
    x, X, Y, dx = make_grid(cfg.grid, cfg.length)
    separation = float(spec.equilibrium_separation)
    separation_velocity = 0.0
    basis = molecular_basis(X, Y, dx, spec, separation)
    psi, reference = initialize_spinor(X, Y, dx, basis, spec, cfg)
    phase_half = kinetic_phase((cfg.grid, cfg.grid), dx, cfg.dt, cfg.kinetic)

    scalar_history: dict[str, list[float]] = {"time": []}
    last_arrays: dict[str, np.ndarray] = {}

    for step in range(cfg.steps + 1):
        time_value = step * cfg.dt
        light_scalar, _, _ = structured_light(X, Y, time_value, cfg)
        basis = molecular_basis(X, Y, dx, spec, separation)
        scalar, arrays = field_diagnostics(
            psi,
            reference,
            basis,
            X,
            Y,
            dx,
            separation,
            separation_velocity,
            spec,
            light_scalar,
        )
        scalar_history["time"].append(float(time_value))
        _append_row(scalar_history, scalar)
        last_arrays = arrays
        if step == cfg.steps:
            break

        # A compact reduced bond plant. Information and optical forces are
        # illustrative feedback channels, not calibrated molecular forces.
        information_force = cfg.bond_information_force * (
            scalar["transition_proxy"] - 0.10
        )
        optical_force = cfg.bond_light_force * scalar["light_rms"]
        force = (
            morse_force(separation, spec)
            + information_force
            + optical_force
            - cfg.bond_damping * separation_velocity
        )
        separation_velocity += cfg.dt * force / max(spec.reduced_mass, EPS)
        separation += cfg.dt * separation_velocity
        separation = float(
            np.clip(
                separation,
                0.62 * spec.equilibrium_separation,
                1.65 * spec.equilibrium_separation,
            )
        )

        midpoint = time_value + 0.5 * cfg.dt
        midpoint_basis = molecular_basis(X, Y, dx, spec, separation)
        scalar_drive, spin_drive, _ = structured_light(X, Y, midpoint, cfg)
        potential = np.asarray(midpoint_basis["potential"]) + scalar_drive
        psi = split_step(psi, potential, spin_drive, phase_half, cfg)
        reference = split_step(reference, potential, spin_drive, phase_half, cfg)

    diagnostics = {
        key: np.asarray(values, dtype=np.float64)
        for key, values in scalar_history.items()
    }
    distinguishability = diagnostics["density_distinguishability_proxy"]
    diagnostics["memory_current_proxy"] = np.gradient(
        distinguishability, diagnostics["time"], edge_order=1
    )

    chart_names, components, explained, dmd, latent, dmd_error = fit_pca_ridge_dmd(
        diagnostics, cfg
    )
    diagnostics["latent_1"] = latent[:, 0]
    diagnostics["latent_2"] = latent[:, 1]
    diagnostics["latent_3"] = latent[:, 2]
    diagnostics["ridge_dmd_one_step_error"] = dmd_error

    alpha, beta = global_spinor(psi)
    q4, r3, base = hopf_fiber_from_spinor(alpha, beta)

    return SimulationResult(
        config=cfg,
        molecule=spec,
        x=x,
        y=x.copy(),
        dx=dx,
        diagnostics=diagnostics,
        chart_names=chart_names,
        pca_components=components,
        explained_variance_ratio=explained,
        dmd_matrix=dmd,
        final_state=psi,
        final_reference_state=reference,
        final_fields=last_arrays,
        final_bonding=np.asarray(basis["bonding"]),
        final_antibonding=np.asarray(basis["antibonding"]),
        hopf_q4=q4,
        hopf_r3=r3,
        hopf_base=base,
    )


def write_diagnostics_csv(path: Path, diagnostics: Mapping[str, np.ndarray]) -> None:
    names = list(diagnostics)
    lengths = {len(np.asarray(diagnostics[name])) for name in names}
    if len(lengths) != 1:
        raise ValueError("diagnostic columns must have one common length")
    count = lengths.pop()
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for row in range(count):
            writer.writerow(
                [
                    "" if not np.isfinite(float(diagnostics[name][row])) else f"{float(diagnostics[name][row]):.12g}"
                    for name in names
                ]
            )


def write_hopf_csv(path: Path, q4: np.ndarray, r3: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "q0", "q1", "q2", "q3", "x", "y", "z"])
        for index, (four, three) in enumerate(zip(q4, r3)):
            writer.writerow([index, *[f"{float(v):.12g}" for v in four], *[f"{float(v):.12g}" for v in three]])


def write_field_plot(path: Path, result: SimulationResult) -> None:
    fields = result.final_fields
    extent = [result.x[0], result.x[-1], result.y[0], result.y[-1]]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)

    images = [
        ("density ρ", fields["density"], "viridis"),
        ("phase arg(ψ↑+ψ↓)", fields["phase_masked"], "twilight"),
        ("signed quadrature Qρ=ρ cosφ", fields["signed_quadrature"], "coolwarm"),
        ("spin-texture curvature |∇n|²", fields["curvature"], "magma"),
        ("topological-charge density q", fields["topological_charge_density"], "coolwarm"),
    ]
    for axis, (title, values, cmap) in zip(axes.flat[:5], images):
        image = axis.imshow(values.T, origin="lower", extent=extent, cmap=cmap, aspect="equal")
        axis.set_title(title)
        axis.set_xlabel("x (reduced units)")
        axis.set_ylabel("y (reduced units)")
        fig.colorbar(image, ax=axis, shrink=0.82)

    time_values = result.diagnostics["time"]
    axes[1, 2].plot(time_values, result.diagnostics["bond_order_proxy"], label="bond order")
    axes[1, 2].plot(
        time_values,
        result.diagnostics["meao_mutual_information_proxy"],
        label="MEAO MI proxy",
    )
    axes[1, 2].plot(
        time_values,
        result.diagnostics["meao_concurrence_proxy"],
        label="concurrence proxy",
    )
    axes[1, 2].set_title("reduced bonding/information diagnostics")
    axes[1, 2].set_xlabel("reduced time")
    axes[1, 2].set_ylabel("diagnostic value")
    axes[1, 2].legend(loc="best", fontsize=8)

    fig.suptitle(
        f"MEAO Spinor Lab concept demo: {result.molecule.name}, {result.config.light.upper()} drive\n"
        "Illustrative reduced-order output, not ab-initio electronic structure",
        fontsize=14,
    )
    fig.savefig(path, dpi=160)
    plt.close(fig)


def write_latent_hopf_plot(path: Path, result: SimulationResult) -> None:
    fig = plt.figure(figsize=(14, 6), constrained_layout=True)
    latent_axis = fig.add_subplot(1, 2, 1, projection="3d")
    hopf_axis = fig.add_subplot(1, 2, 2, projection="3d")

    z1 = result.diagnostics["latent_1"]
    z2 = result.diagnostics["latent_2"]
    z3 = result.diagnostics["latent_3"]
    time_values = result.diagnostics["time"]
    latent_axis.plot(z1, z2, z3, linewidth=1.5)
    latent_axis.scatter([z1[0]], [z2[0]], [z3[0]], marker="o", label="start")
    latent_axis.scatter([z1[-1]], [z2[-1]], [z3[-1]], marker="x", label="end")
    latent_axis.set_title("train-only PCA diagnostic trajectory")
    latent_axis.set_xlabel("PC1")
    latent_axis.set_ylabel("PC2")
    latent_axis.set_zlabel("PC3")
    latent_axis.legend(loc="best")

    hopf_axis.plot(
        result.hopf_r3[:, 0],
        result.hopf_r3[:, 1],
        result.hopf_r3[:, 2],
        linewidth=1.6,
    )
    hopf_axis.set_title("one complete Hopf U(1) fiber\nS³ stereographically projected to R³")
    hopf_axis.set_xlabel("x")
    hopf_axis.set_ylabel("y")
    hopf_axis.set_zlabel("z")

    fig.suptitle(
        "Derived diagnostic coordinates and state-space topology\n"
        "PCA axes are not gamma matrices; the fiber is not a physical-space Hopfion"
    )
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def write_outputs(result: SimulationResult, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    write_diagnostics_csv(out / "diagnostics.csv", result.diagnostics)
    write_hopf_csv(out / "hopf_fiber.csv", result.hopf_q4, result.hopf_r3)
    write_field_plot(out / "field_and_diagnostics.png", result)
    write_latent_hopf_plot(out / "latent_and_hopf.png", result)

    final_mass = float(result.diagnostics["mass"][-1])
    initial_mass = float(result.diagnostics["mass"][0])
    metadata = {
        "schema": "meao-spinor-lab-public-concept-v1",
        "scientific_status": "educational reduced-order concept demo",
        "private_solver_reproduced": False,
        "ab_initio": False,
        "config": asdict(result.config),
        "molecule": asdict(result.molecule),
        "normalization": {
            "mode": result.config.normalization_mode,
            "target_mass": target_mass(result.config),
            "initial_mass": initial_mass,
            "final_mass": final_mass,
            "relative_mass_drift": (final_mass - initial_mass) / max(abs(initial_mass), EPS),
        },
        "latent": {
            "chart_names": result.chart_names,
            "pca_components": result.pca_components,
            "explained_variance_ratio": result.explained_variance_ratio,
            "ridge_dmd_matrix": result.dmd_matrix,
            "interpretation": "derived linear diagnostic coordinates; not an intrinsic manifold or gamma-matrix dynamics",
        },
        "hopf": {
            "base_s2": result.hopf_base,
            "fiber_samples": int(len(result.hopf_q4)),
            "interpretation": "complete U(1) preimage circle in spinor state space; not a physical-space Hopfion",
        },
        "claim_boundaries": [
            "MEAO quantities are illustrative proxies, not orbital-RDM entanglement measurements.",
            "The memory current differentiates a density-field distance, not the canonical optimized BLP measure.",
            "The Clifford map is a state encoding; PCA and ridge-DMD produce the latent trajectory.",
            "The public split-step kernel is simplified and is not the private finite-difference-symbol solver.",
            "Pulsed one-photon and two-photon channels are not implemented in this demo.",
        ],
    }
    (out / "metadata.json").write_text(
        json.dumps(_json_safe(metadata), indent=2) + "\n", encoding="utf-8"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", choices=sorted(MOLECULES), default="H2")
    parser.add_argument(
        "--light", choices=["off", "gaussian", "hg", "lg", "mixed"], default="lg"
    )
    parser.add_argument("--grid", type=int, default=64)
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=0.004)
    parser.add_argument("--steps", type=int, default=180)
    parser.add_argument("--seed", type=int, default=13424)
    parser.add_argument("--light-amplitude", type=float, default=0.28)
    parser.add_argument("--light-waist", type=float, default=2.8)
    parser.add_argument("--lg-l", type=int, default=1)
    parser.add_argument("--hg", nargs=2, type=int, metavar=("M", "N"), default=(1, 0))
    vortex_group = parser.add_mutually_exclusive_group()
    vortex_group.add_argument("--vortex", action="store_true", help="imprint an integer spin vortex")
    vortex_group.add_argument("--no-vortex", action="store_true", help="explicitly disable the vortex")
    parser.add_argument("--vortex-charge", type=int, default=1)
    parser.add_argument(
        "--normalization-mode", choices=["unit", "mean-density"], default="unit"
    )
    parser.add_argument("--mean-density", type=float, default=0.05)
    parser.add_argument("--out", type=Path, default=Path("demo_output/mini_spinor_lab"))
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    cfg = DemoConfig(
        molecule=args.molecule,
        light=args.light,
        grid=args.grid,
        length=args.length,
        dt=args.dt,
        steps=args.steps,
        seed=args.seed,
        light_amplitude=args.light_amplitude,
        light_waist=args.light_waist,
        lg_l=args.lg_l,
        hg_m=int(args.hg[0]),
        hg_n=int(args.hg[1]),
        vortex=bool(args.vortex and not args.no_vortex),
        vortex_charge=args.vortex_charge,
        normalization_mode=args.normalization_mode,
        mean_density=args.mean_density,
    )
    result = run_simulation(cfg)
    write_outputs(result, args.out)
    print(f"Wrote MEAO Spinor Lab concept outputs to {args.out.resolve()}")
    print(
        "Final mass drift: "
        f"{(result.diagnostics['mass'][-1] - result.diagnostics['mass'][0]) / max(abs(result.diagnostics['mass'][0]), EPS):.3e}"
    )
    print(f"Final Hopf base: {result.hopf_base}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
