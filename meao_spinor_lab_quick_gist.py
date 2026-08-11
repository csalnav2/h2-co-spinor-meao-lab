#!/usr/bin/env python3
"""Compact public gist of the MEAO Spinor Lab architecture.

This script is intentionally small and generic. It is not the private solver,
an ab-initio calculation, or a time-resolved MEAO measurement. It demonstrates:

analytic H2/CO-like shapes -> two-component spinor -> HG/LG drive ->
unitary split-step evolution -> reduced diagnostics -> Clifford coordinates ->
PCA/ridge-DMD -> one Hopf U(1) fiber.

Dependencies: NumPy and Matplotlib only.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EPS = 1e-12
MOLECULES = {
    "H2": dict(zeta=(1.0, 1.0), coeff=(1.0, 1.0), R=1.45, polarity=0.0),
    "CO": dict(zeta=(1.1, 1.55), coeff=(0.78, 1.0), R=1.85, polarity=0.22),
}


def grid(n: int, length: float):
    x = np.linspace(-length / 2, length / 2, n, endpoint=False)
    dx = length / n
    X, Y = np.meshgrid(x, x, indexing="ij")
    return x, X, Y, dx


def inner(a, b, dx):
    return np.sum(np.conj(a) * b) * dx * dx


def normalize(phi, dx):
    return phi / np.sqrt(np.sum(np.abs(phi) ** 2) * dx * dx + EPS)


def molecular_model(X, Y, dx, name):
    p = MOLECULES[name]
    R = p["R"]
    xa, xb = -R / 2, R / 2
    za, zb = p["zeta"]
    ca, cb = p["coeff"]
    a = normalize(np.exp(-za * np.sqrt((X - xa) ** 2 + Y**2 + 0.12**2)), dx)
    b = normalize(np.exp(-zb * np.sqrt((X - xb) ** 2 + Y**2 + 0.12**2)), dx)
    bond = normalize(ca * a + cb * b, dx)
    anti = cb * a - ca * b
    anti = normalize(anti - inner(bond, anti, dx) * bond, dx)
    ra = np.sqrt((X - xa) ** 2 + Y**2 + 0.30**2)
    rb = np.sqrt((X - xb) ** 2 + Y**2 + 0.30**2)
    pol = p["polarity"]
    V = -0.85 * (0.5 * (1 - pol) / ra + 0.5 * (1 + pol) / rb)
    V += 0.025 * pol * X
    return a, b, bond, anti, V


def initial_spinor(bond, X, Y, dx, seed, vortex):
    rng = np.random.default_rng(seed)
    rel = 0.24 * np.tanh(Y) + 0.002 * rng.standard_normal(X.shape)
    psi = np.stack(
        [bond * np.exp(0.06j * np.sin(Y)), bond * np.exp(1j * rel)], axis=0
    ) / np.sqrt(2.0)
    if vortex:
        theta = np.arctan2(Y, X)
        core = np.tanh(np.sqrt(X**2 + Y**2) / 0.50)
        psi[0] *= core * np.exp(1j * theta)
        psi[1] *= core * np.exp(-1j * theta)
    return psi / np.sqrt(np.sum(np.abs(psi) ** 2) * dx * dx + EPS)


def light_field(X, Y, t, mode, amplitude, waist, ell):
    if mode == "off":
        z = np.zeros_like(X)
        return z, z
    r2 = X**2 + Y**2
    g = np.exp(-r2 / waist**2)
    hg = (np.sqrt(2) * X / waist) * g * np.exp(-3.2j * t)
    r = np.sqrt(r2 + EPS)
    lg = (np.sqrt(2) * r / waist) ** abs(ell) * g
    lg = lg * np.exp(1j * (ell * np.arctan2(Y, X) - 3.2 * t))
    optical = {"gaussian": g * np.exp(-3.2j * t), "hg": hg, "lg": lg}.get(
        mode, (hg + lg) / np.sqrt(2)
    )
    optical /= np.max(np.abs(optical)) + EPS
    return amplitude * optical.real, 0.38 * amplitude * optical.imag


def kinetic_half_phase(n, dx, dt, kappa):
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    return np.exp(-0.5j * dt * kappa * (KX**2 + KY**2))


def step(psi, V, spin_drive, half_phase, dt, g, omega0):
    p = np.fft.fft2(psi, axes=(-2, -1), norm="ortho")
    psi = np.fft.ifft2(p * half_phase[None], axes=(-2, -1), norm="ortho")
    rho = np.sum(np.abs(psi) ** 2, axis=0)
    phase = np.exp(-1j * dt * (V + g * rho))
    angle = dt * (omega0 + spin_drive)
    c, s = np.cos(angle), np.sin(angle)
    up, dn = psi
    psi = phase[None] * np.stack([c * up - 1j * s * dn, c * dn - 1j * s * up])
    p = np.fft.fft2(psi, axes=(-2, -1), norm="ortho")
    return np.fft.ifft2(p * half_phase[None], axes=(-2, -1), norm="ortho")


def bloch_and_topology(psi, dx):
    up, dn = psi
    rho = np.abs(up) ** 2 + np.abs(dn) ** 2
    n = np.stack(
        [
            2 * np.real(np.conj(up) * dn) / (rho + EPS),
            2 * np.imag(np.conj(up) * dn) / (rho + EPS),
            (np.abs(up) ** 2 - np.abs(dn) ** 2) / (rho + EPS),
        ]
    )
    grad = lambda f, axis: (np.roll(f, -1, axis) - np.roll(f, 1, axis)) / (2 * dx)
    dnx = np.stack([grad(c, 0) for c in n])
    dny = np.stack([grad(c, 1) for c in n])
    curvature = np.sum(dnx**2 + dny**2, axis=0)
    cross = np.cross(np.moveaxis(dnx, 0, -1), np.moveaxis(dny, 0, -1))
    q = np.sum(np.moveaxis(n, 0, -1) * cross, axis=-1) / (4 * np.pi)
    mask = rho > 1e-3 * np.max(rho)
    valid = mask & np.roll(mask, 1, 0) & np.roll(mask, -1, 0)
    valid &= np.roll(mask, 1, 1) & np.roll(mask, -1, 1)
    return rho, n, np.where(valid, curvature, 0.0), np.where(valid, q, 0.0)


def global_spinor_and_clifford(psi):
    rho = np.sum(np.abs(psi) ** 2, axis=0)
    w = rho / (np.sum(rho) + EPS)
    alpha, beta = np.sum(w * psi[0]), np.sum(w * psi[1])
    norm = np.sqrt(abs(alpha) ** 2 + abs(beta) ** 2)
    if norm < 1e-10:
        alpha, beta, norm = 1.0 + 0.0j, 0.0 + 0.0j, 1.0
    alpha, beta = alpha / norm, beta / norm
    # Coordinates of Phi = a + b e12 + c e23 + d e31 in Cl+(3,0).
    clifford = np.array([alpha.real, alpha.imag, beta.real, beta.imag])
    rotor_angle = 2 * np.arctan2(np.linalg.norm(clifford[1:]), abs(clifford[0]) + EPS)
    return alpha, beta, clifford, rotor_angle


def diagnostics(psi, ref, a, b, bond, anti, dx):
    nb = sum(abs(inner(bond, c, dx)) ** 2 for c in psi)
    na = sum(abs(inner(anti, c, dx)) ** 2 for c in psi)
    bond_order = (nb - na) / (nb + na + EPS)
    pa = sum(abs(inner(a, c, dx)) ** 2 for c in psi)
    pb = sum(abs(inner(b, c, dx)) ** 2 for c in psi)
    pa /= pa + pb + EPS
    overlap = np.clip(abs(inner(a, b, dx)), 0, 1)
    entropy = -(pa * np.log(pa + EPS) + (1 - pa) * np.log(1 - pa + EPS)) / np.log(2)
    mi_proxy = entropy * overlap
    concurrence_proxy = 2 * np.sqrt(max(pa * (1 - pa), 0)) * overlap
    rho, n, curvature, q = bloch_and_topology(psi, dx)
    rho_ref = np.sum(np.abs(ref) ** 2, axis=0)
    distance = 0.5 * np.sum(np.abs(rho - rho_ref)) * dx * dx
    alpha, beta, clifford, rotor = global_spinor_and_clifford(psi)
    return {
        "bond_order_proxy": float(bond_order),
        "meao_mi_proxy": float(mi_proxy),
        "concurrence_proxy": float(concurrence_proxy),
        "density_distance_proxy": float(distance),
        "curvature_mean": float(np.mean(curvature)),
        "topological_charge": float(np.sum(q) * dx * dx),
        "clifford_a": float(clifford[0]),
        "clifford_b": float(clifford[1]),
        "clifford_c": float(clifford[2]),
        "clifford_d": float(clifford[3]),
        "rotor_angle": float(rotor),
    }, (rho, n, curvature, q), (alpha, beta)


def latent_model(rows):
    names = [k for k in rows[0] if k != "time"]
    X = np.array([[row[k] for k in names] for row in rows])
    split = max(6, min(len(X) - 2, int(0.65 * len(X))))
    mu, sd = X[:split].mean(0), X[:split].std(0)
    Z = (X - mu) / np.where(sd > 1e-9, sd, 1)
    _, _, vt = np.linalg.svd(Z[:split], full_matrices=False)
    W = vt[:3]
    latent = Z @ W.T
    A = np.linalg.solve(
        latent[: split - 1].T @ latent[: split - 1] + 1e-4 * np.eye(3),
        latent[: split - 1].T @ latent[1:split],
    )
    return names, latent, A


def hopf_fiber(alpha, beta, samples=180):
    phase = np.exp(1j * np.linspace(0, 2 * np.pi, samples, endpoint=False))
    z1, z2 = phase * alpha, phase * beta
    q4 = np.c_[z1.real, z1.imag, z2.real, z2.imag]
    # Fixed S3 rotation followed by an ordinary stereographic projection.
    c, s = np.cos(0.37), np.sin(0.37)
    rotated = q4.copy()
    rotated[:, 0] = c * q4[:, 0] - s * q4[:, 3]
    rotated[:, 3] = s * q4[:, 0] + c * q4[:, 3]
    den = np.maximum(1 - rotated[:, 3], 1e-8)
    r3 = rotated[:, :3] / den[:, None]
    base = np.array(
        [2 * np.real(np.conj(alpha) * beta), 2 * np.imag(np.conj(alpha) * beta), abs(alpha) ** 2 - abs(beta) ** 2]
    )
    return q4, r3, base


def run(args):
    if args.grid < 24 or args.steps < 8 or args.length <= 0 or args.dt <= 0:
        raise ValueError("grid>=24, steps>=8, length>0, and dt>0 are required")
    x, X, Y, dx = grid(args.grid, args.length)
    a, b, bond, anti, Vmol = molecular_model(X, Y, dx, args.molecule)
    psi = initial_spinor(bond, X, Y, dx, args.seed, args.vortex)
    ref = psi.copy()
    ref[0] *= np.exp(0.018j * np.sin(1.3 * X - 0.7 * Y))
    ref[1] *= np.exp(-0.018j * np.sin(1.3 * X - 0.7 * Y))
    half = kinetic_half_phase(args.grid, dx, args.dt, 0.075)
    rows, final_fields, global_ab = [], None, None
    for k in range(args.steps + 1):
        t = k * args.dt
        row, final_fields, global_ab = diagnostics(psi, ref, a, b, bond, anti, dx)
        row["time"] = t
        rows.append(row)
        if k == args.steps:
            break
        scalar, spin = light_field(X, Y, t + args.dt / 2, args.light, 0.28, 2.8, args.ell)
        psi = step(psi, Vmol + scalar, spin, half, args.dt, 0.28, 0.075)
        ref = step(ref, Vmol + scalar, spin, half, args.dt, 0.28, 0.075)
    memory = np.gradient([r["density_distance_proxy"] for r in rows], [r["time"] for r in rows])
    for row, value in zip(rows, memory):
        row["memory_current_proxy"] = float(value)
    names, latent, A = latent_model(rows)
    q4, r3, base = hopf_fiber(*global_ab)
    return x, rows, final_fields, latent, A, names, q4, r3, base


def write_outputs(out, x, rows, fields, latent, A, names, q4, r3, base, args):
    out.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0])
    with (out / "diagnostics.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)
    np.savetxt(out / "hopf_fiber.csv", np.c_[q4, r3], delimiter=",", header="q0,q1,q2,q3,x,y,z", comments="")
    rho, _, curvature, q = fields
    fig, ax = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    extent = [x[0], x[-1], x[0], x[-1]]
    for axy, title, data, cmap in [
        (ax[0, 0], "density ρ", rho, "viridis"),
        (ax[0, 1], "curvature |∇n|²", curvature, "magma"),
        (ax[0, 2], "topological density q", q, "coolwarm"),
    ]:
        im = axy.imshow(data.T, origin="lower", extent=extent, cmap=cmap)
        axy.set_title(title)
        fig.colorbar(im, ax=axy, shrink=0.8)
    t = np.array([r["time"] for r in rows])
    ax[1, 0].plot(t, [r["bond_order_proxy"] for r in rows], label="bond order")
    ax[1, 0].plot(t, [r["meao_mi_proxy"] for r in rows], label="MEAO MI proxy")
    ax[1, 0].plot(t, [r["concurrence_proxy"] for r in rows], label="concurrence proxy")
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].set_title("reduced diagnostics")
    ax[1, 1].plot(latent[:, 0], latent[:, 1])
    ax[1, 1].set_title("PCA diagnostic trajectory")
    ax[1, 1].set_xlabel("PC1")
    ax[1, 1].set_ylabel("PC2")
    ax[1, 2].remove()
    hopf_ax = fig.add_subplot(2, 3, 6, projection="3d")
    hopf_ax.plot(r3[:, 0], r3[:, 1], r3[:, 2])
    hopf_ax.set_title("Hopf U(1) fiber")
    fig.suptitle(f"MEAO Spinor Lab quick gist: {args.molecule}, {args.light.upper()} drive\nnot ab-initio output")
    fig.savefig(out / "quick_gist.png", dpi=150)
    plt.close(fig)
    meta = {
        "status": "abbreviated educational demo",
        "private_solver_reproduced": False,
        "molecule": args.molecule,
        "light": args.light,
        "feature_names": names,
        "ridge_dmd_matrix": A.tolist(),
        "hopf_base_s2": base.tolist(),
        "claim_boundary": "proxies and state-space visualization only",
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--molecule", choices=["H2", "CO"], default="H2")
    p.add_argument("--light", choices=["off", "gaussian", "hg", "lg", "mixed"], default="lg")
    p.add_argument("--ell", type=int, default=1)
    p.add_argument("--vortex", action="store_true")
    p.add_argument("--grid", type=int, default=48)
    p.add_argument("--length", type=float, default=10.0)
    p.add_argument("--dt", type=float, default=0.004)
    p.add_argument("--steps", type=int, default=100)
    p.add_argument("--seed", type=int, default=13424)
    p.add_argument("--out", type=Path, default=Path("demo_output/quick_gist"))
    args = p.parse_args()
    result = run(args)
    write_outputs(args.out, *result, args)
    print(f"Wrote abbreviated demo to {args.out.resolve()}")


if __name__ == "__main__":
    main()
