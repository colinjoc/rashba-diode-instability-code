"""Source-derived continuous Riccati precursor, never physical admission.

This module uses no primary bubble implementation. It evaluates the finite
positive-Matsubara/angular sum of the source force, including its counterterm.
It does not certify continuum, input, or rounding errors. Physical callers must
provide their reviewed native guard; the command line runs synthetic controls
only. Mandatory nonlinear C1/C2 confirmation remains separate.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
from dataclasses import dataclass
from typing import Callable

Matrix = tuple[tuple[complex, complex], tuple[complex, complex]]


def _finite_real(value: float, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def uniform_background(z: complex, gap: float) -> tuple[complex, complex]:
    """Both sectors have a0=b0=D/(z+Omega), not b0=conjugate(a0)."""
    z = complex(z)
    gap = _finite_real(gap, "gap")
    if z.real <= 0 or not all(map(math.isfinite, (z.real, z.imag))):
        raise ValueError("positive finite Matsubara real part required")
    if gap < 0:
        raise ValueError("gap must be nonnegative")
    omega = cmath.sqrt(z * z + gap * gap)
    if omega.real < 0:
        omega = -omega
    if omega.real <= 0:
        raise ArithmeticError("causal root unresolved")
    return omega, gap / (z + omega)


def harmonic_response(z: complex, gap: float, projection: float, sign: int) -> dict:
    """Coefficients for real A/P fields at the signed Fourier harmonic.

    projection = kappa dot (cos(theta), sin(theta)). The same complex z is
    retained in the two oppositely directed Riccati sectors.
    """
    if isinstance(sign, bool) or sign not in (-1, 1):
        raise ValueError("harmonic sign must be -1 or +1")
    projection = _finite_real(projection, "wavevector projection")
    omega, a = uniform_background(z, gap)
    t = a * a
    dp, dm = 2 * omega + 1j * sign * projection, 2 * omega - 1j * sign * projection
    da_a, da_p = (1 - t) / dp, 1j * (1 + t) / dp
    db_a, db_p = (1 - t) / dm, -1j * (1 + t) / dm
    norm = (1 + t) ** 2
    result = dict(omega=omega, a0=a, dp=dp, dm=dm, da_a=da_a,
                  da_p=da_p, db_a=db_a, db_p=db_p,
                  f_a=2 * (da_a - t * db_a) / norm,
                  f_p=2 * (da_p - t * db_p) / norm,
                  fd_a=2 * (db_a - t * da_a) / norm,
                  fd_p=2 * (db_p - t * da_p) / norm)
    if any(not math.isfinite(v.real) or not math.isfinite(v.imag) for v in result.values()):
        raise ArithmeticError("nonfinite harmonic response")
    return result


def force_matrix(z: complex, gap: float, projection: float) -> tuple[Matrix, Matrix]:
    """Return raw B and H: Re spatial_mean(eta* df) = v-dagger H v.

    Physical conjugation reverses harmonic indices. Taking the Hermitian part
    defines the real observable; it is not an independent Hermiticity/Ward test.
    """
    p = harmonic_response(z, gap, projection, 1)
    m = harmonic_response(z, gap, projection, -1)
    b = ((p["f_a"] + m["f_a"], p["f_p"] - 1j * m["f_a"]),
         (m["f_p"] - 1j * p["f_a"], -1j * (p["f_p"] + m["f_p"])))
    h = tuple(tuple((b[i][j] + b[j][i].conjugate()) / 2 for j in range(2)) for i in range(2))
    return b, h


def contract(matrix: Matrix, witness: tuple[complex, complex]) -> float:
    if len(witness) != 2:
        raise ValueError("amplitude/phase witness must have two components")
    v = tuple(complex(x) for x in witness)
    if any(not math.isfinite(x.real) or not math.isfinite(x.imag) for x in v):
        raise ValueError("finite witness required")
    result = sum(v[i].conjugate() * matrix[i][j] * v[j] for i in range(2) for j in range(2))
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ArithmeticError("nonfinite contraction")
    return result.real


@dataclass(frozen=True)
class SourceParameters:
    temperature: float
    field: float
    gradient: float
    gap: float
    kappa: tuple[float, float]
    dos_asymmetry: float = 0.25

    def validate(self) -> None:
        values = [self.temperature, self.field, self.gradient, self.gap, self.dos_asymmetry]
        if len(self.kappa) != 2:
            raise ValueError("two-dimensional kappa required")
        for v in [*values, *self.kappa]:
            _finite_real(v, "source parameter")
        if self.temperature <= 0 or self.gap < 0 or abs(self.dos_asymmetry) >= 1:
            raise ValueError("positive temperature, nonnegative gap and positive DOS weights required")


def finite_grid_diagnostic(
    source: SourceParameters,
    *,
    angular_nodes: int,
    frequencies: int,
    guard: Callable[[], None],
    witness: tuple[complex, complex] | None = None,
    max_evaluations: int = 2_000_000,
) -> dict:
    """Finite full-domain sum; guard must enforce caller authority/resources.

    M = 2 ln(T) I + pi T sum_n,lambda w_lambda <2/omega I-H>.
    Moving the helicity-independent counterterm inside w is valid because
    sum_lambda w_lambda=2, and avoids divergent cancellation between helicities.
    The supplied grid is a numerical diagnostic, not a continuum certificate.
    """
    if not callable(guard):
        raise ValueError("reviewed caller guard required")
    source.validate()
    for value in (angular_nodes, frequencies, max_evaluations):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("positive finite integer grid and allocation required")
    if angular_nodes % 4:
        raise ValueError("inversion-closed angular rule including exact grazing required")
    if 2 * angular_nodes * frequencies > max_evaluations:
        raise ValueError("finite evaluation allocation exceeded before numerical work")
    if witness is not None:
        contract(((0j, 0j), (0j, 0j)), witness)
    guard()
    directions = []
    for j in range(angular_nodes):
        guard()
        # Preserve the source rule theta=2*pi*j/N with exact axial/grazing nodes.
        if j % (angular_nodes // 4) == 0:
            c, s = ((1., 0.), (0., 1.), (-1., 0.), (0., -1.))[j // (angular_nodes // 4)]
        else:
            theta = 2 * math.pi * j / angular_nodes
            c, s = math.cos(theta), math.sin(theta)
        directions.append((c, source.kappa[0] * c + source.kappa[1] * s))
    rows = []
    raw_rows=[]
    tail_terms = [[], [], []]
    for n in range(frequencies):
        guard()
        omega = (2 * n + 1) * math.pi * source.temperature
        terms = [[], [], [], []]
        raw_terms=[[] for _ in range(4)]
        for helicity in (-1, 1):
            weight = 1 + helicity * source.dos_asymmetry
            for c, projection in directions:
                guard()
                delta = (source.gradient / 2 - helicity * source.field) * c
                z=omega+1j*delta
                _, h = force_matrix(z, source.gap, projection)
                r=harmonic_response(z,source.gap,projection,1)
                # Genuine source A/P force derivative from both anomalous sectors,
                # before any Hermitian projection. Angular inversion supplies
                # the physical pair; retain this independent reduction.
                g=((r['f_a']+r['fd_a'],r['f_p']+r['fd_p']),
                   ((r['f_a']-r['fd_a'])/1j,(r['f_p']-r['fd_p'])/1j))
                for i in range(2):
                    for j in range(2):
                        raw_terms[2*i+j].append(weight*((2/omega if i==j else 0)-g[i][j]))
                terms[0].append(weight * (2 / omega - h[0][0].real))
                terms[1].append(weight * (2 / omega - h[1][1].real))
                terms[2].append(-weight * h[0][1].real)
                terms[3].append(-weight * h[0][1].imag)
                if n == 0:
                    common = 2 * delta * delta + projection * projection / 2
                    tail_terms[0].append(weight * (common + 3 * source.gap ** 2))
                    tail_terms[1].append(weight * (common + source.gap ** 2))
                    tail_terms[2].append(weight * 2 * delta * projection)
        row = tuple(math.fsum(t) / angular_nodes for t in terms)
        if not all(map(math.isfinite, row)):
            raise ArithmeticError("nonfinite frequency reduction")
        rows.append(row)
        raw_rows.append(tuple(complex(math.fsum(x.real for x in t),math.fsum(x.imag for x in t))/angular_nodes for t in raw_terms))
    prefactor = math.pi * source.temperature
    aa = 2 * math.log(source.temperature) + prefactor * math.fsum(r[0] for r in rows)
    pp = 2 * math.log(source.temperature) + prefactor * math.fsum(r[1] for r in rows)
    ap = prefactor * complex(math.fsum(r[2] for r in rows), math.fsum(r[3] for r in rows))
    matrix = ((complex(aa), ap), (ap.conjugate(), complex(pp)))
    tail = ((math.fsum(tail_terms[0]) / angular_nodes, 1j * math.fsum(tail_terms[2]) / angular_nodes),
            (-1j * math.fsum(tail_terms[2]) / angular_nodes, math.fsum(tail_terms[1]) / angular_nodes))
    raw_matrix=tuple(tuple((2*math.log(source.temperature) if i==j else 0)+prefactor*complex(math.fsum(r[2*i+j].real for r in raw_rows),math.fsum(r[2*i+j].imag for r in raw_rows)) for j in range(2)) for i in range(2))
    raw_Hermiticity=max(abs(raw_matrix[i][j]-raw_matrix[j][i].conjugate()) for i in range(2) for j in range(2))
    raw_mapping=max(abs(raw_matrix[i][j]-matrix[i][j]) for i in range(2) for j in range(2))
    guard()
    return {
        'raw_source_AP_matrix':[[[x.real,x.imag] for x in row] for row in raw_matrix],
        'raw_source_Hermiticity_defect':raw_Hermiticity,
        'raw_source_to_force_mapping_discrepancy':raw_mapping,
        'raw_source_scope':'Both anomalous sectors, inversion-paired A/P force derivative before Hermitian projection; numerical check not complete error proof',
        "schema": "rashba-continuous-linear-finite-grid-v1",
        "scope": "Finite angular/frequency diagnostic only; no physical sign certification",
        "grid": {"angular_nodes": angular_nodes, "frequencies": frequencies,
                 "evaluations": 2 * angular_nodes * frequencies},
        "matrix": [[[x.real, x.imag] for x in row] for row in matrix],
        "finite_grid_witness_coefficient": None if witness is None else contract(matrix, witness),
        "leading_omega_minus3_tail_matrix": [[[complex(x).real, complex(x).imag] for x in row] for row in tail],
        "tail_semantics": "Multiply leading matrix by pi*T*sum_omitted omega^-3; NOT applied here. Higher orders/remainder unresolved.",
        "unresolved_errors": ["angular", "frequency_tail_remainder", "roundoff", "endpoint_input_coverage", "common_source"],
        "mandatory_followthrough": "Independent source/error/design/admission review and complete nonlinear C1/C2 confirmation",
        "scientific_admission": False, "qualification_pass": False,
        "parent_answered": False, "target_release": False,
    }


def synthetic_controls() -> dict:
    """Closed-form normal/uniform identities at synthetic inputs only."""
    source = SourceParameters(.7, 0., 0., 1.2, (0., 0.))
    result = finite_grid_diagnostic(source, angular_nodes=8, frequencies=12, guard=lambda: None)
    omegas = [(2 * n + 1) * math.pi * source.temperature for n in range(12)]
    expected = [2 * math.log(.7) + 4 * math.pi * .7 * math.fsum(
        1 / w - (w * w / (w * w + 1.2 ** 2) ** 1.5 if i == 0 else 1 / math.sqrt(w * w + 1.2 ** 2))
        for w in omegas) for i in range(2)]
    errors = [abs(result["matrix"][i][i][0] - expected[i]) for i in range(2)]
    return {"scope": "Synthetic uniform BCS finite-sum controls only",
            "absolute_errors": errors, "all_pass": all(x < 1e-12 for x in errors),
            "physical_execution": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-controls", action="store_true", required=True)
    parser.parse_args()
    result = synthetic_controls()
    print(json.dumps(result, allow_nan=False, sort_keys=True))
    if not result["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
