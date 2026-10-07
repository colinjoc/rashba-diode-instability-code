"""Inert binary64 source Riccati pilot implementation; no admission or certificate.

Fourier convention: eta(x)=sum eta[n] exp(i*n*K*x). eta* reverses modes.
Both a and b use z=omega+i*delta, not conjugate z. U=(a-a0)/beta
and V=(b-a0)/beta are solved directly, including the essential beta*D*U².
No defaults contain material parameters. Importing this file executes no solve.
All arithmetic/error diagnostics are empirical, without directed rounding.
"""
import math
import resource
import time
import numpy as np


def _finite_real(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(name + ' must be a finite real scalar')
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(name + ' must be finite')
    return result


def _poly(eta):
    if not isinstance(eta, dict) or not eta:
        raise ValueError('eta must be a nonempty explicit harmonic dictionary')
    terms = {}
    for key, value in eta.items():
        n = int(key)
        if str(n) != str(key) and not isinstance(key, (int, np.integer)):
            raise ValueError('eta harmonic must be an integer')
        if isinstance(key, (bool, np.bool_)):
            raise ValueError('boolean harmonic')
        if isinstance(value, (list, tuple)):
            if len(value) != 2:
                raise ValueError('complex coefficient needs real/imag pair')
            c = complex(_finite_real(value[0], 'eta real'), _finite_real(value[1], 'eta imag'))
        else:
            c = complex(value)
        if not math.isfinite(c.real) or not math.isfinite(c.imag):
            raise ValueError('nonfinite eta coefficient')
        terms[n] = c
    degree = max(abs(n) for n in terms)
    values = np.zeros(2 * degree + 1, dtype=np.complex128)
    for n, c in terms.items():
        values[n + degree] = c
    return values


def _add(a, b):
    size = max(len(a), len(b))
    result = np.zeros(size, dtype=np.complex128)
    for p in (a, b):
        start = (size - len(p)) // 2
        result[start:start + len(p)] += p
    return result


def _multiply(a, b, guard):
    guard()
    return np.convolve(a, b)


def _norm(a):
    return float(np.sum(np.abs(a)))


def _evaluate(p, phase, guard):
    # One harmonic at a time bounds memory independently of Nx*mode count.
    result = np.zeros(len(phase), dtype=np.complex128)
    degree = len(p) // 2
    for j, c in enumerate(p):
        guard()
        if c != 0:
            result += c * np.exp(1j * (j - degree) * phase)
    return result


def _rhs(profile, eta, a0, D, beta, guard):
    star = np.conj(eta[::-1])
    forcing = _add(eta, -a0 * a0 * star)
    nonlinear = _add(2 * beta * a0 * _multiply(star, profile, guard),
                     beta * _multiply(_add(np.array([D], dtype=np.complex128), beta * star),
                                      _multiply(profile, profile, guard), guard))
    return _add(forcing, -nonlinear)


def _sector(eta, a0, Omega, D, beta, signed_uK, M, iterations, tolerance, guard):
    profile = np.zeros(2 * M + 1, dtype=np.complex128)
    modes = np.arange(-M, M + 1)
    denominator = 2 * Omega + 1j * modes * signed_uK
    if np.min(np.abs(denominator)) == 0:
        raise ValueError('unresolved Fourier resolvent pole')
    defect = math.inf
    for step in range(iterations):
        guard()
        rhs = _rhs(profile, eta, a0, D, beta, guard)
        degree = len(rhs) // 2
        nxt = rhs[degree - M:degree + M + 1] / denominator
        if not np.all(np.isfinite(nxt)):
            raise ValueError('nonfinite nonlinear iterate')
        defect = _norm(nxt - profile)
        profile = nxt
        if defect <= tolerance:
            break
    # Entire convolution support is retained here, including modes beyond M.
    residual = _add(denominator * profile, -_rhs(profile, eta, a0, D, beta, guard))
    high = residual.copy()
    rd = len(high) // 2
    high[rd - M:rd + M + 1] = 0
    return profile, {'iterations': step + 1, 'last_step_l1': defect,
                     'step_tolerance_reached': bool(defect <= tolerance),
                     'full_unprojected_residual_l1': _norm(residual),
                     'omitted_mode_residual_l1': _norm(high),
                     'residual_max_mode': rd,
                     'resolvent_denominator_min': float(np.min(np.abs(denominator)))}


def solve_node(omega, delta, D, beta, u, K, eta, *, M, iterations, tolerance, Nx, guard):
    """Return one node's mean Re(eta* (f-f0)/beta), including beta=0.

    guard is mandatory, callable with no arguments, and must raise on deadline,
    STOP or resource breach; it is called before every convolution/iteration and
    spatial harmonic. Caller owns aggregate scheduling/accounting and finite caps.
    Nx is periodic midpoint quadrature in phase. For K=0 the supplied eta still
    defines a formal phase profile; a physical uniform field requires eta[0].
    No frequency/angular/path integration, tail bound or admission is performed.
    """
    started, cpu = time.monotonic(), time.process_time()
    if not callable(guard):
        raise ValueError('mandatory external guard must be callable')
    guard()
    omega, delta, D, beta, u, K = [_finite_real(x, n) for x, n in
                                  zip((omega, delta, D, beta, u, K),
                                      ('omega', 'delta', 'D', 'beta', 'u', 'K'))]
    if omega <= 0:
        raise ValueError('positive Matsubara frequency required')
    for name, value in (('M', M), ('iterations', iterations), ('Nx', Nx)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(name + ' must be positive integer')
    tolerance = _finite_real(tolerance, 'tolerance')
    if tolerance < 0:
        raise ValueError('nonnegative tolerance required')
    E = _poly(eta)
    if len(E) // 2 > M or Nx <= 2 * M:
        raise ValueError('eta modes must fit M and Nx must exceed 2M')
    z = complex(omega, delta)
    Omega = complex(np.sqrt(z * z + D * D))
    if Omega.real < 0:
        Omega = -Omega
    if Omega.real <= 0 or z + Omega == 0:
        raise ValueError('causal background branch unresolved')
    a0 = D / (z + Omega)
    phase = 2 * np.pi * (np.arange(Nx) + 0.5) / Nx
    ex = _evaluate(E, phase, guard)
    exact_grazing = bool(u == 0 or K == 0)
    if exact_grazing:
        # Stable local quadratic root, no subtraction a-a0 or division by beta.
        physical_gap = D + beta * ex
        local_Omega = np.sqrt(z * z + np.abs(physical_gap) ** 2)
        local_Omega = np.where(local_Omega.real < 0, -local_Omega, local_Omega)
        Ua = (ex - a0 * a0 * np.conj(ex)) / (Omega + beta * a0 * np.conj(ex) + local_Omega)
        Vb = (np.conj(ex) - a0 * a0 * ex) / (Omega + beta * a0 * ex + local_Omega)
        ra = 2 * Omega * Ua + 2 * beta * a0 * np.conj(ex) * Ua + beta * (D + beta * np.conj(ex)) * Ua ** 2 - (ex - a0 * a0 * np.conj(ex))
        rb = 2 * Omega * Vb + 2 * beta * a0 * ex * Vb + beta * (D + beta * ex) * Vb ** 2 - (np.conj(ex) - a0 * a0 * ex)
        sa = {'local_algebra_residual_max': float(np.max(np.abs(ra))), 'iterations': 0}
        sb = {'local_algebra_residual_max': float(np.max(np.abs(rb))), 'iterations': 0}
        centre_U_l1 = centre_V_l1 = None
        branch = {'type': 'causal local algebra on each spatial node',
                  'local_ReOmega_min': float(np.min(local_Omega.real)),
                  'continuous_profile_quadrature_required': True}
    else:
        U, sa = _sector(E, a0, Omega, D, beta, u * K, M, iterations, tolerance, guard)
        V, sb = _sector(np.conj(E[::-1]), a0, Omega, D, beta, -u * K, M, iterations, tolerance, guard)
        Ua, Vb = _evaluate(U, phase, guard), _evaluate(V, phase, guard)
        centre_U_l1, centre_V_l1 = _norm(U), _norm(V)
        # A sufficient causal continuation ball diagnostic, uniformly from 0 to beta.
        enorm, A, B, L = _norm(E), abs(a0), abs(beta), 1 / (2 * Omega.real)
        R = 2 * L * enorm * (1 + A * A)
        rhs = L * (enorm * (1 + A * A) + 2 * B * A * enorm * R + B * (abs(D) + B * enorm) * R * R)
        lip = L * (2 * B * A * enorm + 2 * B * (abs(D) + B * enorm) * R)
        branch = {'type': 'empirical Fourier-l1 uniform causal contraction diagnostic',
                  'R': R, 'radius_rhs': rhs, 'contraction_upper_centre': lip,
                  'physical_unit_disk_upper_centre': A + B * R,
                  'passes_empirical_conditions': bool(R > 0 and rhs <= R and lip < 1 and A + B * R < 1 and max(centre_U_l1, centre_V_l1) <= R),
                  'rigorous_certificate': False}
    guard()
    q = 1 + a0 * a0
    second = q + beta * a0 * (Ua + Vb) + beta * beta * Ua * Vb
    denominator = q * second
    if np.min(np.abs(denominator)) == 0:
        raise ValueError('normalized force pole')
    df_over_beta = 2 * (Ua - a0 * a0 * Vb - beta * a0 * Ua * Vb) / denominator
    force_samples = np.real(np.conj(ex) * df_over_beta)
    if not np.all(np.isfinite(force_samples)):
        raise ValueError('nonfinite normalized force')
    force = float(np.mean(force_samples))
    paired_grid_difference = None
    if Nx % 2 == 0:
        paired_grid_difference = abs(float(np.mean(force_samples[::2])) - float(np.mean(force_samples[1::2])))
    physical_a, physical_b = a0 + beta * Ua, a0 + beta * Vb
    normalization_error = float(np.max(np.abs(second - (1 + physical_a * physical_b))))
    return {'scope': 'SUPPORTING_BINARY64_NODE_NOT_TARGET_CONFIRMATION',
            'certified': False, 'node_force': force, 'force_mean': force,
            'exact_grazing_local_algebra': exact_grazing,
            'background': {'Omega': [Omega.real, Omega.imag], 'a0': [a0.real, a0.imag],
                           'same_z_a_b': True},
            'solver_a': sa, 'solver_b': sb, 'branch': branch,
            'normalization': {'q_modulus': abs(q), 'second_denominator_min': float(np.min(np.abs(second))),
                              'full_denominator_min': float(np.min(np.abs(denominator))),
                              'normalization_identity_max_empirical': normalization_error,
                              'physical_a_abs_max': float(np.max(np.abs(physical_a))),
                              'physical_b_abs_max': float(np.max(np.abs(physical_b)))},
            'quadrature': {'Nx': int(Nx), 'rule': 'periodic phase midpoint',
                           'interlaced_mean_difference_empirical': paired_grid_difference,
                           'certified_error': None},
            'centre_U_l1': centre_U_l1, 'centre_V_l1': centre_V_l1,
            'uncertainty': {'binary64_rounding_bound': None, 'input_error_bound': None,
                            'full_residual_force_bound': None, 'tail_bound': None,
                            'interpretation': 'resolution/residual/normalization differences are empirical diagnostics; no certified interval'},
            'accounting': {'self_cpu_seconds': time.process_time() - cpu,
                           'wall_seconds': time.monotonic() - started,
                           'process_maxrss_kib': int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                           'scope': 'current process; RSS is lifetime high water, not per-node delta'}}
