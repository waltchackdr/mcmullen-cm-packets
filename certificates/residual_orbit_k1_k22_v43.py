#!/usr/bin/env python3
"""Exact exclusion certificate for the residual imprimitive pair k={1,22}.

The human-readable completeness proof is in
``docs/RESIDUAL_1_22_EXCLUSION.md``.  This certificate reuses only the exact
field, trace-lattice, determinant, and polygon primitives already regression-
gated in ``residual_orbit_k17_k38_v42``.  Every companion-specific quantity
and every archimedean search bound is rederived here using rational interval
arithmetic.  No binary floating-point value participates in an assertion.
"""

from __future__ import annotations

from fractions import Fraction as Q
from math import floor

import residual_orbit_k17_k38_v42 as core


def outward_cent(value: Q) -> Q:
    """Return the smallest hundredth strictly greater than ``value``."""
    return Q(floor(100 * value) + 1, 100)


def derive_bounds(y2: tuple[Q, ...]):
    """Derive outward rational entry bounds from the v42 polygon constants.

    For each real embedding, ``m`` is a certified lower bound for the minimum
    imaginary part in the target quadrilateral and ``U`` is a certified upper
    bound for |z|.  The standard formulas

        |cz+d| <= sqrt(y/m),       |az+b| <= U sqrt(y/m)

    imply the displayed bounds for a,b,c,d after separating real and
    imaginary parts.  All square roots are enclosed by rational intervals.
    """
    den = []
    num = []
    a_bound = []
    b_bound = []
    c_bound = []
    d_bound = []
    for i in range(core.DEGREE):
        y_lo, y_hi = core.sqrt_interval(
            core.interval_eval(y2, core.ROOT_INTERVALS[i]))
        pi_lo, pi_hi = core.sqrt_interval(
            core.interval_eval(core.pi, core.ROOT_INTERVALS[i]))
        assert y_lo > 0 and pi_lo > 0

        den_actual = core.sqrt_interval(
            (y_hi / core.M_LOWER[i], y_hi / core.M_LOWER[i]))[1]
        c_actual = core.sqrt_interval(
            (1 / (core.M_LOWER[i] * y_lo),
             1 / (core.M_LOWER[i] * y_lo)))[1]
        num_actual = core.U_UPPER[i] * den_actual
        a_actual = num_actual / y_lo
        d_actual = den_actual + c_actual * pi_hi
        b_actual = num_actual + a_actual * pi_hi

        den.append(outward_cent(den_actual))
        num.append(outward_cent(num_actual))
        a_bound.append(outward_cent(a_actual))
        b_bound.append(outward_cent(b_actual))
        c_bound.append(outward_cent(c_actual))
        d_bound.append(outward_cent(d_actual))
    return (tuple(den), tuple(num), tuple(a_bound), tuple(b_bound),
            tuple(c_bound), tuple(d_bound))


def main() -> None:
    print("EXACT RESIDUAL ORBIT CERTIFICATE v43")
    print("pair: imprimitive companion indices {1,22}")

    c1 = core.add(core.add(core.power_trace(1), core.power_trace(2)),
                  core.power_trace(5))
    c22 = core.add(core.add(core.power_trace(22), core.power_trace(23)),
                   core.power_trace(26))
    assert c1 == (Q(-2), Q(6), Q(1), Q(-5), Q(0), Q(1))
    assert c22 == core.neg(c1)
    y2 = core.sub(core.pi, core.scale(Q(1, 4), core.mul(c1, c1)))
    assert y2 == (Q(4), Q(21, 4), Q(-11, 4), Q(-13, 4), Q(1, 2), Q(1, 2))
    assert core.sign_vector(y2) == (1,) * core.DEGREE
    assert core.norm(y2) == Q(1, 4096)
    expected_height_boxes = tuple(
        (Q(lo), Q(hi)) for lo, hi in (
            ("1.12397017", "1.12397018"),
            ("1.11405647", "1.11405648"),
            ("1.77903527", "1.77903528"),
            ("0.04405416", "0.04405417"),
            ("0.21851142", "0.21851143"),
            ("0.72863891", "0.72863892"),
        )
    )
    for i, claimed in enumerate(expected_height_boxes):
        actual = core.sqrt_interval(
            core.interval_eval(y2, core.ROOT_INTERVALS[i]))
        assert claimed[0] < actual[0] < actual[1] < claimed[1]
    print("certified height boxes:", expected_height_boxes)

    bounds = derive_bounds(y2)
    expected_bounds = (
        tuple(map(Q, ("5.30", "2.83", "4.61", "0.51", "1.03", "3.31"))),
        tuple(map(Q, ("13.53", "5.90", "0.82", "0.62", "1.70", "6.32"))),
        tuple(map(Q, ("12.04", "5.30", "0.46", "14.03", "7.76", "8.68"))),
        tuple(map(Q, ("48.46", "11.80", "1.87", "29.44", "6.13", "12.89"))),
        tuple(map(Q, ("4.72", "2.54", "2.59", "11.53", "4.70", "4.54"))),
        tuple(map(Q, ("18.99", "5.66", "10.51", "24.19", "3.71", "6.75"))),
    )
    assert bounds == expected_bounds
    labels = ("denominator", "numerator", "a", "b", "c", "d")
    for label, bound in zip(labels, bounds):
        print(f"{label} outward bounds: {bound}")

    den, num, a_bound, b_bound, c_bound, d_bound = bounds
    core.Y2 = y2
    core.DEN_BOUND = den
    core.NUM_BOUND = num
    core.A_BOUND = a_bound
    core.B_BOUND = b_bound
    core.C_BOUND = c_bound
    core.D_BOUND = d_bound

    a_values = core.enumerate_trace_ellipsoid(
        core.ZERO, sum((x * x for x in a_bound), Q(0)),
        lambda x: core.within_integral_box(x, a_bound))
    c_values = core.enumerate_trace_ellipsoid(
        core.ZERO, sum((x * x for x in c_bound), Q(0)),
        lambda x: core.within_integral_box(x, c_bound))
    d_values = core.enumerate_trace_ellipsoid(
        core.ZERO, sum((x * x for x in d_bound), Q(0)),
        lambda x: core.within_integral_box(x, d_bound))
    print("exact separate entry counts a,c,d:",
          len(a_values), len(c_values), len(d_values))
    assert (len(a_values), len(c_values), len(d_values)) == (2641, 715, 65009)

    result1 = core.enumerate_case(
        1, core.scale(Q(-1, 2), c1), a_values, c_values)
    result22 = core.enumerate_case(
        22, core.scale(Q(-1, 2), c22), a_values, c_values)
    expected_histogram = {
        (0, 0): 1932,
        (0, 1): 5502,
        (0, 2): 966,
        (0, 3): 546,
        (1, 0): 84,
        (1, 1): 42,
    }
    assert result1[:3] == result22[:3] == (1849, 9072, 0)
    assert result1[3] == result22[3] == expected_histogram
    assert sum(expected_histogram.values()) == 9072
    assert result1[4] == "fac4b968a1d907fd25d4877ae7373a930a1ecb23ee6685926a26e9d686cfe22d"
    assert result22[4] == "6a08e191a62fdc81807e1c6447f83c7a0afa78a8ee170aa4edcf92752b73951a"
    print("first-failure histogram k=1:", sorted(result1[3].items()))
    print("first-failure histogram k=22:", sorted(result22[3].items()))
    print("PASS: exact exhaustive exclusion for both imprimitive residual indices")


if __name__ == "__main__":
    main()
