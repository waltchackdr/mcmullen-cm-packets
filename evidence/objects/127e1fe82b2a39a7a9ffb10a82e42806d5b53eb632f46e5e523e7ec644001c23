#!/usr/bin/env python3
"""Bounded exact verification of the September 2026 homology corrections.

Only Python's standard library is used.  No other repository module, stored
output, symbolic algebra system, or floating-point computation is used.
See README_principal_homology_v56.md for the mathematical scope.
"""

from fractions import Fraction
from itertools import combinations, product


def scalar_rank(rows, prime=0):
    """Gaussian elimination over Q or F_prime; the empty matrix has rank 0."""
    a = [[Fraction(x) if not prime else x % prime for x in r] for r in rows]
    if not a:
        return 0
    pivot = 0
    for j in range(len(a[0])):
        k = next((i for i in range(pivot, len(a)) if a[i][j]), None)
        if k is None:
            continue
        a[pivot], a[k] = a[k], a[pivot]
        inv = 1 / a[pivot][j] if not prime else pow(a[pivot][j], -1, prime)
        a[pivot] = [x * inv if not prime else x * inv % prime for x in a[pivot]]
        for i in range(pivot + 1, len(a)):
            q = a[i][j]
            a[i] = [x - q * y if not prime else (x - q * y) % prime
                    for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def determinant(rows):
    """Fraction elimination, independently of any number-field norm formula."""
    a = [list(map(Fraction, row)) for row in rows]
    n = len(a)
    ans = Fraction(1)
    for j in range(n):
        k = next((i for i in range(j, n) if a[i][j]), None)
        if k is None:
            return Fraction(0)
        if k != j:
            a[j], a[k] = a[k], a[j]
            ans = -ans
        q = a[j][j]
        ans *= q
        for i in range(j + 1, n):
            r = a[i][j] / q
            for l in range(j + 1, n):
                a[i][l] -= r * a[j][l]
    return ans


def coordinates(columns, target):
    """Solve a square rational coordinate system with exact pivots."""
    n = len(target)
    assert len(columns) == n
    a = [[Fraction(columns[j][i]) for j in range(n)] + [Fraction(target[i])]
         for i in range(n)]
    for j in range(n):
        k = next(i for i in range(j, n) if a[i][j])
        a[j], a[k] = a[k], a[j]
        q = a[j][j]
        a[j] = [x / q for x in a[j]]
        for i in range(n):
            if i != j:
                q = a[i][j]
                a[i] = [x - q * y for x, y in zip(a[i], a[j])]
    return [row[-1] for row in a]


class Quotient:
    """Dense polynomial quotient with a monic modulus, over Q or F_p."""

    def __init__(self, modulus, prime=0):
        assert modulus[-1] == 1
        self.modulus = tuple(modulus)
        self.degree = len(modulus) - 1
        self.prime = prime

    def __call__(self, x):
        if isinstance(x, Element):
            assert x.ring is self
            return x
        return Element(self, x if isinstance(x, (tuple, list)) else [x])

    def scalar(self, x):
        q = Fraction(x)
        return q if not self.prime else q.numerator * pow(q.denominator, -1, self.prime) % self.prime


class Element:
    def __init__(self, ring, coefficients):
        self.ring = ring
        c = [ring.scalar(x) for x in coefficients]
        n = ring.degree
        c += [ring.scalar(0)] * max(0, n - len(c))
        for k in range(len(c) - 1, n - 1, -1):
            q = c[k]
            for j in range(n):
                c[k - n + j] = ring.scalar(c[k - n + j] - q * ring.modulus[j])
        self.c = tuple(c[:n])

    def __add__(self, other):
        other = self.ring(other)
        return self.ring([x + y for x, y in zip(self.c, other.c)])

    __radd__ = __add__

    def __neg__(self):
        return self.ring([-x for x in self.c])

    def __sub__(self, other):
        return self + -self.ring(other)

    def __rsub__(self, other):
        return self.ring(other) - self

    def __mul__(self, other):
        other = self.ring(other)
        c = [0] * (2 * self.ring.degree - 1)
        for i, x in enumerate(self.c):
            for j, y in enumerate(other.c):
                c[i + j] += x * y
        return self.ring(c)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        answer = self.ring(1)
        base = self
        while n:
            if n & 1:
                answer *= base
            base *= base
            n >>= 1
        return answer

    def __eq__(self, other):
        return self.c == self.ring(other).c

    def __bool__(self):
        return any(self.c)

    def multiplication_columns(self):
        x = self.ring([0, 1])
        return [(self * x**j).c for j in range(self.ring.degree)]

    def inverse(self):
        assert self.ring.prime == 0
        return self.ring(coordinates(self.multiplication_columns(), self.ring(1).c))

    def trace(self):
        return sum(c[j] for j, c in enumerate(self.multiplication_columns()))


def eye(ring):
    return [[ring(1), ring(0)], [ring(0), ring(1)]]


def transpose(m):
    return [list(row) for row in zip(*m)]


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def ma(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(a, q):
    return [[q * x for x in row] for row in a]


def mpow(a, n):
    out = eye(a[0][0].ring)
    # Linear iteration is intentional: it does not reuse the primary matrix
    # binary-power implementation and exponents here are at most 42.
    for _ in range(n):
        out = mm(out, a)
    return out


def act(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def det2(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def flatten(v):
    return v[0].c + v[1].c


def row_rank(vectors, p):
    return scalar_rank([flatten(v) for v in vectors], p)


def wedge(v, w):
    # Literal Eq. (E1): Tr(a'b - ab'), not its negative.
    return w[0] * v[1] - v[0] * w[1]


def integral_matrix(m):
    return all(c.denominator == 1 for row in m for x in row for c in x.c)


def rational_operator(m, basis):
    return transpose([flatten(act(m, v)) for v in basis])


def verify():
    f = [1, 8, 8, -6, -6, 1, 1]
    L = Quotient(f)
    t = L([0, 1])
    one, zero = L(1), L(0)
    I = eye(L)
    O = scale(I, 0)
    J = [[zero, -one], [one, zero]]
    powers = [t**j for j in range(6)]
    rational_basis = [(x, zero) for x in powers] + [(zero, x) for x in powers]
    E = [[wedge(v, w).trace() for w in rational_basis] for v in rational_basis]
    trace_gram = [[(x * y).trace() for y in powers] for x in powers]
    discriminant = determinant(trace_gram)
    assert discriminant == 453789 == 3**3 * 7**5
    assert determinant(E) == discriminant**2

    # Printed McMullen data, then conjugation by diag(t,1).
    printed = {
        "A": [[t**5 - 4*t**3 - t**2 + 3*t,
               -t**5 + t**4 + 6*t**3 - 3*t**2 - 10*t],
              [t**3 - 1, -t**5 + 5*t**3 + t**2 - 6*t]],
        "B": [[zero, -one], [one, t]],
        "C": [[2*t - 4*t**2 - 4*t**3 + t**4 + t**5,
               -5*t - 4*t**2 + t**3 + t**4],
              [-1 + 5*t + 4*t**2 - t**3 - t**4,
               -t + 4*t**2 + 4*t**3 - t**4 - t**5]],
    }
    h = [[t, zero], [zero, one]]
    hi = [[t.inverse(), zero], [zero, one]]
    assert integral_matrix(hi)
    vertices = {name: mm(mm(h, a), hi) for name, a in printed.items()}
    assert all(det2(a) == 1 and integral_matrix(a) for a in vertices.values())

    # Different construction of A's quadratic subfield operator: perform a
    # 12-by-12 rational basis change in the cyclotomic field itself.  Neither
    # a stored coefficient formula nor the quadratic recurrence is used.
    phi42 = [1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1]
    M = Quotient(phi42)
    z = M([0, 1])
    assert z**42 == 1 and all(z**k != 1 for k in (6, 14, 21))
    tz = z + z**41
    assert sum(a * tz**j for j, a in enumerate(f)) == 0
    cyclotomic_basis = [tz**j for j in range(6)] + [z**3 * tz**j for j in range(6)]
    assert scalar_rank([b.c for b in cyclotomic_basis]) == 12

    def substitute_a(x):
        c = coordinates([b.c for b in cyclotomic_basis], x.c)
        a, b = L(c[:6]), L(c[6:])
        return ma(scale(I, a), scale(vertices["A"], b))

    ZA = substitute_a(z)
    assert integral_matrix(ZA) and mpow(ZA, 3) == vertices["A"]
    assert ma(ma(mpow(ZA, 2), scale(ZA, -t)), I) == O
    operators = {"A": substitute_a(1 + 2*z**14)}
    operators.update({name: ma(I, scale(mpow(vertices[name], 14), 2)) for name in ("B", "C")})
    assert operators["A"] == ma(I, scale(mpow(ZA, 14), 2))
    homology = {name: transpose(s) for name, s in operators.items()}
    units = {name: scale(ma(a, scale(I, -1)), Fraction(1, 2)) for name, a in homology.items()}
    Jq = rational_operator(J, rational_basis)
    assert mm(mm(transpose(Jq), E), Jq) == E
    for name, s in operators.items():
        a = homology[name]
        assert s[0][0] + s[1][1] == 0
        assert mpow(s, 2) == scale(I, -3)
        assert integral_matrix(s) and integral_matrix(units[name])
        assert ma(ma(mpow(units[name], 2), units[name]), I) == O
        assert a == mm(mm(J, s), J)  # S^T = -J S J^{-1}, J^{-1}=-J.
        aq = rational_operator(a, rational_basis)
        assert ma(mm(transpose(aq), E), mm(E, aq)) == [[0]*12 for _ in range(12)]
        uq = rational_operator(units[name], rational_basis)
        assert mm(mm(transpose(uq), E), uq) == E
        print(f"{name}: integral S and U; trace 0; S^2=-3; actual action S^T; Rosati=-S^T")

    # Direct Gram calculations on the actual homology basis f_j=(0,t^j)
    # and the old matrix-action basis e_j=(t^j,0).  K-basis status is
    # checked as rational rank 12, not presumed.
    legacy = {"A": -2**6 * 3**6 * 7**5 * 127,
              "B": -2**6 * 3**6 * 7**5,
              "C": -2**6 * 3**6 * 7**5 * 127}
    k_basis = rational_basis[:6]
    for name, s in operators.items():
        a = homology[name]
        sq = rational_operator(s, rational_basis)
        aq = rational_operator(a, rational_basis)
        # Real and imaginary coefficient matrices of H_{S^T}(Jv,Jw)
        # equal those of -conjugate(H_S(v,w)) on all of V, respectively.
        assert mm(mm(transpose(Jq), mm(transpose(aq), E)), Jq) == scale(mm(transpose(sq), E), -1)
        assert mm(mm(transpose(Jq), E), Jq) == E
        old = [[wedge(act(s, v), w).trace() for w in k_basis] for v in k_basis]
        transported_basis = [act(J, v) for v in k_basis]
        transported = [[wedge(act(a, v), w).trace() for w in transported_basis]
                       for v in transported_basis]
        for operator, base in ((s, k_basis), (a, transported_basis)):
            assert scalar_rank([flatten(v) for v in base] +
                               [flatten(act(operator, v)) for v in base]) == 12
            assert all(wedge(v, w) == 0 for v in base for w in base)
        assert transported == transpose(transported)
        assert transported == scale(old, -1)
        assert determinant(old) == determinant(transported) == legacy[name]
        value = determinant(transported)
        ratio = value / legacy[name]
        norm_witness = Fraction(1)
        assert ratio == norm_witness**2  # N_{Q(sqrt(-3))/Q}(norm_witness).
        print(f"{name}: actual f_j=(0,t^j) Gram determinant {value}; legacy {legacy[name]}; ratio N({norm_witness})={ratio}")
    assert 8**2 == 2**6 and 27**2 == 3**6
    assert 2**2 + 3 == 7 and 10**2 + 3*3**2 == 127
    print("J transport checked on all rational basis pairs; all three determinant classes remain -1")

    # Establish the quotient rings used for the local computations from the
    # different generator f'(t), including its trace-dual interpretation.
    derivative = sum(j * f[j] * t**(j-1) for j in range(1, 7))
    derivative_inverse = derivative.inverse()
    dual_columns = transpose([(x * derivative_inverse).c for x in powers])
    dual_pairings = mm(trace_gram, dual_columns)
    assert all(q.denominator == 1 for row in dual_pairings for q in row)
    assert abs(determinant(dual_pairings)) == 1
    assert abs(determinant(transpose(derivative.multiplication_columns()))) == discriminant
    F27 = Quotient([1, 1, -1, 1], 3)
    R = Quotient([0, 0, 0, 0, 0, 1], 7)
    v3, u = F27([0, 1]), R([0, 1])
    assert all((F27([a])**3 - F27([a])**2 + F27([a]) + 1) != 0 for a in range(3))
    assert sum(f[j]*v3**j for j in range(7)) == 0
    assert sum(f[j]*(1+u)**j for j in range(7)) == 0
    assert sum(j*f[j]*v3**(j-1) for j in range(1, 7)) == 0
    assert sum(j*f[j]*(1+u)**(j-1) for j in range(1, 7)) == 0
    assert 3**F27.degree * 7**R.degree == discriminant
    # The maps onto these coprime quotients are surjective: t generates F27
    # and t-1 generates R.  Their product has index disc(L); the inclusion
    # (f') in the joint kernel above is therefore equality.

    def local(x, ring, root):
        return sum(c * root**j for j, c in enumerate(x.c))

    def local_matrix(m, ring, root):
        return [[local(x, ring, root) for x in row] for row in m]

    a7 = {name: local_matrix(a, R, 1+u) for name, a in homology.items()}
    a3 = {name: local_matrix(a, F27, v3) for name, a in homology.items()}
    u7 = {name: local_matrix(a, R, 1+u) for name, a in units.items()}
    u3 = {name: local_matrix(a, F27, v3) for name, a in units.items()}
    i7, i3 = eye(R), eye(F27)
    residue_expected = {"A": [[1, 3], [1, 6]], "B": [[6, 2], [5, 1]], "C": [[5, 5], [0, 2]]}
    h3 = [(v3**j, F27(0)) for j in range(3)]
    assert row_rank(h3, 3) == 3
    for name in operators:
        assert [[x.c[0] for x in row] for row in a7[name]] == residue_expected[name]
        assert a3[name] == scale(i3, 0) and u3[name] == i3
        for v, w in product(h3, repeat=2):
            assert wedge(v, w) == 0
        for v in h3:
            assert row_rank(h3 + [act(u3[name], v)], 3) == 3
        full_vertex = local_matrix(transpose(vertices[name]), F27, v3)
        vertex_trace = full_vertex[0][0] + full_vertex[1][1]
        assert det2(full_vertex) == 1
        vertex_inverse = ma(scale(i3, vertex_trace), scale(full_vertex, -1))
        assert mm(full_vertex, vertex_inverse) == i3
        characteristic_discriminant = vertex_trace**2 - 4
        assert characteristic_discriminant**13 == -1
        # A full-vertex-invariant finite additive subgroup is invariant under
        # the inverse and hence this scalar trace.  It therefore recovers t;
        # no O_L-stability assumption is needed for this local obstruction.
        assert (vertex_trace**9 if name == "A" else vertex_trace) == v3
    print("At p3 all S^T vanish; all full vertex characteristic discriminants satisfy d^13=-1 in F27")
    print("Full vertex plus inverse recovers t at p3: trace^9 for A, trace for B/C")

    # Test both corrected eigenprojectors by their complete F7 image spaces.
    full7 = [(u**j, R(0)) for j in range(5)] + [(R(0), u**j) for j in range(5)]
    for name, a in a7.items():
        eplus = scale(ma(a, scale(i7, 2)), pow(4, -1, 7))
        eminus = scale(ma(scale(i7, 2), scale(a, -1)), pow(4, -1, 7))
        assert ma(eplus, eminus) == i7
        assert mm(eplus, eminus) == scale(i7, 0)
        for sign, projector in ((1, eplus), (-1, eminus)):
            assert mm(projector, projector) == projector
            image = [act(projector, v) for v in full7]
            assert row_rank(image, 7) == 5
            generator = next(v for v in image if v[0].c[0] or v[1].c[0])
            shifts = [(u**j * generator[0], u**j * generator[1]) for j in range(5)]
            assert row_rank(shifts, 7) == row_rank(shifts + image, 7) == 5
            for v in image:
                assert act(a, v) == tuple(2*sign*x for x in v)
                assert row_rank(image + [act(u7[name], v)], 7) == 5
            assert all(wedge(v, w) == 0 for v, w in product(image, repeat=2))
        print(f"{name}: both S^T eigenprojectors have free isotropic R-image of order 7^5; H3 order 27")

    # The omitted nonfree common module.  Its five displayed basis vectors
    # span u^2 R(3,1)+u^3 R^2, using the complement (1,0).
    common = [(3*u**j, u**j) for j in (2, 3, 4)] + [(u**j, R(0)) for j in (3, 4)]
    assert row_rank(common, 7) == 5
    assert all(wedge(v, w) == 0 for v, w in product(common, repeat=2))
    assert all(tuple(u**3*x for x in v) == (R(0), R(0)) for v in common)
    assert any(common[0])  # annihilated by u^3, hence not a free rank-one R-module.
    for v in common:
        assert row_rank(common + [tuple(u*x for x in v)], 7) == 5
        for name in operators:
            assert row_rank(common + [act(a7[name], v)], 7) == 5
            assert row_rank(common + [act(u7[name], v)], 7) == 5
    assert (7**5)**2 == 7**10 and (3**3 * 7**5)**2 == discriminant**2
    print("Common H7=u^2 R(3,1)+u^3 R^2: dimension 5, order 16807, nonfree, isotropic, stable under actual A/B/C and integral U")
    print("Common H3=F27(1,0): order 27; product order 453789 is maximal isotropic")

    # Preserve the original failure witness and the narrower free-line fact.
    wrong = (R(1), R(2))
    wrong_s = act(transpose(a7["B"]), wrong)
    wrong_st = act(a7["B"], wrong)
    assert tuple(x.c[0] for x in wrong_s) == (2, 4)
    assert tuple(x.c[0] for x in wrong_st) == (3, 0)
    assert (wrong[0]*wrong_st[1] - wrong[1]*wrong_st[0]).c[0] == 1
    # All projective residue lines: (m,1) plus (1,0).  These are finite
    # exhaustive searches of 8 lines and then 7 first-order lifts each.
    all_lines = [(R(m), R(1)) for m in range(7)] + [(R(1), R(0))]
    common_lines = [v for v in all_lines if all(wedge(v, act(a, v)).c[0] == 0 for a in a7.values())]
    assert common_lines == [(R(3), R(1))]
    for x, y in combinations(a7.values(), 2):
        assert [v for v in all_lines if wedge(v, act(x, v)).c[0] == wedge(v, act(y, v)).c[0] == 0] == common_lines
    first_lifts = {}
    for name, a in a7.items():
        lifts = [m for m in range(7) if all(c == 0 for c in wedge((3+m*u, R(1)), act(a, (3+m*u, R(1)))).c[:2])]
        assert len(lifts) == 1
        first_lifts[name] = lifts[0]
    assert len(set(first_lifts.values())) == 3
    print("Actual homology free-line first-order lifts (3+a*u,1): " + ", ".join(f"{name}:a={a}" for name, a in first_lifts.items()))
    print("Retained counterexample: original B line (1,2) maps to (3,0) under S_B^T modulo 7")
    print("PASS: bounded correction verification; C18 unrestricted obstruction is false")
    print("NO global family, curve monodromy, exact surviving CM order, or new packet geometry asserted")


if __name__ == "__main__":
    verify()
