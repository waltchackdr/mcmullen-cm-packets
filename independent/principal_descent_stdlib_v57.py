#!/usr/bin/env python3
"""Exact finite support for the common-level descent argument.

Reuses only this implementation's own v56 exact arithmetic.  No primary
certificate is imported.  See README_principal_descent_v57.md for scope,
the dependency boundary, and the human arguments supported by these checks.
"""

from hashlib import sha256
from itertools import product
import json

from principal_homology_stdlib_v56 import (
    Quotient, act, coordinates, det2, determinant, eye, integral_matrix,
    ma, mm, mpow, row_rank, scale, transpose,
)


def local(x, ring, root):
    return sum(c * root**j for j, c in enumerate(x.c))


def local_matrix(m, ring, root):
    return [[local(x, ring, root) for x in row] for row in m]


def in_common_h7(v):
    """Exact membership in u^2 R(3,1)+u^3 R^2, in power coordinates."""
    x, y = v
    return all(c == 0 for c in x.c[:2] + y.c[:2]) and (x.c[2] - 3*y.c[2]) % 7 == 0


def preserves_h7(m, basis):
    return all(in_common_h7(act(m, v)) for v in basis)


def verify_stabilizers(R, u, h7, F27):
    # H7 is a five-dimensional subspace.  The coordinate membership test
    # imposes five independent conditions on F7^10 and contains its basis.
    assert row_rank(h7, 7) == 5 and all(in_common_h7(v) for v in h7)
    # Every matrix with zero residue is an F7-linear combination of these
    # 16 monomial matrices.  Thus all lifts have the same H7-stability test.
    for i, j, k in product(range(2), range(2), range(1, 5)):
        perturbation = scale(eye(R), 0)
        perturbation[i][j] = u**k
        assert preserves_h7(perturbation, h7)

    def line7(x, y):
        assert x % 7 or y % 7
        return (x * pow(y, -1, 7) % 7, 1) if y % 7 else (1, 0)

    matrices7, orbit7, stabilizer7 = [], set(), 0
    all_lines7 = {(m, 1) for m in range(7)} | {(1, 0)}
    for a, b, c, d in product(range(7), repeat=4):
        if (a*d - b*c) % 7 != 1:
            continue
        matrices7.append((a, b, c, d))
        image = line7((3*a+b) % 7, (3*c+d) % 7)
        orbit7.add(image)
        m = [[R(a), R(b)], [R(c), R(d)]]
        stable = preserves_h7(m, h7)
        assert stable == (image == (3, 1))
        stabilizer7 += stable
    assert len(matrices7) == 336 and stabilizer7 == 42
    assert orbit7 == all_lines7 and len(orbit7) == 8
    print("p7: all 336 SL2(F7) matrices; H7 stabilizer 42; line orbit 8")
    print("p7: 16 zero-residue matrix directions preserve H7; its stabilizer is exactly the residue-line preimage")

    # Derive finite-field lookup tables from the quotient implementation;
    # this makes exhaustive enumeration of SL2(F27) small and transparent.
    elements = [F27([i % 3, (i // 3) % 3, i // 9]) for i in range(27)]
    labels = {x.c: i for i, x in enumerate(elements)}
    add = [[labels[(x+y).c] for y in elements] for x in elements]
    mul = [[labels[(x*y).c] for y in elements] for x in elements]
    negative = [labels[(-x).c] for x in elements]
    inverse = {a: next(b for b in range(1, 27) if mul[a][b] == 1) for a in range(1, 27)}

    def line27(x, y):
        assert x or y
        return (mul[x][inverse[y]], 1) if y else (1, 0)

    def determinant27(a, b, c, d):
        return add[mul[a][d]][negative[mul[b][c]]]

    group_count, stabilizer27, orbit27 = 0, 0, set()
    # Every nonzero first column has exactly 27 possible second columns.
    # Parameterizing the determinant-one equation avoids a 27^4 search.
    for a, c in product(range(27), repeat=2):
        if not a and not c:
            continue
        if a:
            second_columns = [(b, mul[add[1][mul[b][c]]][inverse[a]]) for b in range(27)]
        else:
            second_columns = [(negative[inverse[c]], d) for d in range(27)]
        assert len(set(second_columns)) == 27
        for b, d in second_columns:
            assert determinant27(a, b, c, d) == 1
            group_count += 1
            image = line27(a, c)
            orbit27.add(image)
            stabilizer27 += image == (1, 0)
    all_lines27 = {(m, 1) for m in range(27)} | {(1, 0)}
    assert group_count == 19656 and stabilizer27 == 702
    assert orbit27 == all_lines27 and len(orbit27) == 28

    # Verify transitivity already with elementary matrices whose entries
    # are the F3-basis 1,t,t^2; CRT lifts are treated in the human proof.
    todo, reached = [(1, 0)], {(1, 0)}
    while todo:
        x, y = todo.pop()
        for coefficient in (1, 3, 9):
            for pair in ((add[x][mul[coefficient][y]], y),
                         (x, add[y][mul[coefficient][x]])):
                image = line27(*pair)
                if image not in reached:
                    reached.add(image)
                    todo.append(image)
    assert reached == orbit27
    print("p3: all 19656 SL2(F27) matrices; H3 stabilizer 702; line orbit 28")
    joint_orbit = set(product(orbit7, orbit27))
    assert len(joint_orbit) == 224
    print("Product residue orbit: 8*28=224; ambient SL2(O_L) index uses CRT elementary lifts")
    return {"elements": elements, "labels": labels, "add": add,
            "mul": mul, "inverse": inverse}


def verify_curve_orbit(vertices, F27, v3, residue_tables):
    """Orbit for the explicitly normalized vertex group, with rho(g)=g^-T."""
    add, mul, inverse = (residue_tables[k] for k in ("add", "mul", "inverse"))
    labels = residue_tables["labels"]
    generators = {}
    for name, g in vertices.items():
        # The family transport is inverse transpose, distinct from the
        # transpose describing the fixed-fiber CM endomorphism.
        a, b = g[0]
        c, d = g[1]
        rho = [[d, -c], [-b, a]]
        assert mm(transpose(g), rho) == eye(a.ring)
        m3 = local_matrix(rho, F27, v3)
        matrix3 = [[labels[x.c] for x in row] for row in m3]
        matrix7 = [[sum(int(c) for c in x.c) % 7 for x in row] for row in rho]
        generators[name] = matrix3, matrix7

    def image(state, generator):
        (x, y), (v, w) = state
        m3, m7 = generator
        x1 = add[mul[m3[0][0]][x]][mul[m3[0][1]][y]]
        y1 = add[mul[m3[1][0]][x]][mul[m3[1][1]][y]]
        assert x1 or y1
        line3 = (mul[x1][inverse[y1]], 1) if y1 else (1, 0)
        v1 = (m7[0][0]*v + m7[0][1]*w) % 7
        w1 = (m7[1][0]*v + m7[1][1]*w) % 7
        assert v1 or w1
        line7 = (v1*pow(w1, -1, 7) % 7, 1) if w1 else (1, 0)
        return line3, line7

    lines3 = [(m, 1) for m in range(27)] + [(1, 0)]
    lines7 = [(m, 1) for m in range(7)] + [(1, 0)]
    domain = set(product(lines3, lines7))
    assert len(domain) == 224
    start = ((1, 0), (3, 1))
    for generator in generators.values():
        assert {image(state, generator) for state in domain} == domain
        assert image(start, generator)[1] == (3, 1)
    reached, todo = {start}, [start]
    while todo:
        state = todo.pop()
        for generator in generators.values():
            target = image(state, generator)
            if target not in reached:
                reached.add(target)
                todo.append(target)
    projection3 = {state[0] for state in reached}
    projection7 = {state[1] for state in reached}
    print(f"Normalized vertex group rho(g)=g^(-T): orbit {len(reached)}; projections p3={len(projection3)}, p7={len(projection7)}")
    assert reached == set(product(lines3, [(3, 1)]))
    assert len(reached) == 28 and len(projection3) == 28 and len(projection7) == 1
    minus_identity = scale(eye(vertices["A"][0][0].ring), -1)
    assert mpow(vertices["A"], 7) == minus_identity
    for name in ("B", "C"):
        assert mpow(vertices[name], 21) == minus_identity
    # Canonical field labels are c0+3*c1+9*c2.  State ordering and transition
    # arrays are fully specified so another implementation can reproduce it.
    states = sorted(reached)
    indices = {state: i for i, state in enumerate(states)}
    transitions = {name: [indices[image(state, generator)] for state in states]
                   for name, generator in generators.items()}
    record = {"states": states, "permutations": transitions}
    digest = sha256(json.dumps(record, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    print("Normalized curve orbit/transition SHA256: " + digest)
    assert digest == "4ccb0a95451b2c4e77a8b58996051a8a5dbcfcb5a625a2d6863a1348a0c7575e"
    print("Explicit Delta0 stabilizer index 28 (given A/B/C generation); -I is present and projectively trivial")


def reconstruct(L, t):
    one, zero = L(1), L(0)
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
    h, hi = [[t, zero], [zero, one]], [[t.inverse(), zero], [zero, one]]
    vertices = {name: mm(mm(h, x), hi) for name, x in printed.items()}
    assert all(integral_matrix(x) and det2(x) == 1 for x in vertices.values())
    M = Quotient([1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1])
    z = M([0, 1])
    tz = z + z**41
    assert z**42 == 1 and all(z**k != 1 for k in (6, 14, 21))
    assert sum(c*tz**j for j, c in enumerate(L.modulus)) == 0

    def relative_coordinates(x, root):
        basis = [tz**j for j in range(6)] + [root*tz**j for j in range(6)]
        c = coordinates([v.c for v in basis], x.c)
        return L(c[:6]), L(c[6:])

    cyclotomic_integral_basis = [tz**j for j in range(6)] + [z*tz**j for j in range(6)]
    # This identifies Z[z] with O_L+O_L z integrally, not only rationally.
    assert all(c.denominator == 1 for x in cyclotomic_integral_basis for c in x.c)
    assert abs(determinant(transpose([x.c for x in cyclotomic_integral_basis]))) == 1
    root_a, root_b = relative_coordinates(z, z**3)
    ZA = ma(scale(eye(L), root_a), scale(vertices["A"], root_b))
    assert integral_matrix(ZA) and mpow(ZA, 3) == vertices["A"]
    primitive = {"A": ZA, "B": vertices["B"], "C": vertices["C"]}
    for x in primitive.values():
        assert integral_matrix(x) and det2(x) == 1 and x[0][0] + x[1][1] == t
        assert ma(ma(mpow(x, 2), scale(x, -t)), eye(L)) == scale(eye(L), 0)
    return vertices, primitive, relative_coordinates(z**14, z)


def verify_order(L, t, F27, v3, relative):
    a, b = relative
    assert all(c.denominator == 1 for x in (a, b) for c in x.c)
    assert local(b, F27, v3) == 0
    norm_b = determinant(transpose(b.multiplication_columns()))
    assert abs(norm_b) == 27
    print("z^14=a+b*z: a=" + str(tuple(int(c) for c in a.c)) +
          "; b=" + str(tuple(int(c) for c in b.c)))
    print(f"N_L(b)={norm_b}; b mod p3=0; (b)=p3 by the common ideal index 27")
    # A second lattice-index computation in the absolute twelve-dimensional
    # basis checks the order generated by z^14 directly.
    powers = [t**j for j in range(6)]
    columns = [x.c + (L(0)).c for x in powers]
    columns += [(a*x).c + (b*x).c for x in powers]
    assert abs(determinant(transpose(columns))) == 27
    assert abs(determinant(transpose((t**2-4).multiplication_columns()))) == 1
    print("R_K=O_L[O_K]=O_L+p3*O_M: index 27 in O_M; relative discriminant norm 1")


def cyclic_survivors(matrix, limit, F27, v3, R, u, h7):
    identity = eye(matrix[0][0].ring)
    current = identity
    local3, local7 = local_matrix(matrix, F27, v3), local_matrix(matrix, R, 1+u)
    current3, current7 = eye(F27), eye(R)
    surviving3, surviving7, surviving_both = [], [], []
    order = None
    for n in range(limit + 1):
        if n and current == identity:
            order = n
            break
        stable3 = current3[1][0] == 0  # image of F27(1,0)
        stable7 = preserves_h7(current7, h7)
        if stable3:
            surviving3.append(n)
        if stable7:
            surviving7.append(n)
        if stable3 and stable7:
            surviving_both.append(n)
        current = mm(current, matrix)
        current3 = mm(current3, local3)
        current7 = mm(current7, local7)
    assert order is not None, "No exact cyclic order found within the declared finite bound"
    return order, surviving3, surviving7, surviving_both


def verify():
    L = Quotient([1, 8, 8, -6, -6, 1, 1])
    t = L([0, 1])
    F27, R = Quotient([1, 1, -1, 1], 3), Quotient([0, 0, 0, 0, 0, 1], 7)
    v3, u = F27([0, 1]), R([0, 1])
    h7 = [(3*u**j, u**j) for j in (2, 3, 4)] + [(u**j, R(0)) for j in (3, 4)]
    residue_tables = verify_stabilizers(R, u, h7, F27)
    vertices, primitive, relative = reconstruct(L, t)
    verify_curve_orbit(vertices, F27, v3, residue_tables)
    verify_order(L, t, F27, v3, relative)

    for name, primitive_matrix in primitive.items():
        z = transpose(primitive_matrix)
        z3, z7 = local_matrix(z, F27, v3), local_matrix(z, R, 1+u)
        assert preserves_h7(z7, h7)
        assert z3[1][0] != 0
        # Within O_M=O_L+O_L z, enumerate its 729 residue coefficient pairs.
        # Exactly the scalar coefficients preserve H3, which is equivalent
        # to the z coefficient lying in p3.  The p7 condition adds nothing,
        # since z itself preserves H7.
        residue = [F27(c) for c in product(range(3), repeat=3)]
        stabilizer_coefficients = []
        for alpha, beta in product(residue, repeat=2):
            m = ma(scale(eye(F27), alpha), scale(z3, beta))
            if m[1][0] == 0:
                stabilizer_coefficients.append((alpha.c, beta.c))
                assert beta == 0
        assert len(stabilizer_coefficients) == 27
        print(f"{name}: integral CM multipliers modulo p3: 27/729 coefficient pairs preserve H3; all O_M preserves H7")

        primitive_result = cyclic_survivors(z, 42, F27, v3, R, u, h7)
        print(f"{name}: primitive CM generator order {primitive_result[0]}; surviving exponents {primitive_result[3]}")
        # These comparisons are regression checks after the finite search.
        assert primitive_result[0] == 42
        assert primitive_result[2] == list(range(42))
        assert primitive_result[3] == list(range(0, 42, 7))
        vertex_result = cyclic_survivors(transpose(vertices[name]), 42, F27, v3, R, u, h7)
        print(f"{name}: full vertex generator order {vertex_result[0]}; surviving exponents {vertex_result[3]}")
        expected_order = 14 if name == "A" else 42
        assert vertex_result[0] == expected_order
        assert vertex_result[2] == list(range(expected_order))
        assert vertex_result[3] == list(range(0, expected_order, 7))

    print("PASS: ambient index 224; explicit normalized Delta0 index 28; multiplier order index 27; vertex survivors A=mu2, B/C=mu6")
    print("Curve index is before any auxiliary level; no index for an unnamed normalization or new packet geometry asserted")


if __name__ == "__main__":
    verify()
