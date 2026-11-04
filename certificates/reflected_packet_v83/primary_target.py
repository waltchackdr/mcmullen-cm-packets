"""Primary exact target reconstruction for the two holomorphic chambers.

SymPy polynomial arithmetic reconstructs normalized vertex matrices and their
conjugates. The anti-Mobius formula is checked separately. The reflection is
ONLY a target transformation; it is not used as an SL2 source equivalence.
"""
from fractions import Fraction as Q

TARGET_IDS = ("ORIGINAL_Q_V1", "REFLECTED_AB_V1")
EMBEDDINGS = (1, 5, 11, 13, 17, 19)
BOUND_NUMERATORS = {
    "ORIGINAL_Q_V1": ((4006, 13933, 8395, 17094, 20761, 6663),
                      (255285, 208594, 17707, 121724, 165248, 191115)),
    "REFLECTED_AB_V1": ((1368, 23817, 14522, 22806, 10088, 5105),
                        (197767, 354283, 28993, 254603, 165248, 191115)),
}


def polynomial(primary, coefficients):
    return sum(primary.sp.Rational(q.numerator, q.denominator) * primary.ambient.T**i
               for i, q in enumerate(map(Q, coefficients)))


def matrix_from_core(primary, matrix):
    return primary.sp.Matrix([[polynomial(primary, entry) for entry in row]
                              for row in matrix])


def vertex_invariants(primary, matrix):
    """Determinant-one elliptic vertex, not the norm-pi source formula."""
    a = primary.ambient
    assert a.reduce_l(matrix.det()) == 1
    trace = a.reduce_l(matrix.trace())
    lower = matrix[1, 0]
    assert lower != 0
    x = a.reduce_field((matrix[0, 0] - matrix[1, 1]) / (2 * lower))
    y = a.reduce_field((4 - trace**2) / (4 * lower**2))
    modulus = a.reduce_l(x*x + y)
    assert a.sign_vector(y) == a.sign_vector(modulus) == (1,)*6
    return x, y, modulus


def anti_image(primary, reflection, point):
    """Exact invariant check; the unsquared height is -det(R)*h/den > 0."""
    a = primary.ambient
    x, y, _ = point
    r, s, u, v = tuple(reflection)
    determinant = a.reduce_l(reflection.det())
    assert a.sign_vector(determinant) == (-1,)*6
    den = a.reduce_l((u*x+v)**2 + u*u*y)
    assert a.sign_vector(den) == (1,)*6
    image_x = a.reduce_field(((r*x+s)*(u*x+v)+r*u*y)/den)
    image_y = a.reduce_field(determinant**2*y/(den*den))
    image_modulus = a.reduce_field(((r*x+s)**2+r*r*y)/den)
    assert a.reduce_l(image_x*image_x+image_y-image_modulus) == 0
    return image_x, image_y, image_modulus


def derive_target(primary, target_id):
    if target_id not in TARGET_IDS:
        raise ValueError("Unknown exact target chamber")
    a, sp, core = primary.ambient, primary.sp, primary.core
    assert core.EMBEDDINGS == EMBEDDINGS
    # Rebuild from immutable normalized matrices, never a previously installed
    # target's mutable vertex invariants or sides.
    original_actions = tuple(matrix_from_core(primary, value) for value in core.VERTICES)
    original = tuple(vertex_invariants(primary, value) for value in original_actions)
    reflection = determinant = determinant_signs = None
    if target_id == "REFLECTED_AB_V1":
        A, B = original_actions[:2]
        reflection = a.matrix_reduce(A*B-B*A)
        determinant = a.reduce_l(reflection.det())
        determinant_signs = a.sign_vector(determinant)
        assert determinant_signs == (-1,)*6
        assert a.absolute_norm(determinant) == 21952
        assert a.reduce_l(reflection.trace()) == 0
        assert a.matrix_reduce(reflection*reflection + determinant*sp.eye(2)) == sp.zeros(2)
        inverse = a.matrix_reduce(-reflection/determinant)
        assert a.matrix_multiply(reflection, inverse) == sp.eye(2)
        # Conjugate matrices first, using the primary matrix arithmetic.
        actions = tuple(a.matrix_multiply(a.matrix_multiply(reflection, value), inverse)
                        for value in original_actions)
        points = tuple(vertex_invariants(primary, value) for value in actions)
        # Independent invariant formula checks the antiholomorphic height sign.
        assert points == tuple(anti_image(primary, reflection, point) for point in original)
        assert tuple(anti_image(primary, reflection, point) for point in points) == original
        assert points[:2] == original[:2]
    else:
        points = original
    sides = []
    for index in range(4):
        left, right = points[index], points[(index+1) % 4]
        difference = a.reduce_l(right[0]-left[0])
        assert difference != 0  # Both targets have only nonvertical sides.
        beta = a.reduce_field((right[2]-left[2])/(2*difference))
        radius = a.reduce_l((left[0]-beta)**2+left[1])
        assert a.sign_vector(radius) == (1,)*6
        values = [a.reduce_l(modulus-2*beta*x+beta*beta-radius) for x, _, modulus in points]
        assert values[index] == values[(index+1) % 4] == 0
        sign = a.sign_vector(values[(index+2) % 4])
        assert 0 not in sign and sign == a.sign_vector(values[(index+3) % 4])
        sides.append((primary.to_core(beta), primary.to_core(radius), sign))
    m, U = (tuple(Q(value, 100000) for value in values)
            for values in BOUND_NUMERATORS[target_id])
    points_core = tuple(tuple(primary.to_core(value) for value in point) for point in points)
    # Vertex inequalities imply whole convex-polygon bounds: the height
    # horoball and the centered geodesic half-plane |z|<=U are convex.
    for i in range(6):
        assert m[i] > 0 and U[i] > 0
        for _, y, modulus in points_core:
            assert core.sign_at(core.sub(y, core.scale(m[i]*m[i], core.ONE)), i) == 1
            assert core.sign_at(core.sub(core.scale(U[i]*U[i], core.ONE), modulus), i) == 1
    def record(value):
        return [str(q) for q in primary.to_core(value)]
    descriptor = {
        'schema_version': 1, 'target_id': target_id, 'embedding_order': list(EMBEDDINGS),
        'reflection_matrix': [record(v) for v in reflection] if reflection is not None else None,
        'reflection_determinant': record(determinant) if determinant is not None else None,
        'reflection_determinant_signs': list(determinant_signs) if determinant_signs is not None else None,
        'vertices': [{'name': name, 'center': record(x), 'height_squared': record(y),
                      'modulus_squared': record(modulus)}
                     for name, (x, y, modulus) in zip(('A', 'B', 'C', 'D'), points)],
        'sides': [{'center': list(map(str, beta)), 'radius_squared': list(map(str, radius)),
                   'interior_signs': list(sign)} for beta, radius, sign in sides],
        'm_lower': list(map(str, m)), 'U_upper': list(map(str, U)),
    }
    return descriptor, points_core, sides, m, U


def install_target(primary, target_id):
    """Install all target globals before fresh source bounds are constructed."""
    descriptor, points, sides, m, U = derive_target(primary, target_id)
    core = primary.core
    core.VERTEX_INVARIANTS = points
    core.SIDES = sides
    core.M_LOWER = m
    core.U_UPPER = U
    core.TARGET_ID = target_id
    # Fail closed if a caller attempts enumeration before reconstruct().
    core.Y2 = None
    for name in ("DEN", "NUM", "A", "B", "C", "D"):
        setattr(core, name+"_BOUND", None)
    return descriptor


def positive_boundary_controls(primary):
    """Actual classifier regression; no candidate search or source equivalence claim."""
    a, sp, core = primary.ambient, primary.sp, primary.core
    t = a.T
    delta = a.SIGN_UNITS[15]
    v = 1+5*t+8*t*t+2*t**3-4*t**4-2*t**5
    assert a.reduce_field(v*v-delta/(t*(4-t*t))) == 0
    epsilon = (a.reduce_l(2*v), a.reduce_l(-t*v))
    norm = a.reduce_l(epsilon[0]**2+t*epsilon[0]*epsilon[1]+epsilon[1]**2)
    assert a.reduce_field(norm-delta/t) == 0
    H = a.matrix_multiply(a.matrix_multiply(sp.diag(t, 1), a.multiplication_matrix(epsilon)),
                          sp.diag(a.reduce_field(1/delta), 1))
    expected = ((0,-8,0,2,0,0), (2,18,33,9,-16,-8),
                (0,4,0,-1,0,0), (0,-8,-17,-5,8,4))
    assert tuple(primary.to_core(value) for value in H) == expected
    inverse = a.matrix_reduce(H.inv())
    a.assert_integral_matrix(H)
    a.assert_integral_matrix(inverse)
    assert a.reduce_l(H.det()) == 1
    Z = sp.Matrix([[0, -1], [1, t]])
    raw = a.matrix_multiply(a.matrix_multiply(sp.diag(delta, 1), Z),
                            sp.diag(a.reduce_field(1/delta), 1))
    B = matrix_from_core(primary, core.VERTICES[1])
    assert a.matrix_multiply(a.matrix_multiply(H, raw), inverse) == B
    raw_point = vertex_invariants(primary, raw)
    old_y = core.Y2
    try:
        identity = (core.ONE, core.ZERO, core.ZERO, core.ONE)
        accepted_vertices = []
        for name, (x, y, _) in zip(('A', 'B'), core.VERTEX_INVARIANTS[:2]):
            core.Y2 = y
            assert core.in_all_coordinate_polygons(identity, x) == (True, None)
            accepted_vertices.append(name)
        core.Y2 = primary.to_core(raw_point[1])
        assert core.in_all_coordinate_polygons(expected, primary.to_core(raw_point[0])) == (True, None)
        real, modulus = core.image_invariants(expected, primary.to_core(raw_point[0]))
        assert (real, modulus) == (core.VERTEX_INVARIANTS[1][0], core.VERTEX_INVARIANTS[1][2])
    finally:
        core.Y2 = old_y
    return {'target_id': core.TARGET_ID, 'fixed_vertices_accepted': accepted_vertices,
            'C39_B_source': 'OM:pic0:15', 'C39_B_bridge': [list(map(str, row)) for row in expected],
            'C39_B_boundary_accepted': True}
