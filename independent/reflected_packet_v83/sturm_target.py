"""Independent exact target chambers; no imports from the primary arithmetic.

The reflected target changes the searched region, never the source's SL2 group.
All arithmetic comes from the retained standard-library field/Sturm foundation.
"""
from fractions import Fraction as Q

TARGET_IDS = ("ORIGINAL_Q_V1", "REFLECTED_AB_V1")
EMBEDDINGS = (1, 5, 11, 13, 17, 19)
BOUNDS = {
    "ORIGINAL_Q_V1": ((4006, 13933, 8395, 17094, 20761, 6663),
                      (255285, 208594, 17707, 121724, 165248, 191115)),
    "REFLECTED_AB_V1": ((1368, 23817, 14522, 22806, 10088, 5105),
                        (197767, 354283, 28993, 254603, 165248, 191115)),
}


def coefficients(value):
    return [str(Q(q)) for q in value.coefficients]


def matrix_record(value):
    return [coefficients(entry) for row in value for entry in row]


def point_invariants(field, action):
    x, y = field.exact_fixed_invariants(action)
    return x, y, x*x+y


def anti_image(field, reflection, point):
    """Invariant formulas for R conjugate(z), including det(R)^2 in Y'."""
    x, y, _modulus = point
    a, b, c, d = (entry for row in reflection for entry in row)
    determinant = field.mdet(reflection)
    top, bottom = a*x+b, c*x+d
    denominator = bottom*bottom+c*c*y
    center = (top*bottom+a*c*y)/denominator
    height = determinant*determinant*y/(denominator*denominator)
    modulus = (top*top+a*a*y)/denominator
    assert modulus == center*center+height
    return (center, height, modulus), denominator


def side_equation(point, side):
    x, _y, modulus = point
    return modulus-2*side.center*x+side.center*side.center-side.radius_squared


def install_target(retained, target_id):
    """Rebuild target data and invalidate all point-dependent cached bounds."""
    if target_id not in TARGET_IDS:
        raise ValueError("Unknown exact target chamber")
    field, orbit = retained.ambient, retained.orbit
    assert field.EMBEDDINGS == orbit.EMBEDDINGS == EMBEDDINGS
    # Reconstruct from the retained printed matrices each time, not mutable globals.
    original_actions = orbit.normalized_vertices()
    original = tuple(point_invariants(field, action) for action in original_actions)
    assert all(orbit.signs(y) == (1,)*6 for _x, y, _m in original)
    reflection = determinant = None
    determinant_signs = None
    if target_id == "REFLECTED_AB_V1":
        a, b = original_actions[:2]
        reflection = field.msub(field.mmul(a, b), field.mmul(b, a))
        determinant = field.mdet(reflection)
        determinant_signs = orbit.signs(determinant)
        assert determinant_signs == (-1,)*6
        assert determinant.norm() == 21952
        assert field.mtrace(reflection) == field.ZERO
        assert field.mmul(reflection, reflection) == field.mscale(-determinant, field.IDENTITY)
        for vertex in (a, b):
            assert field.mmul(reflection, vertex) == field.mmul(field.minverse(vertex), reflection)
        derived = tuple(anti_image(field, reflection, point) for point in original)
        assert all(orbit.signs(denominator) == (1,)*6 for _, denominator in derived)
        points = tuple(point for point, _ in derived)
        actions = tuple(field.conjugate_by(reflection, action) for action in original_actions)
        assert points == tuple(point_invariants(field, action) for action in actions)
        assert points[:2] == original[:2]
        assert tuple(anti_image(field, reflection, point)[0] for point in points) == original
        # These are explicit vertex-orbit relations, not a full packet permutation.
        assert points[2] == point_invariants(field, field.conjugate_by(b, original_actions[2]))
        assert points[3] == point_invariants(field, field.conjugate_by(field.minverse(a), original_actions[3]))
    else:
        actions, points = original_actions, original
    assert all(orbit.signs(y) == orbit.signs(modulus) == (1,)*6 for _x, y, modulus in points)
    sides = []
    for index in range(4):
        x0, y0, m0 = points[index]
        x1, _y1, m1 = points[(index+1) % 4]
        assert x1 != x0  # These actual two targets have no vertical side.
        center = (m1-m0)/(2*(x1-x0))
        radius = (x0-center)**2+y0
        assert orbit.signs(radius) == (1,)*6
        trial = orbit.Side(center, radius, (0,)*6)
        assert side_equation(points[index], trial) == field.ZERO
        assert side_equation(points[(index+1) % 4], trial) == field.ZERO
        interior = orbit.signs(side_equation(points[(index+2) % 4], trial))
        assert 0 not in interior
        assert interior == orbit.signs(side_equation(points[(index+3) % 4], trial))
        sides.append(orbit.Side(center, radius, interior))
    m_numerators, u_numerators = BOUNDS[target_id]
    lower = tuple(Q(value, 100000) for value in m_numerators)
    upper = tuple(Q(value, 100000) for value in u_numerators)
    assert all(value > 0 for value in lower+upper)
    for i, root in enumerate(orbit.ROOTS):
        for _x, height, modulus in points:
            assert lower[i] < orbit.sqrt_interval(orbit.interval_eval(height, root), 160).lo
            assert upper[i] > orbit.sqrt_interval(orbit.interval_eval(modulus, root), 160).hi
    # Replace target globals, then invalidate caches until reconstruct()/install_point_data.
    orbit.VERTICES = actions
    orbit.VERTEX_DATA = points
    orbit.SIDES = tuple(sides)
    orbit.M_LOWER, orbit.U_UPPER = lower, upper
    orbit.TARGET_ID = target_id
    for name in ("Y2", "Y2_SCALED", "DEN_RADIUS", "NUM_RADIUS"):
        setattr(orbit, name, None)
    for name in ("DEN", "NUM", "A", "B", "C", "D"):
        setattr(orbit, name+"_BOUND", None)
    return {
        "schema_version": 1, "target_id": target_id,
        "embedding_order": list(EMBEDDINGS),
        "reflection_matrix": None if reflection is None else matrix_record(reflection),
        "reflection_determinant": None if determinant is None else coefficients(determinant),
        "reflection_determinant_signs": None if determinant_signs is None else list(determinant_signs),
        "vertices": [{"name": name, "center": coefficients(x), "height_squared": coefficients(y),
                      "modulus_squared": coefficients(modulus)}
                     for name, (x, y, modulus) in zip(("A", "B", "C", "D"), points)],
        "sides": [{"center": coefficients(side.center), "radius_squared": coefficients(side.radius_squared),
                   "interior_signs": list(side.interior_signs)} for side in sides],
        "m_lower": list(map(str, lower)), "U_upper": list(map(str, upper)),
    }
