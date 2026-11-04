"""Finite inputs for Pic(R)=Z/2: Frobenius/GCD and exact resultants.
The full theorem additionally requires the stated class-number, normalized
ideal classification and global CM-unit arguments. No code is run on import.
"""
import argparse, hashlib, json, resource, signal, sys
from pathlib import Path
import sympy as s

def check(directory=None):
    if sys.flags.optimize:
        raise RuntimeError('Assertions must be enabled')
    directory = Path(directory) if directory is not None else Path(__file__).resolve().parent
    pins = {
        'principal_bc_unpolarized_data.json': 'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e',
        'principal_bc_hom_data.json': '9f815f1592c8e22a2a14b15086d8d9cb63092f9ab84bb90a3d5ee4935b3fe22c',
    }
    data = {}
    for path, digest in pins.items():
        raw = (directory/path).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == digest
        data[path] = json.loads(raw)
    x = s.Symbol('x')
    phi = s.Poly(s.cyclotomic_poly(42, x), x, domain=s.QQ)
    h = s.Poly(x**6-x**5+x**4-x**3+x**2-x+1, x, modulus=3)
    assert s.Poly(phi.as_expr(), x, modulus=3) == h*h
    def fp(a):
        return s.rem(s.Poly(a, x, modulus=3), h)
    def mul(a, b):
        return s.rem(a*b, h)
    one, zero, z = fp(1), fp(0), fp(x)
    def power(a, n):
        b = one
        while n:
            if n & 1:
                b = mul(b, a)
            a = mul(a, a)
            n //= 2
        return b
    # Rabin criterion: only prime divisors 2 and 3 of degree 6 are needed.
    assert power(z, 729) == z
    for n in (9, 27):
        assert s.gcd(h, power(z, n)-z).degree() == 0
    assert power(z, 14) == one and power(z, 7) != one and power(z, 2) != one
    assert s.resultant(phi, s.Poly(1-x**14, x)) == 729
    # P lies in the residue kernel and both have index 729, hence are equal.
    assert power(z, 27) == power(z, 13)  # residue conjugation = Frobenius^3
    t = z+power(z, 13)
    assert power(t, 27) == t and power(t, 3) != t
    # Therefore OL's residue image contains the degree-3 fixed subfield F27.
    def ratio(a):
        assert a != zero
        return mul(a, s.invert(power(a, 27), h))
    gamma_q = s.Poly(x**6*(1-x**2), x, domain=s.QQ)
    assert s.resultant(phi, gamma_q) == 1  # global integral unit
    assert s.rem(gamma_q+s.Poly(x**14*(x**36-x**34), x), phi).is_zero
    assert s.rem(gamma_q-s.Poly(sum(s.Rational(c)*x**j for j,c in enumerate(data['principal_bc_hom_data.json']['square_generator'])), x), phi).is_zero
    gamma = fp(gamma_q.as_expr())
    alpha_coeff = data['principal_bc_unpolarized_data.json']['alpha']
    alpha_q = s.Poly(sum(s.Rational(c)*x**j for j,c in enumerate(alpha_coeff)), x)
    assert s.resultant(phi, alpha_q) == 64  # deliberately not a global unit
    alpha = fp(alpha_q.as_expr())
    sqrt7 = fp(sum(x**j for j in (6,12,24))-sum(x**j for j in (18,30,36)))
    assert mul(sqrt7, sqrt7) == fp(-1)
    assert ratio(alpha) == sqrt7
    rz, rg, ra = ratio(z), ratio(gamma), ratio(alpha)
    assert power(rz, 7) == one and rz != one
    assert rg == fp(-1) and rg != one
    assert power(ra, 4) == one and power(ra, 2) == fp(-1)
    assert power(ra, 14) != one
    # In cyclic F729*, Q has order 28. Ratio induces an injection on Q.
    # Global CM-unit ratios reduce into mu14. Explicit unit classes have
    # orders 7 and 2, hence generate all 14 of the possible unit-image classes.
    return {
        'status': 'PRINCIPAL_ORDER_PICARD_FINITE_INPUTS_PASS',
        'method': 'SymPy Frobenius/GCD irreducibility and rational resultants; no field enumeration',
        'input_pins': pins,
        'residue_field_size': 729, 'real_residue_field_size': 27,
        'conductor_generator_absolute_norm': 729,
        'quotient_order': 28, 'z_ratio_order': 7, 'gamma_ratio_order': 2,
        'gamma_absolute_norm': 1, 'alpha_absolute_norm': 64,
        'alpha_ratio_order': 4, 'alpha_outside_mu14': True,
        'conditional_unit_image_order': 14, 'conditional_Picard_order': 2,
        'human_dependencies': ['Pic(OM)=1', 'normalized invertible ideal classification', 'CM-unit ratios lie in mu42'],
    }

def main():
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    signal.alarm(30)
    ap = argparse.ArgumentParser()
    ap.add_argument('--certificates', type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args()
    print(json.dumps(check(args.certificates), indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
