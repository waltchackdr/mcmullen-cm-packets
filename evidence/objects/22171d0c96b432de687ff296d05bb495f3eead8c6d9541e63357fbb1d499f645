"""Exact small-prime input for the separate minimum-polarization human proof.

No numerical roots, external factorization library, or global-unit inference.
Coefficient vectors use increasing powers. The finite search is all rational
primes below 27, with all possible small residue degrees accounted for.
"""
from fractions import Fraction
from itertools import product
import json

F = [1, 8, 8, -6, -6, 1, 1]
BETA = [-2, -8, -2, 6, 1, -1]


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def remainder(a, b, p):
    a = trim([x % p for x in a])
    while len(a) >= len(b) and a != [0]:
        c = a[-1] * pow(b[-1], -1, p) % p
        shift = len(a) - len(b)
        for j, x in enumerate(b):
            a[j + shift] = (a[j + shift] - c * x) % p
        a = trim(a)
    return a


def multiply(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def evaluate(a, x, p=None):
    v = 0
    for c in reversed(a):
        v = v * x + c
        if p:
            v %= p
    return v


def norm(a):
    columns = []
    for j in range(6):
        b = [0] * j + list(a)
        b += [0] * max(0, 6 - len(b))
        for k in range(len(b) - 1, 5, -1):
            for i, c in enumerate(F):
                b[k - 6 + i] -= b[k] * c if i < 6 else 0
            b[k] = 0
        columns.append(b[:6])
    m = [[Fraction(x) for x in row] for row in zip(*columns)]
    answer = Fraction(1)
    for i in range(6):
        pivots = [j for j in range(i, 6) if m[j][i]]
        if not pivots:
            return 0
        j = pivots[0]
        if i != j:
            m[i], m[j] = m[j], m[i]
            answer = -answer
        v = m[i][i]
        answer *= v
        for j in range(i + 1, 6):
            q = m[j][i] / v
            for k in range(i, 6):
                m[j][k] -= q * m[i][k]
    assert answer.denominator == 1
    return int(answer)


def check():
    if not __debug__:
        raise RuntimeError('Enabled assertions required')
    primes = [n for n in range(2, 27)
              if all(n % d for d in range(2, n))]
    assert primes == [2, 3, 5, 7, 11, 13, 17, 19, 23]
    # Any reducible degree-six polynomial has a monic factor of degree <=3.
    trials = 0
    for d in range(1, 4):
        for coeffs in product(range(2), repeat=d):
            assert remainder(F, list(coeffs) + [1], 2) != [0]
            trials += 1
    assert trials == 14
    factors = {3: [[1, 1, -1, 1], [1, 1, -1, 1]],
               5: [[-2, -1, 0, 1], [2, 0, 1, 1]],
               7: [[-1, 1]] * 6}
    for p, fs in factors.items():
        prod = [1]
        for q in fs:
            prod = multiply(prod, q)
            if len(q) == 4:
                assert all(evaluate(q, x, p) for x in range(p))
        assert [x % p for x in prod] == [x % p for x in F]
    # For p>=11, p^2>27: only linear factors can give norm below27.
    no_roots = {}
    for p in primes[4:]:
        no_roots[str(p)] = [evaluate(F, x, p) for x in range(p)]
        assert all(no_roots[str(p)])
    assert norm([-1, 1]) == evaluate(F, 1) == 7
    assert norm(BETA) == -27
    return {'status': 'EXACT_SMALL_PRIME_INPUT_PASS', 'f': F,
            'primes_exhausted_below27': primes,
            'mod2_monic_divisor_trials': trials,
            'factorizations_mod3_5_7': factors,
            'nonzero_root_evaluations': no_roots,
            'only_prime_ideal_norm_below27': 7,
            'unique_prime_above3_norm': 27,
            'norm_t_minus1': 7, 'norm_beta': -27,
            'scope': 'Prime-ideal conclusions use the accepted O_L=Z[t]. '
                     'No unit-sign, positivity, or geometric theorem proved by this loop.'}


if __name__ == '__main__':
    import resource
    import signal
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    signal.alarm(30)
    print(json.dumps(check(), indent=2, sort_keys=True))
