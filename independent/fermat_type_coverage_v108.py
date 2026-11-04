#!/usr/bin/env python3
"""Ordered-character, fractional-age cross-check of CM-type coverage.

No Jacobian isogeny or Hodge algebraicity conclusion is asserted here.
"""
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path
import json
import resource
import signal

ROOT = Path(__file__).resolve().parents[1]
CENSUS = 'experiments/fixed_field_cm_census_v98/primary.json'
PIN = '2135edfe6a09b1a2c51a254f80e4570a761674e5be3208c15c549087f225917e'


def check():
    raw = (ROOT / CENSUS).read_bytes()
    if sha256(raw).hexdigest() != PIN:
        raise ValueError('Census digest mismatch')
    representatives = [set(r['representative']) for r in json.loads(raw)['orbits']]
    units = [g for g in range(1, 42) if gcd(g, 42) == 1]
    orbit_lookup = {}
    for index, representative in enumerate(representatives):
        for g in units:
            key = tuple(sorted(g * t % 42 for t in representative))
            if key in orbit_lookup and orbit_lookup[key] != index:
                raise ValueError('Overlapping census orbits')
            orbit_lookup[key] = index
    counts = [0] * 8
    conductors = [set() for _ in range(8)]
    actual_types, rows = set(), []
    # Parameterize ALL nonzero ordered characters, including the age-two half.
    for first in range(1, 42):
        for second in range(1, 42):
            third = (-(first + second)) % 42
            if not third:
                continue
            triple = (first, second, third)
            selected = []
            for g in units:
                fractions = [Fraction(g * x, 42) for x in triple]
                age = sum(x - int(x) for x in fractions)
                if age not in (1, 2):
                    raise ValueError('Unexpected character age')
                if age == 1:
                    selected.append(g)
            key = tuple(selected)
            if len(key) != 6 or any((42 - g) in key for g in key):
                raise ValueError('Not a CM type')
            index = orbit_lookup[key]
            conductor = 42 // gcd(42, gcd(first, gcd(second, third)))
            counts[index] += 1
            conductors[index].add(conductor)
            actual_types.add(key)
            rows.append({'character': list(triple), 'type': list(key),
                         'orbit_index': index, 'conductor': conductor})
    if len(rows) != 1640 or len(actual_types) != 64 or set(actual_types) != set(orbit_lookup):
        raise ValueError('Exact domain or coverage failed')
    if (ROOT / CENSUS).read_bytes() != raw:
        raise ValueError('Census changed during check')
    return {'status': 'EXACT_FINITE_CHARACTER_COVERAGE_CROSSCHECK_NOT_THEOREM_PROMOTION',
            'input_sha256': {CENSUS: PIN}, 'ordered_character_count': len(rows),
            'distinct_types': len(actual_types), 'orbit_counts': counts,
            'conductors_by_orbit': [sorted(c) for c in conductors],
            'rows': rows, 'new_hodge_algebraicity_claim': False,
            'jacobian_isogeny_interpretation_verified_by_this_code': False}


if __name__ == '__main__':
    resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
    signal.alarm(30)
    result = check()
    target = ROOT / 'experiments/fermat_type_coverage_v108/independent.json'
    target.parent.mkdir(exist_ok=True)
    with target.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('status', 'ordered_character_count', 'distinct_types', 'orbit_counts', 'conductors_by_orbit')}))
