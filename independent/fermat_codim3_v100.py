"""Verify supplied Shioda codimension-three witnesses; never search partitions.

Only the input atom list's sixty mixed degree-six supports are covered here.
Completeness of the atom list is a separate proof obligation.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import math

ROOT = Path(__file__).resolve().parents[1]
G = tuple(g for g in range(42) if math.gcd(g, 42) == 1)
TYPES = {"B": {1, 5, 11, 13, 17, 19}, "D": {1, 5, 11, 13, 19, 25}}


def shift_support(support, g):
    return tuple(sorted((f, a*g % 42) for f, a in support))


def verify_block(block, triple_positions):
    assert len(block) in (2, 4, 6)
    assert all(0 < a < 42 for a in block)
    sums = [sum(a*g % 42 for a in block) for g in G]
    assert sums == [21*len(block)]*12
    if len(block) == 6:
        assert triple_positions is not None
        assert len(triple_positions) == len(set(triple_positions)) == 3
        assert set(triple_positions) <= set(range(6))
        left = [block[i] for i in triple_positions]
        right = [a for i, a in enumerate(block) if i not in triple_positions]
        assert sum(left) % 42 == sum(right) % 42 == 0
    else:
        assert triple_positions is None
    return sums


def check(probe_path=None, atom_path=None):
    probe_path = Path(probe_path or ROOT/'experiments/mixed_products_v100/codim3_source_probe.json')
    atom_path = Path(atom_path or ROOT/'experiments/mixed_products_v100/primitive_pair_atoms.json')
    probe, atoms = json.loads(probe_path.read_text()), json.loads(atom_path.read_text())
    assert probe['input_sha256'][str(atom_path.relative_to(ROOT))] == sha256(atom_path.read_bytes()).hexdigest()
    target = {tuple(sorted(tuple(x) for x in atom['support'])) for atom in atoms['atoms']
              if len(atom['support']) == 6 and {x[0] for x in atom['support']} == {'B','D'}}
    seen, orbit_sizes, presentations, block_sizes = set(), [], Counter(), Counter()
    checks = []
    for entry in probe['orbits']:
        rep = tuple(sorted(tuple(x) for x in entry['representative']))
        assert len(set(rep)) == 6 and rep in target
        successful = [a for a in entry['attempts'] if a['blocks'] is not None]
        assert entry['success'] and len(successful) == 1
        cert = successful[0]
        triples = {'B': tuple(cert['B_triple']), 'D': tuple(cert['D_triple'])}
        for f, triple in triples.items():
            assert len(triple) == 3 and all(0 < a < 42 for a in triple)
            assert sum(triple) % 42 == 0
            phi = {g for g in G if sum(a*g % 42 for a in triple) == 42}
            assert phi == TYPES[f]
            assert len({tuple(a*g % 42 for a in triple) for g in G}) == 12
        orbit = {shift_support(rep, g) for g in G}
        assert orbit == {tuple(sorted(tuple(x) for x in s)) for s in entry['orbit']}
        assert orbit <= target and not seen.intersection(orbit)
        seen.update(orbit)
        orbit_sizes.append(len(orbit))
        presentations[str(triples['B'])] += 1
        gamma = [a*label % 42 for f, label in rep for a in triples[f]]
        assert cert['gamma'] == sorted(gamma)
        assert Counter(gamma) == Counter(a for b in cert['blocks'] for a in b)
        assert len(cert['blocks']) == len(cert['block_checks'])
        for block, witness in zip(cert['blocks'], cert['block_checks']):
            assert witness['block'] == block
            sums = verify_block(block, witness['triple_positions'])
            assert sums == witness['all_unit_sums']
            block_sizes[len(block)] += 1
        for g in G:
            translated_gamma = [a*label % 42 for f, label in shift_support(rep, g) for a in triples[f]]
            translated_blocks = [[a*g % 42 for a in block] for block in cert['blocks']]
            assert Counter(translated_gamma) == Counter(a for b in translated_blocks for a in b)
            for block, witness in zip(translated_blocks, cert['block_checks']):
                verify_block(block, witness['triple_positions'])
        checks.append({'representative': rep, 'orbit_size':len(orbit),
                       'B_triple':triples['B'], 'D_triple':triples['D'],
                       'block_lengths':[len(b) for b in cert['blocks']]})
    assert seen == target and len(target) == 60
    assert orbit_sizes == [12, 6, 12, 12, 6, 12]
    return {'status':'PASS_ALL_SUPPLIED_MIXED_DEGREE6_ATOM_WITNESSES',
            'mixed_support_count':len(target), 'Galois_orbit_sizes':orbit_sizes,
            'presentations_by_orbit':dict(presentations),
            'representative_block_counts_by_length':dict(block_sizes),
            'orbit_checks':checks,
            'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()
                            for p in (probe_path,atom_path)},
            'atom_completeness_independently_established_here':False,
            'all_powers_inference':False}


if __name__ == '__main__':
    import resource
    resource.setrlimit(resource.RLIMIT_CPU, (18, 19))
    assert __debug__
    print(json.dumps(check(), sort_keys=True, separators=(',', ':')))
