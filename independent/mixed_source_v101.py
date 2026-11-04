#!/usr/bin/env python3
"""Direct independent verification of V101 Fermat source witnesses.

Does not import the worker, enumerate partitions, or rerun Graver search.
The complete Graver list is consumed only after checking its accepted
independent-review digest. Missing/negative/timeout results remain unresolved.
Usage: python3 independent/mixed_source_v101.py OUTPUT_DIR
       python3 independent/mixed_source_v101.py --manifest-only
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import gcd
from pathlib import Path
import resource
import signal

ROOT = Path(__file__).resolve().parents[1]
G = tuple(a for a in range(1, 42) if gcd(a, 42) == 1)
R = (1, 5, 11, 13, 17, 19)
TYPES = {'B': frozenset((1, 5, 11, 13, 17, 19)),
         'D': frozenset((1, 5, 11, 13, 19, 25))}
FACTORS = ('B', 'D')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def file_hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def vector_multiset(z):
    require(len(z) == 12 and all(type(a) is int for a in z), 'bad vector')
    return Counter({(FACTORS[i // 6], R[i % 6] if a > 0 else 42-R[i % 6]): abs(a)
                    for i, a in enumerate(z) if a})


def multiset_vector(counts):
    out = []
    for f in FACTORS:
        for a in R:
            require(not (counts[f, a] and counts[f, 42-a]), 'opposite signs coexist')
            out.append(counts[f, a]-counts[f, 42-a])
    return tuple(out)


def translated(z, h):
    return multiset_vector(Counter({(f, a*h % 42): n
                                    for (f, a), n in vector_multiset(z).items()}))


def canonical_labels(z):
    return sorted((f, copy, a) for (f, a), n in vector_multiset(z).items()
                  for copy in range(n))


def balanced(entries):
    return all(sum((h*a) % 42 for a in entries) == 21*len(entries) for h in G)


def triple_type(t):
    require(len(t) == 3 and all(type(a) is int and 0 < a < 42 for a in t), 'bad triple residues')
    require(sum(t) % 42 == 0, 'triple not in A_42^1')
    return frozenset(h for h in G if sum(h*a % 42 for a in t) == 42)


def exact_presentations(f):
    # Each sorted residue triple occurs once: choose its first two entries,
    # recover the third from the zero-sum equation, and retain sorted order.
    found = set()
    for a in range(1, 42):
        for b in range(a, 42):
            c = (-a-b) % 42
            if b <= c < 42 and triple_type((a, b, c)) == TYPES[f]:
                found.add((a, b, c))
    return found


def check_case(case, copy_counts):
    key = []
    for f in FACTORS:
        ts = case[f+'_copy_triples']
        require(len(ts) == copy_counts[f], 'copy count differs')
        for t in ts:
            require(triple_type(t) == TYPES[f], 'exact marked CM type differs')
            conductor = 42 // gcd(42, gcd(t[0], gcd(t[1], t[2])))
            phi = sum(gcd(i, conductor) == 1 for i in range(1, conductor))
            require(phi == 12, 'factor does not have dimension six')
            require(len({tuple(h*a % 42 for a in t) for h in G}) == 12,
                    'character orbit is not twelve-dimensional in H1')
        key.append(tuple(map(tuple, ts)))
    return tuple(key)


def verify_manifest(path):
    manifest = json.loads(path.read_bytes())
    require(manifest['modulus'] == 42 and manifest['units'] == list(G), 'wrong modulus/units')
    require(manifest['column_order'] == [[f, a] for f in FACTORS for a in R], 'wrong column marks')
    worker = ROOT/'experiments/mixed_powers_v101/source_worker.py'
    require(file_hash(worker) == manifest['worker_sha256'], 'worker source pin mismatch')
    for name, pin in manifest['input_sha256'].items():
        require(file_hash(ROOT/name) == pin, 'input pin mismatch: '+name)
    review = json.loads((ROOT/'experiments/mixed_powers_v101/graver_review.json').read_bytes())
    kernel_path = ROOT/'experiments/mixed_powers_v101/kernel_box.json'
    require(file_hash(kernel_path) == review['input_sha256'][str(kernel_path.relative_to(ROOT))],
            'kernel differs from reviewed object')
    all_g = sorted(map(tuple, json.loads(kernel_path.read_bytes())['graver_vectors_up_to_sign']))
    digest = sha256(json.dumps(all_g, separators=(',', ':')).encode()).hexdigest()
    require(digest == review['independent_output_summary']['graver_digest'],
            'Graver list differs from independently reviewed digest')
    signed = {tuple(sign*a for a in z) for z in all_g if max(map(abs, z)) == 2 for sign in (-1, 1)}
    require(len(signed) == 168, 'wrong repeated-pattern domain')
    presentations = {f: exact_presentations(f) for f in FACTORS}
    require(len(presentations['B']) == 12 and presentations['D'] == {(1, 14, 27)},
            'unexpected exact-type presentation domain')
    require(manifest['task_count'] == len(manifest['tasks']) == 14, 'wrong task count')
    require(manifest['signed_repeated_vector_count'] == len(signed), 'wrong coverage declaration')
    covered = set(); case_count = 0
    for i, task in enumerate(manifest['tasks']):
        require(task['task_index'] == i, 'noncanonical task index')
        z = tuple(task['representative_vector'])
        require(z in signed, 'representative outside reviewed Graver list')
        labels = canonical_labels(z)
        require(sorted(map(tuple, task['embedding_labels'])) == labels, 'wrong embedding multiplicity/copy placement')
        require(len(labels) == len(set(labels)) == task['degree'], 'exterior label repeated or degree wrong')
        counts = {f: max((c+1 for ff, c, a in labels if ff == f), default=0) for f in FACTORS}
        require(counts == task['source_copy_counts'], 'source copy count wrong')
        # Direct Hodge balance from marked types, independently of weight matrices.
        require(all(sum(1 if h*a % 42 in TYPES[f] else -1 for f, c, a in labels) == 0 for h in G),
                'representative fails synchronized Hodge balance')
        orbit = {translated(z, h) for h in G}
        require(len(orbit) == 12 and orbit == set(map(tuple, task['signed_Galois_orbit'])), 'wrong orbit')
        require(not covered.intersection(orbit), 'duplicate orbit task')
        covered.update(orbit)
        actual = set()
        for j, case in enumerate(task['cases']):
            require(case['case_index'] == j, 'noncanonical case index')
            actual.add(check_case(case, counts)); case_count += 1
        expected = {(b, d) for b in product(sorted(presentations['B']), repeat=counts['B'])
                    for d in product(sorted(presentations['D']), repeat=counts['D'])}
        require(actual == expected and len(actual) == len(task['cases']), 'case list coverage/duplicate mismatch')
    require(covered == signed, 'manifest does not cover exactly all168 signed repeated generators')
    return manifest, {'signed_pattern_count': len(covered), 'orbit_count': 14,
                      'case_count': case_count, 'exact_B_presentations': sorted(presentations['B']),
                      'exact_D_presentations': sorted(presentations['D']),
                      'manifest_sha256': file_hash(path), 'worker_sha256': file_hash(worker)}


def check_block(block):
    entries = block['entries']
    require(len(entries) in (2, 4, 6), 'unsupported block length')
    require(all(type(a) is int and 0 < a < 42 for a in entries), 'bad block residues')
    require(balanced(entries), 'block is not Galois-balanced')
    positions = block.get('triple_positions')
    if len(entries) == 6:
        require(isinstance(positions, list) and len(positions) == 3
                and all(type(i) is int and 0 <= i < 6 for i in positions)
                and len(set(positions)) == 3, 'bad sextuple split positions')
        left = [entries[i] for i in positions]
        right = [entries[i] for i in range(6) if i not in positions]
        require(sum(left) % 42 == sum(right) % 42 == 0, 'sextuple not in A1*A1')
    else:
        require(positions is None, 'unexpected triple split on nonsix block')


def verify_result(path, manifest, task):
    result = json.loads(path.read_bytes())
    require(result['manifest_sha256'] == manifest['_hash'], 'returned manifest pin differs')
    require(result['worker_sha256'] == manifest['worker_sha256'], 'returned worker pin differs')
    for name in ('task_index', 'representative_vector', 'embedding_labels', 'source_copy_counts'):
        require(result[name] == task[name], 'result/task mismatch: '+name)
    require(result['case_count'] == len(task['cases']), 'returned case count wrong')
    found = []; statuses = Counter()
    for j, attempt in enumerate(result['attempts']):
        require(j < len(task['cases']), 'too many attempts')
        case = task['cases'][j]
        for key, value in case.items():
            require(attempt[key] == value, 'returned case metadata differs')
        triples = {f: case[f+'_copy_triples'] for f in FACTORS}
        gamma = sorted(g*a % 42 for f, c, g in task['embedding_labels'] for a in triples[f][c])
        require(gamma == attempt['gamma'] and balanced(gamma), 'gamma mismatch or imbalance')
        status = attempt['status']; statuses[status] += 1
        require(status in ('FOUND', 'INCOMPLETE', 'NO_PARTITION_THIS_PRESENTATION'), 'unknown attempt status')
        if status != 'FOUND':
            require(attempt.get('blocks') is None, 'nonpositive result has unexpected blocks')
            continue
        blocks = attempt['blocks']
        require(isinstance(blocks, list) and blocks, 'missing positive witness')
        require(Counter(a for b in blocks for a in b['entries']) == Counter(gamma), 'partition multiplicity mismatch')
        for block in blocks:
            check_block(block)
        # Validate all translates including copy labels, not just claimed orbit size.
        for h in G:
            labels = sorted((f, c, g*h % 42) for f, c, g in task['embedding_labels'])
            require(labels == canonical_labels(translated(task['representative_vector'], h)),
                    'Galois transfer changes canonical copy assignment')
            translated_gamma = Counter(g*a % 42 for f, c, g in labels for a in triples[f][c])
            require(translated_gamma == Counter(h*a % 42 for a in gamma), 'translated presentation mismatch')
            for block in blocks:
                check_block(dict(block, entries=[h*a % 42 for a in block['entries']]))
        found.append(j)
    if found:
        require(len(found) == 1 and found[0] == len(result['attempts'])-1, 'positive result did not stop once')
        require(result['terminal'] is True and result['status'] == 'POSITIVE_WITNESS_REQUIRES_INDEPENDENT_VALIDATION',
                'positive witness terminal status mismatch')
        require(result['witness_case_index'] == found[0], 'wrong witness case index')
    else:
        require(result['status'] != 'POSITIVE_WITNESS_REQUIRES_INDEPENDENT_VALIDATION', 'false positive declaration')
    return {'task_index': task['task_index'], 'file_sha256': file_hash(path),
            'positive_witness_validated': bool(found), 'witness_case_index': found[0] if found else None,
            'attempt_status_counts': dict(statuses), 'reported_status': result['status'],
            'reported_terminal': result['terminal'],
            'negative_or_incomplete_search_claims_independently_verified': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir', type=Path, nargs='?')
    parser.add_argument('--manifest', type=Path, default=ROOT/'experiments/mixed_powers_v101/job_manifest.json')
    parser.add_argument('--manifest-only', action='store_true')
    args = parser.parse_args()
    require(args.manifest_only or args.output_dir is not None, 'output directory required')
    resource.setrlimit(resource.RLIMIT_CPU, (18, 19)); signal.alarm(30)
    manifest, audit = verify_manifest(args.manifest)
    manifest['_hash'] = file_hash(args.manifest)
    report = {'status': 'SCIENTIFIC_MANIFEST_VALIDATED', 'manifest_audit': audit,
              'validator_sha256': file_hash(__file__), 'results': [], 'missing_tasks': [],
              'marked_isogeny_argument': 'Each source slot has exactly the prescribed marked CM type over M, even when the two B slots use different Fermat triples. Rank-one M rational H1 spaces of identical type are M-linearly Hodge-isomorphic; full lattices are commensurable and clearing denominators gives an M-linear isogeny on each slot. The product isogeny preserves factor, embedding, and copy labels up to nonzero scalars. Pullback/pushforward preserves rational algebraicity. No polarization equality is asserted.',
              'primary_source_scope': 'Shioda1981 Lemmas4.1-4.2 and the supportwise proof of Theorem4.3, printed pp70-72, as pinned in second_primitive_scope.json. A balanced sextuple alone is insufficient: an actual split into two zero-sum nonzero triples is required and checked. Individual split eigenlines mean membership in the extension of rational algebraic classes.',
              'scope': 'Positive witness verification only, with exact168-pattern orbit coverage. Does not rerun partition search or infer nonalgebraicity from absence/failure. All-powers consequences require the separately reviewed Graver reduction and rational descent.'}
    if not args.manifest_only:
        require(args.output_dir.is_dir(), 'output directory does not exist')
        expected = {f'task-{i:02d}.json' for i in range(14)}
        require({p.name for p in args.output_dir.glob('task-*.json')} <= expected, 'unexpected task result file')
        for task in manifest['tasks']:
            path = args.output_dir/f"task-{task['task_index']:02d}.json"
            if path.exists():
                report['results'].append(verify_result(path, manifest, task))
            else:
                report['missing_tasks'].append(task['task_index'])
        positive = [r['task_index'] for r in report['results'] if r['positive_witness_validated']]
        report['validated_positive_orbits'] = positive
        report['validated_signed_patterns'] = 12*len(positive)
        report['unresolved_signed_patterns'] = 168-12*len(positive)
        report['status'] = 'ALL_168_SIGNED_SOURCE_WITNESSES_VALIDATED' if len(positive) == 14 else 'PARTIAL_POSITIVE_WITNESSES_VALIDATED_REMAINDER_UNRESOLVED'
    print(json.dumps(report, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
