#!/usr/bin/env python3
"""Reconstruct cyclotomic weights and verify a displayed integral basis.

No primary pivot chart, Graver search, or HNF implementation is imported.
"""
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def check():
    names=['experiments/mixed_powers_v101/lattice_generation.json',
           'experiments/fixed_field_cm_census_v98/independent.json']
    data,census=[json.loads((ROOT/n).read_text()) for n in names]
    units=[a for a in range(42) if s.gcd(a,42)==1]
    reps=[1,5,11,13,17,19]
    types=[set(census['orbits'][i]['representative']) for i in (0,2)]
    W=s.Matrix([[1 if h*a%42 in phi else -1 for phi in types for a in reps] for h in units])
    B=s.Matrix(data['basis_vectors']).T
    free=data['free_columns']
    P=B[free,:]
    assert W.rank()==6 and B.rank()==6 and W*B==s.zeros(12,6)
    assert P.det()==1
    inverse=P.inv()
    assert all(x.q==1 for x in inverse)
    assert inverse.tolist()==data['integer_inverse']
    assert [sum(abs(x) for x in B[:,j]) for j in range(6)]==[6,4,4,4,4,4]
    assert list(B[:6,0])==[0]*6
    assert all(abs(x)<=1 for x in B)
    distinct_counts=[]
    for start in (0,6):
        directions=[tuple(W[:,i]) for i in range(start,start+6)]
        all_weights=directions+[tuple(-x for x in v) for v in directions]
        assert len(set(all_weights))==12 and all(any(v) for v in all_weights)
        distinct_counts.append(len(set(all_weights)))
    supports=[]
    for j in range(6):
        supports.append([[('B' if i<6 else 'D'),reps[i%6] if B[i,j]>0 else 42-reps[i%6]]
                         for i in range(12) if B[i,j]])
    return {'status':'INDEPENDENT_INTEGER_LATTICE_BASIS_VERIFIED',
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        'weight_rank':6,'basis_rank':6,'projected_basis_determinant':1,
        'distinct_nonzero_weights_per_factor':distinct_counts,
        'integer_inverse':inverse.tolist(),'basis_embedding_supports':supports,
        'argument':'The displayed columns lie in the rank-six rational kernel. Their projection to the six named coordinates is unimodular. Every integer kernel vector therefore has an integer coefficient vector in this basis; subtracting it leaves a kernel vector with zero free projection and hence zero. No saturation is inferred from rank alone.',
        'algebraicity_inputs':'First support is pure-D degree6 covered by C45; all five other supports are invariant squarefree quartets on B×D covered by C47. Algebraic contraction is a separate proof.'}

if __name__=='__main__':
    import resource,signal
    resource.setrlimit(resource.RLIMIT_CPU,(20,22));signal.alarm(30)
    print(json.dumps(check(),indent=2,sort_keys=True,default=int))
