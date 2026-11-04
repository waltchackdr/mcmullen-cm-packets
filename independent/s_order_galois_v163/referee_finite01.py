"""Separate finite audit in actual six-dimensional coordinates.

Enumerates semilinear maps by their image of a cyclic vector, using the
actual root matrix; no polynomial-field producer, RREF or generator search.
Subspaces are full 27-point sets. No producer modules are imported.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[2]
def sha(b):return hashlib.sha256(b).hexdigest()
def load(p):return json.loads((ROOT/p).read_bytes())
V=list(itertools.product(range(3),repeat=6))
CODE={v:i for i,v in enumerate(V)}
ZERO=(0,)*6
E=tuple(tuple(int(i==j) for i in range(6)) for j in range(6))
# All matrices below are tuples of columns, not rows.
def act(m,v):return tuple(sum(m[j][i]*v[j] for j in range(6))%3 for i in range(6))
def comp(m,n):return tuple(act(m,col) for col in n)
def power(m,n):
    out=E
    for _ in range(n):out=comp(m,out)
    return out
def span(rows):
    return frozenset(CODE[tuple(sum(c[j]*rows[j][i] for j in range(3))%3 for i in range(6))] for c in itertools.product(range(3),repeat=3))

def run():
    packet=load('experiments/s_order_galois_v163/review_inputs01.json')
    for p,h in packet['pins'].items():assert sha((ROOT/p).read_bytes())==h,p
    seed=load('experiments/s_order_models_v161/model.json')
    catalog=load('experiments/s_order_models_v161/catalog_candidate.json')['models']
    root=tuple(tuple(int(seed['order7_action'][i][j])%3 for i in range(6)) for j in range(6))
    omega=[[int(x)%3 for x in row] for row in seed['discriminant_gram']]
    pair=lambda x,y:sum(x[i]*omega[i][j]*y[j] for i in range(6) for j in range(6))%3
    assert power(root,7)==E and root!=E
    cyclic=tuple(act(power(root,i),E[0]) for i in range(6))
    decode={act(cyclic,v):v for v in V};assert len(decode)==729
    basis_coeffs=[decode[e] for e in E]
    # The complete scalar-times-Frobenius map is fixed by h(e0)=w and
    # h(R^i e0)=R^(3^j i)w. Test against the actual alternating form.
    group={};counts={};bad=0
    for j in range(6):
        rp=power(root,pow(3,j,7));operators=[power(rp,i) for i in range(6)]
        for w in V[1:]:
            images=tuple(act(op,w) for op in operators)
            m=tuple(act(images,c) for c in basis_coeffs)
            factors=[eps for eps in (1,2) if all(pair(m[a],m[b])==eps*omega[a][b]%3 for a in range(6) for b in range(a+1,6))]
            if not factors:bad+=1;continue
            assert len(factors)==1 and comp(m,root)==comp(rp,m)
            group[m]=(j,factors[0]);counts[j,factors[0]]=counts.get((j,factors[0]),0)+1
    assert len(group)==336 and all(counts[j,e]==28 for j in range(6) for e in (1,2))
    normal={power(root,i) for i in range(7)}
    normal|={tuple(tuple(-x%3 for x in col) for col in m) for m in list(normal)}
    assert len(normal)==14 and normal<=set(group)
    remaining=set(group);cosets=[]
    while remaining:
        m=min(remaining);coset={comp(n,m) for n in normal}
        assert len(coset)==14 and coset<=remaining
        remaining-=coset;cosets.append((m,coset))
    assert len(cosets)==24
    # Full point sets avoid the producer's polynomial coordinates and RREF.
    lookup={};representatives=[]
    for model in catalog:
        representatives.append(model['lagrangian_rows'])
        for rows in model['root_orbit']:
            h=span(rows);assert len(h)==27 and h not in lookup
            assert all(pair(x,y)==0 for x in rows for y in rows)
            lookup[h]=model['index']
    assert len(lookup)==1120
    # Precompute ternary addition solely to make point-set construction small.
    add=[[CODE[tuple((x+y)%3 for x,y in zip(a,b))] for b in V] for a in V]
    neg=[CODE[tuple(-x%3 for x in v)] for v in V]
    def smallspan(rows):
        a,b,c=[CODE[v] for v in rows]
        return frozenset(add[add[x][y]][z] for x in (0,a,neg[a]) for y in (0,b,neg[b]) for z in (0,c,neg[c]))
    permutations=[]
    for m,coset in cosets:
        perm=tuple(lookup[smallspan([act(m,tuple(row)) for row in rows])] for rows in representatives)
        assert sorted(perm)==list(range(160))
        assert all(catalog[i]['endomorphism_order']==catalog[k]['endomorphism_order'] for i,k in enumerate(perm))
        permutations.append(perm)
    assert len(set(permutations))==24
    perms=set(permutations)
    assert all(tuple(p[q[i]] for i in range(160)) in perms for p in perms for q in perms)
    unseen=set(range(160));orbits=[]
    while unseen:
        i=min(unseen);orb=sorted({p[i] for p in perms});assert set(orb)<=unseen
        unseen-=set(orb);orbits.append(orb)
    r=[o for o in orbits if catalog[o[0]]['endomorphism_order']=='R']
    ss=[o for o in orbits if catalog[o[0]]['endomorphism_order']=='S']
    assert len(r)==1 and r[0]==[8,38,118,153]
    assert sorted(map(len,ss))==[4,8,12,12,12,12,24,24,24,24]
    phi=seed['actual_tangent_type'];units=[x for x in range(42) if math.gcd(x,42)==1]
    stabilizer=[x for x in units if {x*a%42 for a in phi}==set(phi)]
    assert stabilizer==[1]
    inverses=[pow(a,-1,42) for a in phi]
    exponents=[next(j for j in range(6) if pow(3,j,7)==a%7) for a in inverses]
    residue_power=sum(3**j for j in exponents)
    assert exponents==[0,1,2,3,5,1] and residue_power==286
    # Enumerate the abstract cyclic image, including the root quotient, rather
    # than trusting a single generator or the displayed gcd claim.
    image={(286*i)%728 for i in range(728)};roots={(52*i)%728 for i in range(14)}
    assert len(image)==28 and roots<image and len({min((v-r)%728 for r in roots) for v in image})==2
    candidate=load('experiments/s_order_galois_v163/finite_action01.json')
    expected=candidate['orbits']['full_semilinear_similitude_candidate']
    assert sorted(orbits)==sorted(expected)
    return {'status':'PASS','method':'Actual-coordinate exhaustive semilinear maps by cyclic-vector image, direct pairing, complete point-set subspaces; no polynomial arithmetic/RREF/producer imports',
            'semilinear_candidates_tested':6*728,'non_similitudes_rejected':bad,'similitude_group_order':336,'lifts_per_j_epsilon':28,
            'normal_root_subgroup_order':14,'faithful_quotient_permutation_order':24,'catalog_subspaces_checked':1120,
            'orbits':orbits,'R_orbit':r[0],'S_orbit_sizes':list(map(len,ss)),'S_representatives':[o[0] for o in ss],
            'type_stabilizer':stabilizer,'inverse_embeddings':inverses,'inverse_residue_frobenius_exponents':exponents,
            'reflex_residue_exponent':286,'reflex_residue_image_order':28,'reflex_image_mod_roots_order':2,
            'exact_partition_matches_candidate':True,'geometric_galois_bridge_proved_by_this_loop':False,'pins':packet['pins']}

if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
    start=time.monotonic();data=run();data['wall_seconds']=time.monotonic()-start
    print(json.dumps(data,indent=2,sort_keys=True))
