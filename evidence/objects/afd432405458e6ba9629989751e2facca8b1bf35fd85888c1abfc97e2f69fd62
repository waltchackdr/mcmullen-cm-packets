#!/usr/bin/env python3
"""Portable independent cyclotomic verifier for the full 42-trace packet.

Reconstructs all 3,520 matrices, 368 points and 736 root bridges with the
separate degree-12 arithmetic path. No quadratic-path implementation imports,
experiment/build dependencies, geometry searches or output-file writes.
Printed PASS concerns this arithmetic check; promotion is ledger-governed.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from pathlib import Path
import json
import runpy
import sys

if not __debug__:
    raise RuntimeError("Assertions must be enabled")
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/"certificates/full_packet_v62"
ARITHMETIC = ROOT/"certificates/d3_completeness_v60.py"
SOURCE_SHA256 = {
    "certificates/d3_completeness_v60.py": "7341c4dddc68e06e3288cda31372e6ad3a7617074a01cb2292526eba17baaf17",
    # These three files are historical discovery provenance, not imported code.
    "independent/d3_completeness_quadratic_v60.py": "7b34e8c676a0045b8c02d1ae065b395003bb48191d1cddf8d085a458db4f0363",
    "independent/packet_order_coverage_stdlib_v58.py": "ed328d11d87d2cc9e93322bddb2b3207b1edf8c5d8dc66dea2026f0a07ea07c3",
    "independent/quadrilateral_k17_stdlib_v41.py": "c5bd45f334a028820a071b5c758acfb8e48afc1e8454b18a3174dd04ca58a38c",
}
ARTIFACT_SHA256 = {
    "all_42_matrix_classes.jsonl": "08492a608ce13b99b6e9470241379a1a01746e78539b126aeb785a15af9b9761",
    "ambient_points.json": "c0e84da6c266c087cbce65a6e8aca8e3f30350ce56b5a1089c4ae4978a381f36",
    "local_quotients.json": "9ef01be9df80366e398b75a7aad39b4705ea06c9a84ff884d6b4a08b69f05613",
    "module_templates.json": "09ccff217511a67e946567135f9ffb74da020843621d0895ecce8746d25c245d",
    "root_bridges.jsonl": "13dcf16cd44e5359659b307dc9f0eb45d8d0864b8dda9c3aabb4291bb620ece6",
}
HISTORICAL_SUMMARY_SHA256 = "839341e7b52f3d8726fe9f64d16685327315aa4bb24b3e4f3fc5d3ae3f2a956c"
DISCOVERY_IMPLEMENTATION_SHA256 = "d1877461017b7de42068f718ab1ecbb05fa8b10b0445a0589c8cb602d75eb3d5"
ENTRY_SHA256 = sha256(Path(__file__).read_bytes()).hexdigest()


def check_inputs():
    for relative, expected in SOURCE_SHA256.items():
        assert sha256((ROOT/relative).read_bytes()).hexdigest() == expected, relative
    for filename, expected in ARTIFACT_SHA256.items():
        assert sha256((BASE/filename).read_bytes()).hexdigest() == expected, filename
    assert sha256((BASE/"summary.json").read_bytes()).hexdigest() == HISTORICAL_SUMMARY_SHA256
    assert sha256(Path(__file__).read_bytes()).hexdigest() == ENTRY_SHA256, "Entry changed during execution"


check_inputs()
arithmetic = runpy.run_path(str(ARITHMETIC))
assert Path(arithmetic["__file__"]).resolve() == ARITHMETIC.resolve()
L, M, relative = arithmetic["L"], arithmetic["M"], arithmetic["relative_coordinates"]
L.divide_integrally = cache(L.divide_integrally)
O, Z, T = L.one, L.zero, L.x
U = L.reduce((4, 0, -1))
PI = L.reduce((5, -2, -4, 1, 1))
neg = lambda a: L.scale(a, -1)
sub = lambda a, b: L.add(a, neg(b))


def mm(a, b):
    return (L.add(L.mul(a[0], b[0]), L.mul(a[1], b[2])),
            L.add(L.mul(a[0], b[1]), L.mul(a[1], b[3])),
            L.add(L.mul(a[2], b[0]), L.mul(a[3], b[2])),
            L.add(L.mul(a[2], b[1]), L.mul(a[3], b[3])))


def det(a):
    return sub(L.mul(a[0], a[3]), L.mul(a[1], a[2]))


def barred(a):
    return a[3], neg(a[1]), neg(a[2]), a[0]


def matrix(raw):
    return tuple(tuple(c) for c in raw)


@cache
def reciprocal(a):
    columns = [L.mul(a, L.reduce((0,)*j+(1,))) for j in range(6)]
    inverse = arithmetic["inverse_matrix"]([[columns[j][i] for j in range(6)] for i in range(6)])
    value = tuple(inverse[i][0] for i in range(6))
    assert L.mul(value, a) == O
    return value


def eval_mod(a, r, p):
    return sum(c*pow(r,j,p) for j,c in enumerate(a)) % p


def order_mod(a, p):
    value = 1
    for n in range(1, p):
        value = value*a % p
        if value == 1:
            return n
    raise AssertionError("Not a nonzero residue")


# At 127 use the split evaluation-ratio model, rather than projective pairs.
p, r, z1, z2 = 127, 94, 40, 54
assert eval_mod(L.modulus, r, p) == 0
assert [(a*a-r*a+1) % p for a in (z1,z2)] == [0,0]
ratio = lambda a,b: (a+b*z1)*pow((a+b*z2) % p,-1,p) % p
eps, zz, beta = ratio(2,-r), ratio(0,1), ratio(1,9)
assert (order_mod(eps,p),order_mod(zz,p),order_mod(beta,p)) == (42,21,126)
assert pow(eps,22,p) == zz
image = {pow(eps,j,p) for j in range(42)}
ratios = [ratio(1,s) for s in (0,9,100)]
cosets = [{a*v % p for v in image} for a in ratios]
assert len(set.union(*cosets)) == sum(map(len,cosets)) == 126
assert ratios[2] == pow(ratios[1],-1,p)
assert (-(9)*pow(1+9*r,-1,p)) % p == 100

# At 41 use u/bar(u) in the norm-one subgroup of the quadratic finite field.
p41, r41 = 41, 22
assert not [a for a in range(p41) if (a*a-r41*a+1) % p41 == 0]
assert eval_mod(L.modulus,r41,p41) == 0


def mul41(a,b):
    return ((a[0]*b[0]-a[1]*b[1]) % p41,
            (a[0]*b[1]+a[1]*b[0]+r41*a[1]*b[1]) % p41)


def conj41(a):
    return ((a[0]+r41*a[1]) % p41, -a[1] % p41)


def inv41(a):
    n=mul41(a,conj41(a)); assert n[1] == 0 and n[0]
    ni=pow(n[0],-1,p41)
    return tuple(v*ni % p41 for v in conj41(a))


def normone41(a):
    return mul41(a,inv41(conj41(a)))


def powers41(a):
    result, current = [], (1,0)
    while current not in result:
        result.append(current); current=mul41(current,a)
    assert current == (1,0)
    return result


quotient41={normone41((a,b)) for a in range(41) for b in range(41) if (a,b)!=(0,0)}
epsilon41=normone41((2,-22))
z41=normone41((0,1))
epsilon_powers41=powers41(epsilon41)
assert len(quotient41) == len(epsilon_powers41) == 42
assert set(epsilon_powers41) == quotient41
assert epsilon_powers41[22] == z41 and len(powers41(z41)) == 21

z=M.x
alpha0=M.add(M.add(M.one,z),M.pow(z,4))
alpha=[relative(M.mul(M.pow(z,k),alpha0)) for k in range(42)]
conductors=[]
for a,b in alpha:
    assert L.add(L.add(L.mul(a,a),L.mul(T,L.mul(a,b))),L.mul(b,b)) == PI
    norm=abs(L.norm(b))
    if norm == 1:
        key=(1,0,0)
    else:
        prime, exponent=(7,2) if norm==49 else (norm,1)
        assert prime in (7,41,43,83,127)
        roots=[v for v in range(prime) if eval_mod(L.modulus,v,prime)==0 and eval_mod(b,v,prime)==0]
        assert len(roots)==1
        key=(prime,roots[0],exponent)
        if prime==7:
            quotient=L.divide_integrally(b,L.pow(L.reduce((-1,1)),exponent))
            assert abs(L.norm(quotient))==1
    conductors.append(key)
assert len(set(conductors))==9
strata={"OM":((1,0,0),O,1,(0,)),
        "R7":((7,1,1),alpha[0][1],0,(0,)),
        "R7sq":((7,1,2),alpha[12][1],12,(0,)),
        "R41_12":((41,12,1),alpha[14][1],14,(0,)),
        "R41_22":((41,22,1),alpha[5][1],5,(0,)),
        "R43_20":((43,20,1),alpha[4][1],4,(0,)),
        "R43_31":((43,31,1),alpha[8][1],8,(0,)),
        "R83":((83,36,1),alpha[2][1],2,(0,3)),
        "R127":((127,94,1),alpha[10][1],10,(0,9,100))}
generators=tuple(map(L.reduce,((1,0,-3,0,1),(-1,-2,1,1),(1,3,0,-4,0,1),(0,-2,0,1),(-1,),(4,0,-1))))
units=[]
for mask in range(64):
    value=O
    for i,g in enumerate(generators):
        if mask & (1<<i): value=L.mul(value,g)
    units.append(value)
unit_inverses=[L.divide_integrally(O,value) for value in units]


def module_basis(name,j):
    _,g,_,slopes=strata[name]
    return O,Z,L.scale(O,slopes[j]),g


@cache
def action(k,name,j,mask):
    _,g,_,slopes=strata[name]; s=slopes[j]
    a,b=alpha[k]
    multiplication=(a,neg(b),b,L.add(a,L.mul(T,b)))
    basis=module_basis(name,j)
    n=L.add(L.scale(O,1+s*s),L.scale(T,s))
    result=(sub(a,L.scale(b,s)),neg(L.mul(L.mul(units[mask],g),b)),
            L.mul(L.mul(L.divide_integrally(b,g),n),unit_inverses[mask]),
            L.add(a,L.mul(L.add(T,L.scale(O,s)),b)))
    # Undo the unit twist; directly verify the full multiplication basis equation.
    plain=(result[0],L.mul(result[1],unit_inverses[mask]),
           L.mul(result[2],units[mask]),result[3])
    assert mm(basis,plain)==mm(multiplication,basis)
    assert det(result)==PI
    return result


def partner(name,j,mask):
    if name=="R127" and j:
        return 3-j,mask^16
    return j,mask^(48 if name=="R83" and j else 16)


def canonical(name,j,mask):
    return min((j,mask),partner(name,j,mask))


summary=json.loads((BASE/"summary.json").read_text())
for filename,expected in summary["artifact_sha256"].items():
    assert sha256((BASE/filename).read_bytes()).hexdigest()==expected
assert summary["implementation_sha256"] == DISCOVERY_IMPLEMENTATION_SHA256
assert summary["artifact_sha256"] == ARTIFACT_SHA256
for filename,expected in summary["source_sha256"].items():
    assert sha256((ROOT/filename).read_bytes()).hexdigest()==expected
local=json.loads((BASE/"local_quotients.json").read_text())
assert (local['R41_22']['quotient_order'],local['R41_22']['epsilon_order'],local['R41_22']['picard_order'])==(42,42,1)
assert (local['R127']['quotient_order'],local['R127']['epsilon_order'],local['R127']['picard_order'])==(126,42,3)

templates=json.loads((BASE/"module_templates.json").read_text())
assert len(templates)==12
assert {(row['stratum'],row['picard_class']) for row in templates} == {
    (name,j) for name,(_,_,_,slopes) in strata.items() for j in range(len(slopes))}
epsilon=M.add(M.one,M.scale(M.pow(z,2),-1))
eta83=relative(M.mul(epsilon,M.pow(z,28)))
for item in templates:
    name,j=item['stratum'],item['picard_class']
    key,g,k,slopes=strata[name]
    basis=module_basis(name,j)
    assert matrix(item['basis_columns_in_one_z'])==basis
    assert det(basis)==g and abs(L.norm(g))==item['index_norm']
    assert tuple(item['conductor'])==key
    assert matrix(item['norm_pi_action'])==action(k,name,j,0)
    cj=partner(name,j,0)[0]
    assert item['conjugate_picard_class']==cj
    c=matrix(item['anti_intertwiner'])
    expected_det=neg(U if name=='R83' and j else O)
    assert det(c)==expected_det==tuple(item['anti_determinant'])
    eta=eta83 if name=='R83' and j else (O,Z)
    a,b=eta
    multeta=(a,neg(b),b,L.add(a,L.mul(T,b)))
    assert mm(module_basis(name,cj),c)==mm(mm(multeta,(O,T,Z,neg(O))),basis)
    assert mm(c,action(k,name,j,0))==mm(barred(action(k,name,cj,0)),c)
    inverse_det=L.divide_integrally(O,det(c))
    cinverse=tuple(L.mul(inverse_det,x) for x in (c[3],neg(c[1]),neg(c[2]),c[0]))
    assert mm(c,cinverse)==(O,Z,Z,O)

bridges={}
for line in (BASE/'root_bridges.jsonl').read_text().splitlines():
    row=json.loads(line);name,j,m=row['stratum'],row['picard_class'],row['orientation_mask']
    cj,cm=partner(name,j,m)
    assert (row['target_picard_class'],row['target_orientation_mask'])==(cj,cm)
    key=(name,j,m); assert key not in bridges
    G=matrix(row['SL2_root_bridge']); assert det(G)==O
    bridges[key]=G
assert len(bridges)==736
seen,incidences,per_trace=set(),Counter(),Counter()
for line in (BASE/'all_42_matrix_classes.jsonl').read_text().splitlines():
    row=json.loads(line);k,name,j,m=row['trace_index'],row['stratum'],row['picard_class'],row['orientation_mask']
    key=(k,name,j,m); assert key not in seen;seen.add(key)
    actual=action(k,name,j,m);assert matrix(row['matrix_row_major'])==actual
    cj,cm=partner(name,j,m);assert row['root_partner']==[cj,cm]
    G=bridges[name,j,m]
    assert mm(G,actual)==mm(barred(action(k,name,cj,cm)),G)
    jj,mmask=canonical(name,j,m)
    label=f'{name}:pic{jj}:{mmask:02d}'
    assert row['canonical_point_label']==label
    incidences[label]+=1;per_trace[k]+=1
expected_domain=set()
for k,(p,r,e) in enumerate(conductors):
    for name,(key,g,ref,slopes) in strata.items():
        if name!='OM' and not(key[:2]==(p,r) and key[2]<=e):continue
        for j in range(len(slopes)):
            for m in range(32 if name=='OM' else 64): expected_domain.add((k,name,j,m))
assert seen==expected_domain and len(seen)==3520

points=json.loads((BASE/'ambient_points.json').read_text())
point_keys=set();tuple_values=set();counts=Counter()
for row in points:
    name,j,m=row['stratum'],row['picard_class'],row['orientation_mask']
    assert row['point_label']==f'{name}:pic{j}:{m:02d}'
    assert canonical(name,j,m)==(j,m)
    key=(name,j,m);assert key not in point_keys;point_keys.add(key)
    _,g,k,slopes=strata[name];s=slopes[j]
    q=L.mul(g,units[m]);n=L.add(L.scale(O,1+s*s),L.scale(T,s));ni=reciprocal(n)
    x2=neg(L.mul(L.mul(q,L.add(T,L.scale(O,2*s))),ni))
    h4=L.mul(L.mul(L.mul(q,q),U),L.mul(ni,ni))
    mod2=L.mul(L.mul(q,q),ni)
    for field,value in [('twice_real_part',x2),('four_times_height_squared',h4),('modulus_squared',mod2)]:
        assert tuple(map(Q,row[field]))==value
    assert L.add(L.mul(x2,x2),h4)==L.scale(mod2,4)
    assert matrix(row['norm_pi_matrix'])==action(k,name,j,m)
    assert row['root_partner']==list(partner(name,j,m))
    assert (x2,h4) not in tuple_values;tuple_values.add((x2,h4))
    counts[name]+=1
assert len(points)==len(point_keys)==len(incidences)==368
assert dict(counts)==summary['point_counts_by_multiplier']
assert all(incidences[row['point_label']]==summary['incidences_per_point'][row['stratum']] for row in points)
assert counts['R127']==96
assert sum(row['stratum']=='R127' and row['picard_class']==1 for row in points)==64
assert not any(row['stratum']=='R127' and row['picard_class']==2 for row in points)

output={"status":"PASS_INDEPENDENT_FULL_PACKET_CYCLOTOMIC_V62",
        "geometry_search_run":False,"imported_proposed_v62_or_quadratic_code":False,
        "matrix_classes":3520,"canonical_ambient_points":368,"exact_root_bridges":736,
        "literal_canonical_tuples_pairwise_distinct":True,
        "local_models":{"R41_22":"u/bar(u) norm-one subgroup in F41^2",
                        "R127":"evaluation ratio at z=40 and54 in F127*"},
        "R127_epsilon_ratio":eps,"R127_z_ratio":zz,"R127_beta_ratio":beta,
        "R127_picard_cosets":[sorted(v) for v in cosets],
        "point_counts":dict(counts),"matrix_counts_per_trace":dict(per_trace),
        "historical_discovery_implementation_sha256":summary['implementation_sha256'],
        "reviewed_artifact_sha256":summary['artifact_sha256'],
        "arithmetic_source_sha256":sha256(ARITHMETIC.read_bytes()).hexdigest(),
        "entrypoint_sha256":ENTRY_SHA256,
        "source_sha256":SOURCE_SHA256,
        "historical_summary_sha256":HISTORICAL_SUMMARY_SHA256}
check_inputs()
print(json.dumps({k:v for k,v in output.items() if k!='R127_picard_cosets'},indent=2,sort_keys=True))
