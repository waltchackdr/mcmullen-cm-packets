"""Independent V140 check on L^2 homology, using only the standard library.

No primary arithmetic imported. Constructs the actual transpose Z_A action
on the common principal basis; constructs sigma29 using a cyclic vector;
proves saturation by integer Euclidean row operations on 144 constraints.
Checks the trace Gram physically, then enumerates with reversed coordinates.
The geometric CM identification remains an adopted human-proof input.
"""
import hashlib
import json
import math
import resource
import signal
import sys
import time
from fractions import Fraction as F
from pathlib import Path

resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
signal.alarm(30)
if not __debug__:
    raise RuntimeError('Assertions required')
start = time.monotonic()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
record_path = ROOT / 'certificates/a_unmarked_polarization_data.json'
record_bytes = record_path.read_bytes()
record = json.loads(record_bytes)

# Arithmetic in Q[t]/(t^6+t^5-6t^4-6t^3+8t^2+8t+1).
mod = [1, 8, 8, -6, -6, 1]
def elt(v): return tuple(F(x) for x in list(v)+[0]*(6-len(v)))
zero, one, t = elt([]), elt([1]), elt([0, 1])
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def scale(a,c): return tuple(x*c for x in a)
def sub(a,b): return add(a,scale(b,-1))
def mul(a,b):
    v=[F(0)]*11
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if x and y: v[i+j]+=x*y
    for k in range(10,5,-1):
        for j,c in enumerate(mod): v[k-6+j]-=v[k]*c
    return tuple(v[:6])
def eye(n): return [[F(i==j) for j in range(n)] for i in range(n)]
def transpose(A): return [list(x) for x in zip(*A)]
def mm(A,B):
    out=[[F(0)]*len(B[0]) for _ in A]
    for i,row in enumerate(A):
        for k,x in enumerate(row):
            if x:
                for j,y in enumerate(B[k]):
                    if y: out[i][j]+=x*y
    return out
def inverse(A):
    n=len(A);R=[[F(x) for x in a]+b for a,b in zip(A,eye(n))]
    for j in range(n):
        p=next(i for i in range(j,n) if R[i][j]);R[p],R[j]=R[j],R[p]
        x=R[j][j];R[j]=[v/x for v in R[j]]
        for i in range(n):
            if i!=j and R[i][j]:
                x=R[i][j];R[i]=[v-x*w for v,w in zip(R[i],R[j])]
    return [r[n:] for r in R]
def reg(a): return transpose([mul(a,elt([0]*j+[1])) for j in range(6)])
def inv(a): return tuple(row[0] for row in inverse(reg(a)))
tp=[one]
for j in range(10): tp.append(mul(tp[-1],t))
def poly(c):
    ans=zero
    for j,x in enumerate(c): ans=add(ans,scale(tp[j],x))
    return ans
def div(a,b):return mul(a,inv(b))
def blocks(A):
    rs=[[reg(x) for x in row] for row in A]
    return [rs[i][0][j]+rs[i][1][j] for i in range(2) for j in range(6)]
def integer(A):
    assert all(F(x).denominator==1 for row in A for x in row)
    return [[int(x) for x in row] for row in A]
def trace(A):return sum(A[i][i] for i in range(len(A)))

u=sub(t,one);beta=poly([-2,-8,-2,6,1,-1]);a=inv(mul(beta,mul(u,u)));b=inv(mul(mul(u,u),u));delta=mul(a,b)
T=blocks([[a,scale(b,3)],[zero,b]])
Ti=inverse(T)
assert mm(T,Ti)==eye(12)
At=[[poly([0,3,-1,-4,0,1]),mul(t,poly([0,-10,-3,6,1,-1]))],
    [div(poly([-1,0,0,1]),t),poly([0,-6,1,5,0,-1])]]
ZA=[[div(add(At[i][j],t if i==j else zero),sub(tp[2],one)) for j in range(2)] for i in range(2)]
assert add(ZA[0][0],ZA[1][1])==t
assert sub(mul(ZA[0][0],ZA[1][1]),mul(ZA[0][1],ZA[1][0]))==one
ZH=blocks(transpose(ZA))
Zpowers=[eye(12)]
for j in range(42): Zpowers.append(mm(Zpowers[-1],ZH))
assert Zpowers[42]==eye(12)
# The cyclic vector (0,1) maps to 1 in the primary field realization,
# but here all calculations stay on the original rational L^2 homology.
cyclic=transpose([[row[6] for row in Zpowers[j]] for j in range(12)])
conjugate_cyclic=transpose([[row[6] for row in Zpowers[(29*j)%42]] for j in range(12)])
sigma=mm(conjugate_cyclic,inverse(cyclic))
assert mm(sigma,sigma)==eye(12)
assert mm(sigma,ZH)==mm(Zpowers[29],sigma)
ops=[mm(mm(Ti,Zpowers[j]),T) for j in range(12)]
ops += [mm(mm(mm(Ti,Zpowers[j]),sigma),T) for j in range(12)]
phi={1,5,17,19,29,31}
stabilizer=[g for g in range(42) if math.gcd(g,42)==1 and {g*i%42 for i in phi}==phi]
assert stabilizer==[1,29]

# Reconstruct E from the L-trace alternating form in the common T basis.
E=[[F(0)]*12 for _ in range(12)]
for i in range(6):
    for j in range(6):
        v=trace(reg(mul(delta,tp[i+j])))
        E[i][j+6]=-v;E[j+6][i]=v
E=integer(E)
data58_bytes=(ROOT/'certificates/principal_unmarked_bc_data.json').read_bytes()
data58_hash=hashlib.sha256(data58_bytes).hexdigest()
assert data58_hash=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
data58=json.loads(data58_bytes)
assert E==[[F(x) for x in row] for row in data58['data']['principal_alternating_matrix']]
Ei=integer(inverse(E))
q=sub(elt([4]),tp[2]);Qq=integer(blocks([[q,zero],[zero,q]]));Qi=integer(inverse(Qq))
assert mm(transpose(Qq),E)==mm(E,Qq)

# Treat the recorded coefficient columns only as candidates. Prove they
# constitute the entire integral lattice by a distinct saturation check.
L=[[F(x) for x in row] for row in record['coefficient_basis_a_then_b']]
U=[[F(x) for x in row] for row in record['U']]
coeff=mm(L,U)
maps=[]
for j in range(24):
    matrix=[[sum(coeff[k][j]*ops[k][i][h] for k in range(24)) for h in range(12)] for i in range(12)]
    maps.append(integer(matrix))
B=[[maps[j][i][h] for j in range(24)] for i in range(12) for h in range(12)]
R=[row[:] for row in B]
pivots=[]
for j in range(24):
    p=next(i for i in range(j,len(R)) if R[i][j]);R[p],R[j]=R[j],R[p]
    for i in range(j+1,len(R)):
        while R[i][j]:
            quotient=R[j][j]//R[i][j]
            R[j]=[x-quotient*y for x,y in zip(R[j],R[i])]
            R[j],R[i]=R[i],R[j]
    if R[j][j]<0:R[j]=[-x for x in R[j]]
    pivots.append(R[j][j])
assert all(not any(row) for row in R[24:])
assert math.prod(pivots)==1
# Integer row lattice = Z^24 implies Bv integral iff v integral.
# The rank24 physical maps span the full human-identified Hodge algebra.

weighted=[mm(mm(mm(Qi,Ei),transpose(V)),E) for V in maps]
H=[[sum(weighted[i][k][l]*maps[j][l][k] for k in range(12) for l in range(12)) for j in range(24)] for i in range(24)]
assert H==transpose(H)
assert H==[[F(x) for x in row] for row in record['reduced_Gram']]

# Recompute LDL after reversing all coordinates. Thus the branch order is
# opposite to the primary checker; no primary LDL or bounds are trusted.
G=[list(reversed(row)) for row in reversed(H)]
n=24;ell=[[F(0)]*n for _ in range(n)];ds=[]
for i in range(n):
    ell[i][i]=F(1)
    for j in range(i):ell[i][j]=(G[i][j]-sum(ell[i][k]*ell[j][k]*ds[k] for k in range(j)))/ds[j]
    ds.append(G[i][i]-sum(ell[i][k]**2*ds[k] for k in range(i)))
assert all(x>0 for x in ds)
Gi=inverse(G)
bounds=[math.isqrt((12*Gi[i][i]).numerator//(12*Gi[i][i]).denominator) for i in range(n)]
stack=[(n-1,F(12),[0]*n)];nodes=0;leaves=[]
while stack:
    i,remaining,x=stack.pop();nodes+=1
    if nodes>20000:raise RuntimeError('Independent enumeration node cap exceeded')
    if i<0:
        value=sum(x[j]*G[j][k]*x[k] for j in range(n) for k in range(n))
        assert value==12-remaining
        leaves.append((x,value));continue
    center=sum(ell[j][i]*x[j] for j in range(i+1,n))
    for v in range(-bounds[i],bounds[i]+1):
        term=ds[i]*(v+center)**2
        if term<=remaining:
            child=x[:];child[i]=v;stack.append((i-1,remaining-term,child))
assert leaves==[([0]*24,F(0))]
print(json.dumps({'status':'INDEPENDENT_FULL_A_OBSTRUCTION_PASS','method':'rational L^2 homology; cyclic-vector sigma29; Euclidean row-saturation; physical Rosati Gram; reversed-coordinate LDL enumeration','dimension':24,'stabilizer':stabilizer,'constraint_row_lattice_index':math.prod(pivots),'all_pivots':pivots,'nodes':nodes,'nonzero_closed_trace12_ball_vectors':0,'elapsed_seconds':time.monotonic()-start,'python_version':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'primary_data_sha256':hashlib.sha256(record_bytes).hexdigest(),'adopted_principal_data_sha256':data58_hash,'scope':'Independent finite full-lattice obstruction; actual geometric model and Hodge type remain adopted human inputs.'},indent=2,sort_keys=True))
