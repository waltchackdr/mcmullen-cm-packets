"""Independent exact companion-matrix audit; imports no producer arithmetic."""
import hashlib, json, math, resource, signal, time
from pathlib import Path
import sympy as s

if not __debug__: raise RuntimeError("Enabled assertions required")
resource.setrlimit(resource.RLIMIT_CPU, (20,20))
signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError()))
signal.alarm(30)
start=time.monotonic()
path=Path(__file__).resolve().parents[1]/'certificates/principal_unmarked_bc_data.json'
raw=path.read_bytes()
assert hashlib.sha256(raw).hexdigest()=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
record=json.loads(raw); d=record['data']
assert record['schema']=='PRINCIPAL_UNMARKED_BC_FINITE_NORM_V122_1'
assert record['actual_tangent_types']=={'B':[1,5,23,25,29,31],'C':[1,5,11,13,17,19]}
assert record['galois_exponent_B_to_C']==31
x=s.Symbol('x'); f=s.Poly(s.cyclotomic_poly(42,x),x)
C=s.zeros(12)
for i in range(11): C[i+1,i]=1
for i in range(12): C[i,11]=-f.nth(i)
I=s.eye(12); e=I[:,0]; powers=[C**i for i in range(12)]
def V(v): return s.Matrix([s.Rational(a) for a in v])
def M(v): return sum((v[i]*powers[i] for i in range(12)),s.zeros(12))
def inv(v): return M(v).inv()*e
def mul(v,w): return M(v)*w
def sigma(k): return s.Matrix.hstack(*[(C**(k*j))*e for j in range(12)])
conj=sigma(41); twist=sigma(31)
t=(C+C**41)*e; T=M(t); u=t-e
beta=(-T**5+T**4+6*T**3-2*T**2-8*T-2*I)*e
a=inv(M(beta)*M(u)**2*e); b=inv(M(u)**3*e)
cB=inv(t); dB=t
cC=(-I+5*T+4*T**2-T**3-T**4)*inv(t)
dC=(-T+4*T**2+4*T**3-T**4-T**5)*e
delta=mul(a,b)
tr=s.Matrix([[p.trace() for p in powers]])
forms={}; ideals={}; xis={}
for name,c,dd in [('B',cB,dB),('C',cC,dC)]:
    w=mul(C*e-dd,inv(c))
    q=mul(a,w); p=mul(b,e+3*w)
    B=s.Matrix.hstack(*[T**j*q for j in range(6)],*[T**j*p for j in range(6)])
    assert B==s.Matrix.hstack(*map(V,d['I_'+name+'_basis']))
    xi=-mul(c,inv((C-C**41)*e)); xis[name]=xi; ideals[name]=B
    assert xi==V(d['xi_'+name])
    gram=s.Matrix(12,12,lambda i,j:(tr*M(xi)*M(B[:,i])*conj*B[:,j])[0])
    assert gram==s.Matrix(d['principal_alternating_matrix']).applyfunc(s.Rational)
    expected=s.Matrix(12,12,lambda i,j:0 if (i<6)==(j<6) else (-1 if i<6 else 1)*(tr*T**((i%6)+(j%6))*delta)[0]/2)
    assert gram==expected and gram.det()==1
    forms[name]=gram
units=[k for k in range(42) if math.gcd(k,42)==1]
Btype={1,5,23,25,29,31}; Ctype={1,5,11,13,17,19}
assert [g for g in units if {g*k%42 for k in Ctype}==Btype]==[31]
assert all([g for g in units if {g*k%42 for k in typ}==typ]==[1] for typ in (Btype,Ctype))
r=mul(twist*xis['B'],inv(xis['C']))
assert r==V(d['relative_norm_target_r']) and conj*r==r
J=M(inv(twist*ideals['B'][:,0]))*ideals['C']
assert J==s.Matrix.hstack(*map(V,d['search_overlattice_J_basis']))
weighted_trace=tr*M(inv(r))
Q=s.Matrix(12,12,lambda i,j:(weighted_trace*M(J[:,i])*conj*J[:,j])[0])
assert Q==s.Matrix(d['positive_definite_trace_Gram_Q']).applyfunc(s.Rational) and Q==Q.T
lattice=record['lattice']; U=s.Matrix(lattice['U_columns_old_coordinates'])
assert all(v.q==1 for v in U) and abs(U.det())==1
H=U.T*Q*U
assert H==s.Matrix(lattice['reduced_Gram_H']).applyfunc(s.Rational)
assert all(H[:k,:k].det()>0 for k in range(1,13))
assert all(v.q==1 for v in H) and all(H[i,i]%2==0 for i in range(12))
Hi=H.inv()
assert all(1<=12*Hi[i,i]<4 for i in range(12))
assert lattice['reduced_coordinate_bounds']==[1]*12
# Full Cartesian box evaluated via six-coordinate split quadratic forms.
# No producer code, LDL decomposition or recursive pruning is used.
import itertools
h=[[int(v) for v in row] for row in H.tolist()]
vectors=list(itertools.product((-1,0,1),repeat=6))
left=[sum(v[i]*h[i][j]*v[j] for i in range(6) for j in range(6)) for v in vectors]
right=[sum(v[i]*h[i+6][j+6]*v[j] for i in range(6) for j in range(6)) for v in vectors]
visited=0;ball=[]
for a,qa in zip(vectors,left):
    cross=[2*sum(a[i]*h[i][j+6] for i in range(6)) for j in range(6)]
    for b,qb in zip(vectors,right):
        value=qa+qb+sum(c*x for c,x in zip(cross,b));visited+=1
        if value<=12:ball.append((a+b,value))
assert visited==3**12 and ball==[((0,)*12,0)]
assert H[0,0]==14
print(json.dumps({'status':'INDEPENDENT_PRINCIPAL_UNMARKED_BC_NORM_OBSTRUCTION_PASS',
 'data_sha256':hashlib.sha256(raw).hexdigest(),'integer_box_points':visited,
 'nonzero_vectors_of_norm_at_most12':0,'global_lattice_minimum':14,
 'arithmetic':'Rational cyclotomic companion matrices and complete split-box enumeration; no producer imports'},sort_keys=True))
