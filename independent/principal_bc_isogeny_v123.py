"""V123 independent companion field and Schur-prefix enumeration certificate.
Reads frozen C58 and V123 data only; imports no repository algorithms.
Actual geometric lattice and type identifications are inherited proof inputs.
"""
import argparse, hashlib, json, math, resource, signal, time
from pathlib import Path
import sympy as s

if not __debug__: raise RuntimeError("Enabled assertions required")
resource.setrlimit(resource.RLIMIT_CPU, (20,20))
signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError()))
signal.alarm(30)
start=time.monotonic()
parser=argparse.ArgumentParser()
parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'certificates/principal_bc_isogeny_data.json')
parser.add_argument('--c58-data',type=Path,default=Path(__file__).resolve().parents[1]/'certificates/principal_unmarked_bc_data.json')
args=parser.parse_args()
new_raw=args.data.read_bytes()
assert hashlib.sha256(new_raw).hexdigest()=='4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee'
new=json.loads(new_raw)
assert new['schema']=='PRINCIPAL_BC_MULTIPLIER2_V123_1'
assert new['C58_data_sha256']=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
path=args.c58_data
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
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
# Complete breadth expansion; real-prefix Schur minima prune whole subtrees.
bounds=[0]*12
for i in range(12):
 while (bounds[i]+1)**2<=24*Hi[i,i]:bounds[i]+=1
assert bounds==[1,1,1,1,2,1,1,1,1,1,1,1]
frontier={()};counts=[1]
for k in range(1,13):
 prefix=12-k
 schur=H[prefix:,prefix:] if prefix==0 else H[prefix:,prefix:]-H[prefix:,:prefix]*H[:prefix,:prefix].inv()*H[:prefix,prefix:]
 denominator=s.ilcm(1,*[v.q for v in schur])
 mat=[[int(v*denominator) for v in row] for row in schur.tolist()]
 successors=set()
 for tail in frontier:
  for value in range(-bounds[prefix],bounds[prefix]+1):
   v=(value,)+tail
   q=sum(v[i]*mat[i][j]*v[j] for i in range(k) for j in range(k))
   if q<=24*denominator:successors.add(v)
 frontier=successors;counts.append(len(frontier))
assert counts==[1,3,7,17,33,49,79,103,125,133,145,163,175]
ball={v:int((s.Matrix(v).T*H*s.Matrix(v))[0]) for v in frontier}
shell={v for v,q in ball.items() if q==24}
ci=ideals['C'].inv();sb=twist*ideals['B'];E=forms['B'];assert E==forms['C']
accepted=[];norm_rejects=integral_rejects=0
for v in sorted(shell):
 old=U*s.Matrix(v);m=J*old
 if M(m)*conj*m!=2*r:norm_rejects+=1;continue
 A=ci*M(m)*sb
 if any(x.q!=1 for x in A):integral_rejects+=1;continue
 assert A.T*E*A==2*E and A.det()==64
 snf=smith_normal_form(A,domain=ZZ)
 assert [abs(int(snf[i,i])) for i in range(12)]==[1]*6+[2]*6
 adjoint=E.inv()*A.T*E
 assert all(x.q==1 for x in adjoint) and adjoint*A==2*s.eye(12)
 accepted.append((m,A))
assert len(shell)==84 and norm_rejects==72 and integral_rejects==6 and len(accepted)==6
assert new['expected']=={'automorphism_order':6,'degree':64,'integrality_rejected':6,'isogenies':6,'multiplier':2,'nodes':1033,'norm_rejected':72,'trace_shell':84}
assert sum(counts)==new['expected']['nodes']
expected={tuple(V(v)) for v in new['expected_scalar_set']}
assert len(new['expected_scalar_set'])==len(expected)==6
assert {tuple(m) for m,A in accepted}==expected
m0=V(new['witness_m']);A0=s.Matrix(new['witness_matrix'])
assert any(m==m0 and A==A0 for m,A in accepted)
groups={}
for name in ('B','C'):
 G=ideals[name].inv()*C**7*ideals[name]
 assert all(x.q==1 for x in G) and G.det()==1 and G.T*E*G==E
 assert G**6==s.eye(12) and all(G**j!=s.eye(12) for j in range(1,6))
 groups[name]=G
assert {tuple(C**(7*j)*m0) for j in range(6)}==expected
for j in range(6):
 mj=C**(7*j)*m0;Aj=groups['C']**j*A0
 assert any(m==mj and A==Aj for m,A in accepted)
for m,A in accepted:assert A*groups['B']==groups['C']*A
signal.alarm(0)
print(json.dumps({'status':'INDEPENDENT_PRINCIPAL_BC_MULTIPLIER2_COMPLETE_PASS',
 'arithmetic':'Rational cyclotomic companion matrices; exact real-prefix Schur enumeration; no repository algorithm imports',
 'data_sha256':hashlib.sha256(new_raw).hexdigest(),'C58_data_sha256':hashlib.sha256(raw).hexdigest(),
 'closed_ball_vectors':len(ball),'trace24_vectors':len(shell),'full_norm_rejected':norm_rejects,
 'integrality_rejected':integral_rejects,'isogenies':len(accepted),'multiplier':2,'degree_each':64,
 'kernel_each':'(Z/2Z)^6','target_mu6_orbit':True,'source_target_mu6_integral_polarized':True,
 'scope':'Specified principal B/C; minimal positive scalar multiplier uses C58 exclusion of1; no minimum arbitrary unpolarized degree.'},sort_keys=True))
