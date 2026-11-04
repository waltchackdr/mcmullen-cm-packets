"""Independent C60 companion/mod3/membership certificate.
No primary algorithms, coefficient arithmetic or HNF imported.
Supplied order/saturation bases are proved exact by independent finite criteria.
Human CM/local-unit argument turns these identities into nonisomorphism.
"""
from pathlib import Path
import argparse,json,hashlib,resource,signal,time,sys
import sympy as s
resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30);start=time.monotonic()
if sys.flags.optimize:raise RuntimeError('Enabled assertions required')
HERE=Path(__file__).resolve().parents[1]/'certificates'
ap=argparse.ArgumentParser()
ap.add_argument('--data',type=Path,default=HERE/'principal_bc_unpolarized_data.json')
ap.add_argument('--c58-data',type=Path,default=HERE/'principal_unmarked_bc_data.json')
ap.add_argument('--c59-data',type=Path,default=HERE/'principal_bc_isogeny_data.json')
args=ap.parse_args()
newraw=args.data.read_bytes()
assert hashlib.sha256(newraw).hexdigest()=='f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e'
new=json.loads(newraw);assert new['schema']=='PRINCIPAL_BC_UNPOLARIZED_V125_1'
assert new['C58_data_sha256']=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
assert new['C59_data_sha256']=='4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee'
raw=args.c58_data.read_bytes();assert hashlib.sha256(raw).hexdigest()=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
old=json.loads(raw);d=old['data']
assert old['actual_tangent_types']=={'B':[1,5,23,25,29,31],'C':[1,5,11,13,17,19]}
assert [g for g in (1,5,11,13,17,19,23,25,29,31,37,41) if {g*i%42 for i in old['actual_tangent_types']['C']}==set(old['actual_tangent_types']['B'])]==[31]
def rref3(rows):
 A=[[int(v)%3 for v in row] for row in rows];r=0;pivots=[]
 for col in range(len(A[0])):
  pivot=next((i for i in range(r,len(A)) if A[i][col]),None)
  if pivot is None:continue
  A[r],A[pivot]=A[pivot],A[r];inv=pow(A[r][col],-1,3);A[r]=[(v*inv)%3 for v in A[r]]
  for i in range(len(A)):
   if i!=r and A[i][col]:
    a=A[i][col];A[i]=[(v-a*w)%3 for v,w in zip(A[i],A[r])]
  pivots.append(col);r+=1
  if r==len(A):break
 return A[:r],pivots

def mod3(M):return M.applyfunc(lambda v:int(v.p)*pow(int(v.q),-1,3)%3)
def rank3(M):return len(rref3(M.tolist())[1])
x=s.Symbol('x');poly=s.Poly(s.cyclotomic_poly(42,x),x);C=s.zeros(12)
for i in range(11):C[i+1,i]=1
for i in range(12):C[i,11]=-poly.nth(i)
e=s.eye(12)[:,0];powers=[C**j for j in range(12)]
def mat(v):return sum((v[j]*powers[j] for j in range(12)),s.zeros(12))
def vec(v):return s.Matrix(v).applyfunc(s.Rational)
conj=s.Matrix.hstack(*[C**(41*j)*e for j in range(12)]);twist=s.Matrix.hstack(*[C**(31*j)*e for j in range(12)])
t=(C+C**41)*e;T=mat(t);OL=s.Matrix.hstack(*[T**j*e for j in range(6)])
records={};Rs={};constraints={}
for name in ('B','C'):
 I=s.Matrix.hstack(*[vec(v) for v in d['I_'+name+'_basis']]);II=I.inv()
 coefficients=s.Matrix.hstack(*[s.Matrix(list(II*P*I)) for P in powers])
 assert coefficients.shape==(144,12) and all(v.q in (1,3) for v in coefficients)
 K=s.Matrix(rref3((3*coefficients).tolist())[0]);assert K.shape==(3,12)
 R=s.Matrix(new['R_C_basis']).applyfunc(s.Rational)
 assert all(v.q==1 for v in R) and abs(R.det())==27
 assert all(v.q==1 for v in coefficients*R)
 assert rank3(R)==9 and all(int(v)%3==0 for v in K*R)
 assert all(v.q==1 for v in R.inv()*OL)
 # Any multiplier is an algebraic integer (integral action on full lattice),
 # hence in OM=Z[z]. Kernel of these3 independent mod3 conditions is complete.
 Rs[name]=R;constraints[name]=K
 records[name]={'constraint_rank_mod3':3,'index_in_OM':27,'constraints_mod3':[[int(v) for v in row] for row in K.tolist()],'basis_matches_full_multiplier_ring_by_equal_index':True}
change=Rs['C'].inv()*twist*Rs['B'];assert all(v.q==1 for v in change) and abs(change.det())==1
assert Rs['B']==Rs['C']
# Conductor: x*OM subset R iff K*C^j*x=0mod3 for allj.
D=s.Matrix.vstack(*[constraints['C']*P for P in powers]);DR=s.Matrix(rref3(D.tolist())[0]);assert DR.shape==(6,12)
h=s.Poly(x**6-x**5+x**4-x**3+x**2-x+1,x,modulus=3)
assert s.rem(s.Poly(poly.as_expr(),x,modulus=3),h*h).is_zero
# Exact irreducibility criterion for degree6, characteristic3.
def power_mod(n):
 base=s.Poly(x,x,modulus=3);a=s.Poly(1,x,modulus=3)
 while n:
  if n&1:a=s.rem(a*base,h)
  base=s.rem(base*base,h);n//=2
 return a
xp=s.Poly(x,x,modulus=3)
assert power_mod(3**6)==xp
assert s.gcd(h,power_mod(3**3)-xp).degree()==0 and s.gcd(h,power_mod(3**2)-xp).degree()==0
reductions=[s.rem(s.Poly(x**j,x,modulus=3),h) for j in range(12)]
P=s.Matrix(6,12,lambda i,j:int(reductions[j].nth(i))%3)
assert s.Matrix(rref3(P.tolist())[0])==DR
# P=conductor ideal, quotient OM/P=F_(3^6), since h irreducible.
# Let f=(3,g(t)) in OL. It has norm27 and f*OM=P.
gt=T**3-T**2+T+s.eye(12);g=s.Poly(x**3-x**2+x+1,x,modulus=3)
assert all(int(g.eval(a))%3 for a in range(3)) # cubic without root is irreducible
assert s.rem(s.Poly(x**6+x**5-6*x**4-6*x**3+8*x**2+8*x+1,x,modulus=3),g*g).is_zero
assert all(int(v)%3==0 for v in P*gt) and rank3(gt)==6
assert rank3(P*OL)==3
# R is preimage of residue OL image: dimensions9 and constraints vanish on OL+P.
assert all(int(v)%3==0 for v in constraints['C']*gt)
assert rank3(s.Matrix.hstack(OL,gt))==9
assert all(int(v)%3==0 for v in P*(conj-s.eye(12))*OL)
# Quadratic Gauss sum sqrt(-7); alpha=(1+s)/2 is in OM with relative norm2.
q=C**6;sqrt7=(q+q**2+q**4-q**3-q**5-q**6)*e
assert mat(sqrt7)*sqrt7==-7*e and conj*sqrt7==-sqrt7
alpha=(e+sqrt7)/2;assert all(v.q==1 for v in alpha)
assert mat(alpha)*conj*alpha==2*e
ratio=mat(alpha)*alpha/2;assert all(int(v)%3==0 for v in mod3(P*(ratio-sqrt7)))
sres=mod3(P*sqrt7);assert sres!=s.zeros(6,1)
assert mod3(P*(mat(sqrt7)*sqrt7))==mod3(P*(-e))
assert mod3(P*(mat(sqrt7)**4*e))==mod3(P*e)
# Reduction of z has exact order14; every root of unity in M is mu42.
assert all(int(v)%3==0 for v in P*(C**14*e-e))
assert all(any(int(v)%3 for v in P*(C**j*e-e)) for j in (1,2,7))
# Source/target conductor is the same and fixed by sigma31 and conjugation.
assert s.Matrix(rref3((P*twist).tolist())[0])==DR
assert s.Matrix(rref3((P*conj).tolist())[0])==DR
# Independently certify proposed saturation bases by containment and mod3 generation.
# This avoids rational row-lattice duality and Hermite normal form.
Ss={};sat_records={}
for name in ('B','C'):
 I=s.Matrix.hstack(*[vec(v) for v in d['I_'+name+'_basis']]);S=s.Matrix(new['S_'+name+'_basis']).applyfunc(s.Rational);SI=S.inv();A=SI*I
 assert all(v.q==1 for v in A) and abs(A.det())==27
 assert all(v.q==1 for v in SI*C*S) and all(v.q==1 for v in 3*A.inv())
 generators=s.Matrix.hstack(*[SI*power*I for power in powers]);assert all(v.q==1 for v in generators)
 assert rank3(generators)==12
 # Generated lattice contains3S (already inI) and maps ontoS/3S, so equalsS.
 Ss[name]=S;sat_records[name]={'contains_actual_I':True,'index_S_over_I':27,'stable_under_z':True,'3S_subset_I':True,'saturation_generators_rank_mod3':12,'exact_saturation_without_HNF':True}
recraw=args.c59_data.read_bytes();assert hashlib.sha256(recraw).hexdigest()=='4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee'
rec=json.loads(recraw);m0=vec(rec['witness_m']);image=mat(m0)*twist*Ss['B'];target=mat(alpha)*Ss['C'];change=target.inv()*image
assert all(v.q==1 for v in change) and change.det()==1
bara=conj*alpha;wrong=(mat(bara)*Ss['C']).inv()*image;assert any(v.q!=1 for v in wrong)
assert alpha==vec(new['alpha']) and sqrt7==vec(new['sqrt_minus7'])
# Prime pi=(1-z^14), agree with earlier mod3-kernel conductor by index and containment.
F=s.eye(12)-C**14;assert F.det()==729 and all(int(v)%3==0 for v in P*F)
# Actual F0 is a local lattice isomorphism at rational3: integral T, inverse denominators powers2.
IB=s.Matrix.hstack(*[vec(v) for v in d['I_B_basis']]);IC=s.Matrix.hstack(*[vec(v) for v in d['I_C_basis']]);T0=IC.inv()*mat(m0)*twist*IB
assert T0==s.Matrix(rec['witness_matrix']) and T0.det()==64 and all(v.q==1 for v in T0)
assert all(int(v.q)%3!=0 for v in T0.inv())
# Every exact multiplier basis is fixed by conjugation moduloP, the only residue property needed.
assert all(int(v)%3==0 for v in P*(conj-s.eye(12))*Rs['C'])
# Recheck exact polarization and local isomorphism of the inherited m0 witness.
E=s.Matrix(d['principal_alternating_matrix']).applyfunc(s.Rational)
assert T0.T*E*T0==2*E
# Contradiction itself is human: eta/bar(eta) is in mu42 but would reduce to this order4 ratio.
assert mod3(P*(mat(ratio)**14*e+e))==s.zeros(6,1)
assert mod3(P*(mat(ratio)**14*e-e))!=s.zeros(6,1)
signal.alarm(0)
print(json.dumps({'status':'INDEPENDENT_PRINCIPAL_BC_UNPOLARIZED_OBSTRUCTION_INPUTS_PASS',
 'method':'Cyclotomic companion arithmetic, mod3 multiplier/conductor kernels and saturation membership/generation; no primary imports or HNF',
 'full_multiplier_order_indices':[27,27],'saturation_indices':[27,27],
 'saturation_identity':'m0*sigma31(SB)=alpha*SC','modulus_index':729,
 'RC_residues_conjugation_fixed':True,'alpha_ratio_residue_order':4,
 'mu42_residue_order':14,'CM_unit_contradiction':'Order4 cannot divide14; human CM-unit and local-isomorphism deductions required',
 'data_sha256':hashlib.sha256(newraw).hexdigest(),'C58_data_sha256':hashlib.sha256(raw).hexdigest(),
 'C59_data_sha256':hashlib.sha256(recraw).hexdigest(),
 'scope':'Specified C56 principal B/C underlying varieties; no claim for arbitrary isogenous models or original trace lattices'},sort_keys=True))
