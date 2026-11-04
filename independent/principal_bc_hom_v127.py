"""Independent V127 complete Hom ideal, inverse/class-square and pullback inputs.
Reuses our adopted C59 companion reconstruction, not its primary arithmetic.
"""
import argparse,hashlib,json,math,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
assert __debug__
resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30);start=time.monotonic()
def load(p,h):
 raw=Path(p).read_bytes();assert hashlib.sha256(raw).hexdigest()==h;return json.loads(raw)
HERE=Path(__file__).resolve().parents[1]/'certificates'
parser=argparse.ArgumentParser()
parser.add_argument('--data',type=Path,default=HERE/'principal_bc_hom_data.json')
parser.add_argument('--c58-data',type=Path,default=HERE/'principal_unmarked_bc_data.json')
parser.add_argument('--c59-data',type=Path,default=HERE/'principal_bc_isogeny_data.json')
parser.add_argument('--c60-data',type=Path,default=HERE/'principal_bc_unpolarized_data.json')
parser.add_argument('--c61-data',type=Path,default=HERE/'principal_bc_min_degree_data.json')
args=parser.parse_args()
new=load(args.data,'9f815f1592c8e22a2a14b15086d8d9cb63092f9ab84bb90a3d5ee4935b3fe22c')
assert new['schema']=='PRINCIPAL_BC_HOM_IDEAL_V127_1'
record=load(args.c58_data,'b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f');d=record['data']
c59=load(args.c59_data,'4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee')
c60=load(args.c60_data,'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e')
c61=load(args.c61_data,'1e6b331975adac3175523d6a0a1de7152e4dc230be440e7e96a71d2243e95a6a')
for name,p in [('principal_unmarked_bc_data.json',args.c58_data),('principal_bc_isogeny_data.json',args.c59_data),('principal_bc_unpolarized_data.json',args.c60_data),('principal_bc_min_degree_data.json',args.c61_data)]:
 assert new['input_pins'][name]==hashlib.sha256(p.read_bytes()).hexdigest()
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

s7=(C**6+C**12+C**24-C**18-C**30-C**36)*e
assert M(s7)*s7==-7*e and conj*s7==-s7
alpha=(e+s7)/2;assert alpha==V(c60['alpha']) and M(alpha).det()==64
m0=V(c59['witness_m']);A0=ideals['C'].inv()*M(m0)*twist*ideals['B']
assert A0==s.Matrix(c59['witness_matrix']) and A0.det()==64
SB=s.Matrix(c60['S_B_basis']).applyfunc(s.Rational);SC=s.Matrix(c60['S_C_basis']).applyfunc(s.Rational)
U=(M(alpha)*SC).inv()*M(m0)*twist*SB
assert all(v.q==1 for v in U) and U.det()==1
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

J=s.Matrix(new['J_basis']).applyfunc(s.Rational);JI=J.inv();R=s.Matrix(c60['R_C_basis']).applyfunc(s.Rational)
assert all(v.q==1 for v in J) and J.det()==27
OL=s.Matrix.hstack(*[T**j*e for j in range(6)])
Pi=I-C**14;pi=Pi*e
G=s.Matrix.hstack(M(alpha)*OL,Pi)
assert all(v.q==1 for v in JI*G)
assert all(v.q==1 for v in JI*(3*I)) and rank3(G)==9
# G contains3OM because pi*bar(pi)=3 and bar(pi) is integral.
assert Pi*M(conj*pi)==3*I
# Since generated G contains3OM and has9-dimensional image modulo3,
# its index is27, hence G generates exactly the supplied J.
factor=M(alpha).inv()*m0
assert M(factor)*conj*factor==r
Btw=twist*ideals['B'];CI=ideals['C'].inv()
forms_eta=s.Matrix.hstack(*[s.Matrix(list(CI*M(factor)*P*Btw)) for P in powers])
assert forms_eta.shape==(144,12) and all(v.q in (1,3) for v in forms_eta)
constraints=s.Matrix(rref3((3*forms_eta).tolist())[0]);assert constraints.shape==(3,12)
assert all(v.q==1 for v in forms_eta*J)
# All rational integral maps normalize to algebraic-integer eta by the
# accepted saturation bridge; these3 independent mod3 constraints are exhaustive.
h=s.Poly(x**6-x**5+x**4-x**3+x**2-x+1,x,modulus=3)
residues=[s.rem(s.Poly(x**j,x,modulus=3),h) for j in range(12)]
P=s.Matrix(6,12,lambda i,j:int(residues[j].nth(i))%3)
fixed=mod3(P*(conj-I)*M(alpha).inv())
assert rref3(fixed.tolist())[0]==rref3(constraints.tolist())[0]
# This proves eta modP lies in alpha F27 using the inherited fixed-field identification.
assert rank3(P*OL)==3 and rank3(P*M(alpha)*OL)==3
assert all(v.q==1 for v in JI*s.Matrix.hstack(*[M(R[:,i])*J for i in range(12)]))
products=s.Matrix.hstack(*[M(J[:,i])*J for i in range(12)])
mixed=s.Matrix.hstack(*[M(J[:,i])*conj*J for i in range(12)])
assert all(v.q==1 for v in R.inv()*mixed) and rank3(mixed)==9
# JJbar contains pi*barpi*OM=3OM. Inclusion inR and equal mod3 rank
# therefore proves exact equality JJbar=R, soJ is invertible with inversebarJ.
gamma=(C**6-C**8)*e;assert gamma==V(new['square_generator']);GM=M(gamma);assert GM.det()==1 and all(v.q==1 for v in GM.inv())
assert all(v.q==1 for v in (GM*R).inv()*products) and rank3(products)==9
assert all(v.q==1 for v in (Pi*Pi/3).inv()) and (Pi*Pi/3).det()==1
# JJ contains pi²OM=3OM; gammaR also contains3OM. Both mod3rank9,
# inclusion proves J²=gammaR without HNF or row-lattice dualization.
# Polarization is pulled back to SOURCE coordinates, requiring inverse twist19.
untwist=sigma(19);assert untwist*twist==I
# General identity: untwist(xi_C * factor*barfactor)=xi_B.
assert untwist*(M(xis['C'])*M(factor)*conj*factor)==xis['B']
# Coefficient comparison of the quadratic identity for symbolic
# eta=sum X_j J_j, including every mixed X_jX_k term.
source_matrices=[M(M(factor)*J[:,j]) for j in range(12)]
Jmatrices=[M(J[:,j]) for j in range(12)]
for j in range(12):
 for k in range(j,12):
  target_coefficient=source_matrices[j]*conj*M(factor)*J[:,k]
  eta_coefficient=Jmatrices[j]*conj*J[:,k]
  if k!=j:
   target_coefficient+=source_matrices[k]*conj*M(factor)*J[:,j]
   eta_coefficient+=Jmatrices[k]*conj*J[:,j]
  assert untwist*M(xis['C'])*target_coefficient==M(xis['B'])*untwist*eta_coefficient
for j in range(12):
 eta=J[:,j];m=M(factor)*eta;A=CI*M(m)*Btw
 assert all(v.q==1 for v in A)
 multiplier=untwist*M(eta)*conj*eta
 expected=s.Matrix(12,12,lambda k,l:(tr*M(xis['B'])*M(multiplier)*M(ideals['B'][:,k])*conj*ideals['B'][:,l])[0])
 assert A.T*forms['C']*A==expected
m7=V(c61['m7']);T7=CI*M(m7)*Btw;assert T7==s.Matrix(c61['matrix']) and T7.det()==7
SNF=smith_normal_form(T7.T*forms['C']*T7,domain=ZZ)
assert [abs(int(SNF[i,i])) for i in range(12)]==[1]*10+[7]*2
assert new['source_pullback_twist']==19 and new['degree7_polarization_type']==[1,1,1,1,1,7]
signal.alarm(0)
print(json.dumps({'status':'INDEPENDENT_PRINCIPAL_BC_HOM_IDEAL_AND_PULLBACK_INPUTS_PASS','seconds':time.monotonic()-start,'data_sha256':hashlib.sha256(args.data.read_bytes()).hexdigest(),'method':'Companion actual lattices, all144 normalized integrality forms modulo3, fixed-residue kernel, membership plus equal mod3 dimensions; no HNF or row-dual algorithms','constraint_rank_mod3':3,'J_index_in_OM':27,'full_normalized_Hom_equals_J':True,'J_equals':'alpha OL + pi OM','J_inverse':'bar(J)','J_times_conjugate':'R','J_square':'z^6(1-z^2)R','square_generator_norm':1,'picard_class_order':'Exactly2 using adopted C60 nonisomorphism to exclude principalJ','source_pullback_multiplier':'sigma19(eta*bar(eta))','source_pullback_checked_on_12_generators':True,'symbolic_quadratic_coefficients_checked':78,'universal_source_identity':'sigma19(xi_C*(m0/alpha)*bar(m0/alpha))=xi_B','degree':'Norm_M/Q(eta)','degree7_pullback_polarization_type':[1,1,1,1,1,7],'no_surjectivity_onto_totally_positive_elements_claimed':True,'no_total_Picard_group_order_claimed':True},indent=2,sort_keys=True))
