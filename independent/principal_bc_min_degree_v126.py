"""Independent V126 actual-lattice companion verification; no producer imports.
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
parser.add_argument('--data',type=Path,default=HERE/'principal_bc_min_degree_data.json')
parser.add_argument('--c58-data',type=Path,default=HERE/'principal_unmarked_bc_data.json')
parser.add_argument('--c59-data',type=Path,default=HERE/'principal_bc_isogeny_data.json')
parser.add_argument('--c60-data',type=Path,default=HERE/'principal_bc_unpolarized_data.json')
args=parser.parse_args()
new=load(args.data,'1e6b331975adac3175523d6a0a1de7152e4dc230be440e7e96a71d2243e95a6a')
assert new['schema']=='PRINCIPAL_BC_MIN_DEGREE_V126_1'
record=load(args.c58_data,'b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f');d=record['data']
c59=load(args.c59_data,'4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee')
c60=load(args.c60_data,'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e')
for name,p in [('C58',args.c58_data),('C59',args.c59_data),('C60',args.c60_data)]:
 assert new[name+'_data_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
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
beta7=(I+C-C**3)*e;assert beta7==V(new['beta_coefficients']);assert new['root_exponent']==5;eta=C**5*beta7
assert M(beta7).det()==7 and M(eta).det()==7
m=M(eta)*M(alpha).inv()*m0
assert m==V(new['m7'])==s.Matrix([-40,-41,-55,-9,62,46,31,23,12,53,68,-51])/127
A=ideals['C'].inv()*M(m)*twist*ideals['B']
assert A==s.Matrix(new['matrix']) and all(v.q==1 for v in A) and A.det()==7
assert all(v.q==1 for v in 7*A.inv())
SNF=smith_normal_form(A,domain=ZZ)
assert [abs(int(SNF[i,i])) for i in range(12)]==[1]*11+[7]
E=forms['B'];pullback=A.T*forms['C']*A
ratios={pullback[i,j]/E[i,j] for i in range(12) for j in range(12) if E[i,j]}
assert len(ratios)>1
# Equivalent Hermitian obstruction: m*bar(m)/r is not rational scalar.
hermitian=M(inv(r))*M(m)*conj*m
assert any(hermitian[i] for i in range(1,12))
# Direct enumeration in F4, F3 and F5, no polynomial gcd/factor routine.
# F4=F2[a]/(a^2+a+1); integers 0..3 encode constant/linear bits.
def mul4(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a&4:a^=7
 return out
coeff=[int(v) for v in f.all_coeffs()]
values4=[]
for v in range(4):
 out=0
 for c in coeff:out=mul4(out,v)^(c%2)
 values4.append(out)
assert all(values4) # excludes roots overF2 and irreducible quadratic factors.
prime_roots={}
for p in (3,5):
 values=[]
 for v in range(p):
  out=0
  for c in coeff:out=(out*v+c)%p
  values.append(out)
 assert all(values);prime_roots[str(p)]=values
# OM=Z[z] inherited: prime ideal quotients of size2,4 embed inF4,
# and those of size3,5 areF3,F5. Each would give an excluded root.
assert new['expected_degree']==7 and new['kernel']=='Z/7Z'
signal.alarm(0)
result={'status':'INDEPENDENT_PRINCIPAL_BC_MINIMUM_DEGREE7_INPUTS_PASS','seconds':time.monotonic()-start,'degree':7,'kernel':'Z/7Z','actual_lattices_and_forms_reconstructed':True,'unique_hodge_twist':31,'eta_formula':'z^5*(1+z-z^3)','eta_norm':7,'saturation_identity_rechecked':True,'not_scalar_polarized':True,'finite_field_evaluations':{'F4':values4,**prime_roots},'excluded_prime_ideal_norms':[2,3,4,5],'data_sha256':hashlib.sha256(args.data.read_bytes()).hexdigest(),'map_count_claimed':False,'human_dependencies':'C58 unique twist and actual geometric lattices; C60 exact saturations/nonisomorphism; integral action implies eta in OM; equal saturation indices imply degree=Neta; small finite-field root exclusions rule out degrees2..6.'}
print(json.dumps(result,sort_keys=True))
