"""Exact four-pair B/C polarization/isogeny witness, pending review.
Verifies the two primes above2 via irreducible factors, all84 norm2 family
members per pair, six integral Hom survivors, actual degree64 polarization
matrices and complementary kernels. Completeness is the named human proof
in experiments/bc_polarization_square_v141/PROPOSAL.md.
"""
import resource,signal,importlib.util,json,hashlib,sys
from pathlib import Path
if sys.flags.optimize:raise RuntimeError("Enabled assertions are required")
resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
import sympy as s
ROOT=Path(__file__).resolve().parents[1];p=ROOT/'certificates/principal_unmarked_bc_v122.py'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='a73582098061839d7405a1fe00576b8bc50b703d0d830bb7d3977ecad446fdd3'
spec=importlib.util.spec_from_file_location('c58',p);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
K=base.K;z=K([0,1]);zp=[z**i for i in range(12)];t=z+z**41
mat=lambda vs:s.Matrix.hstack(*(s.Matrix(v.c) for v in vs))
mm=lambda v:mat([v*w for w in zp])
def auto(v,g):return sum((c*z**((g*i)%42) for i,c in enumerate(v.c)),K())
bar=lambda v:auto(v,41)
PINS={'principal_unmarked_bc_data.json':'b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f','principal_bc_isogeny_data.json':'4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee','principal_bc_unpolarized_data.json':'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e','principal_bc_hom_data.json':'9f815f1592c8e22a2a14b15086d8d9cb63092f9ab84bb90a3d5ee4935b3fe22c'}
def load(n):
 p=ROOT/'certificates'/n
 assert hashlib.sha256(p.read_bytes()).hexdigest()==PINS[n]
 return json.loads(p.read_text())
x=s.Symbol('x');f1=s.Poly(x**6+x**4+x**2+x+1,x,modulus=2);f2=s.Poly(x**6+x**5+x**4+x**2+1,x,modulus=2)
phi=s.Poly(x**12+x**11-x**9-x**8+x**6-x**4-x**3+x+1,x,modulus=2)
assert f1*f2==phi and f1.is_irreducible and f2.is_irreducible
c58=load('principal_unmarked_bc_data.json')['data'];c59=load('principal_bc_isogeny_data.json');c60=load('principal_bc_unpolarized_data.json');c62=load('principal_bc_hom_data.json')
IB=[K(v) for v in c58['I_B_basis']];IC=[K(v) for v in c58['I_C_basis']];BM=mat(IB);CM=mat(IC);CI=CM.inv()
E=s.Matrix(c58['principal_alternating_matrix']);alpha=K(c60['alpha']);factor=K(c59['witness_m'])*alpha.inverse();J=s.Matrix(c62['J_basis']);JI=J.inv()
gamma=z**6*(1-z**2);q=gamma*bar(gamma);assert q==4-t**2
assert alpha+bar(alpha)==1 and alpha*bar(alpha)==2
assert all(c.denominator==1 for v in [gamma,gamma.inverse(),auto(gamma,31),auto(gamma.inverse(),31)] for c in v.c)
qB=BM.inv()*mm(q)*BM;qC=CI*mm(q)*CM
integral=lambda A:all(v.q==1 for v in A)
enc=lambda M:[[str(v) for v in M.row(i)] for i in range(M.rows)]
models={}
for e in (0,1):
 for f in (0,1):
  u=auto(gamma,31)**e * (gamma.inverse()**f)
  normtarget=2*auto(q,31)**e * (q.inverse()**f)
  hits=[]
  for name,al in [('alpha',alpha),('bar_alpha',bar(alpha))]:
   for j in range(42):
    eta=al*u*z**j
    assert eta*bar(eta)==normtarget
    if integral(JI*s.Matrix(eta.c)):hits.append((name,j,eta))
  chosen='alpha' if e==f else 'bar_alpha'
  assert [(name,j) for name,j,eta in hits]==[(chosen,j) for j in range(0,42,7)]
  eta=hits[0][2];T=CI*mat([factor*eta*auto(v,31) for v in IB])
  assert integral(T) and T.det()==64 and T.T*E*(qC**f)*T==2*E*(qB**e)
  models[str(e)+str(f)]={'source_twist':e,'target_twist':f,'eta':[str(c) for c in eta.c],'branch':chosen,'root_exponents':[j for name,j,eta in hits],'matrix':enc(T),'degree':64}
# On B[2], the kernel consists of x in F2^12 with T*x=0.
def rank2(A):
 a=[[int(v)%2 for v in A.row(i)] for i in range(A.rows)];r=0
 for j in range(A.cols):
  i=next((i for i in range(r,len(a)) if a[i][j]),None)
  if i is None:continue
  a[r],a[i]=a[i],a[r]
  for i in range(len(a)):
   if i!=r and a[i][j]:a[i]=[x^y for x,y in zip(a[i],a[r])]
  r+=1
 return r
mats={key:s.Matrix(v['matrix']) for key,v in models.items()}
assert all(rank2(T)==6 for T in mats.values())
assert rank2(mats['00'].col_join(mats['11']))==6
assert rank2(mats['01'].col_join(mats['10']))==6
assert rank2(mats['00'].col_join(mats['01']))==12
out={'status':'FOUR_BC_POLARIZATION_PAIRS_PASS','models':models,'norm2_candidates_per_pair':84,'integral_Hom_survivors_per_pair':6,'kernel_classes':'equal for00/11 and01/10; two complementary order64 subgroups of B[2]','scope':'Witnesses and exact membership enumeration among supplied norm2 families; completeness/minimality require structural proof.'}
print(json.dumps(out,indent=2,sort_keys=True))
