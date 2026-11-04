"""Exact full induced-A polarization obstruction, pending review/adoption.
Reconstructs the full M+M sigma29 endomorphism lattice and exhausts its
necessary trace12 ball. Geometric model/full-End proof is HUMAN-PROOF in
experiments/a_unmarked_polarization_v140/PROPOSAL.md. Reuses C58 arithmetic.
"""
import resource,signal,importlib.util,hashlib,json,time,math,sys
from pathlib import Path
from fractions import Fraction as F
if sys.flags.optimize:raise RuntimeError("Enabled assertions are required")
resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form
ROOT=Path(__file__).resolve().parents[1];start=time.monotonic()
p=ROOT/'certificates/principal_unmarked_bc_v122.py'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='a73582098061839d7405a1fe00576b8bc50b703d0d830bb7d3977ecad446fdd3'
spec=importlib.util.spec_from_file_location('c58',p);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
K=base.K;z=K([0,1]);zp=[z**i for i in range(12)];t=z+z**41
mat=lambda vs:sp.Matrix.hstack(*(sp.Matrix(v.c) for v in vs))
mm=lambda v:mat([v*w for w in zp])
def auto(v,g):return sum((c*z**((g*i)%42) for i,c in enumerate(v.c)),K())
bar=lambda v:auto(v,41)
trace_coeff=[sp.trace(mm(v)) for v in zp]
tr=lambda v:sum(c*a for c,a in zip(v.c,trace_coeff))
cs=lambda v:[str(c) for c in v.c]
enc=lambda M:[[str(v) for v in M.row(i)] for i in range(M.rows)]
q=4-t**2;u=t-1;beta=-t**5+t**4+6*t**3-2*t**2-8*t-2
a=(beta*u**2).inverse();b=(u**3).inverse();delta=a*b
c=(t**3-1)*t.inverse()*(t**2-1).inverse()
d=(-t**5+5*t**3+t**2-5*t)*(t**2-1).inverse()
w=(z-d)*c.inverse();v1=a*w;v2=b*(1+3*w)
I=[v1*t**i for i in range(6)]+[v2*t**i for i in range(6)]
IM=mat(I);IMi=IM.inv();xi=-c*(z-z**41).inverse()
E=sp.Matrix(12,12,lambda i,j:tr(xi*I[i]*bar(I[j])))
data58=ROOT/'certificates/principal_unmarked_bc_data.json'
assert hashlib.sha256(data58.read_bytes()).hexdigest()=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
target=sp.Matrix(json.loads(data58.read_text())['data']['principal_alternating_matrix'])
assert E==target and E.det()==1
Phi={1,5,17,19,29,31};assert {29*i%42 for i in Phi}==Phi
sigmaI=[auto(v,29) for v in I]

units=[j for j in range(42) if math.gcd(j,42)==1]
assert [j for j in units if {j*x%42 for x in Phi}==Phi]==[1,29]
# Literal selected A generator, not an arbitrary field presentation.
AA=[[t**5-4*t**3-t**2+3*t,t*(-t**5+t**4+6*t**3-3*t**2-10*t)],[(t**3-1)*t.inverse(),-t**5+5*t**3+t**2-6*t]]
ZZ=[[(AA[i][j]+(t if i==j else 0))*(t**2-1).inverse() for j in range(2)] for i in range(2)]
assert ZZ[1][0]==c and ZZ[1][1]==d
assert ZZ[0][0]+ZZ[1][1]==t and ZZ[0][0]*ZZ[1][1]-ZZ[0][1]*ZZ[1][0]==1
assert w*ZZ[0][0]+ZZ[0][1]==z*w and w*ZZ[1][0]+ZZ[1][1]==z
# The actual Hodge endomorphism coefficients: a+b*sigma29.
ops=[IMi*mat([v*x for x in X]) for X in (I,sigmaI) for v in zp]
constraints=sp.Matrix.hstack(*(sp.Matrix(list(U)) for U in ops))
den=sp.ilcm(*[v.q for v in constraints]);H=hermite_normal_form((den*constraints).applyfunc(int).T);L=den*H.T.inv()
assert L.shape==(24,24) and all(v.q==1 for v in constraints*L)
matrices=[sum((L[i,j]*ops[i] for i in range(24)),sp.zeros(12)) for j in range(24)]
qmatrix=IMi*mat([q*v for v in I]);assert qmatrix.T*E==E*qmatrix
qinv=qmatrix.inv();Einv=E.inv();adj=[U.T for U in matrices]
weighted=[qinv*Einv*U.T*E for U in matrices]
Q=sp.Matrix(24,24,lambda i,j:sum(weighted[i][a,b]*matrices[j][b,a] for a in range(12) for b in range(12)))
assert Q==Q.T

record=json.loads(Path(__file__).with_name('a_unmarked_polarization_data.json').read_text())
assert L==sp.Matrix(record['coefficient_basis_a_then_b']) and Q==sp.Matrix(record['Gram_Q'])
U=sp.Matrix(record['U']);reduced=sp.Matrix(record['reduced_Gram']);mu=sp.Matrix(record['mu']);ds=list(map(sp.Rational,record['diagonal']))
assert all(v.q==1 for v in U) and abs(U.det())==1 and U.T*Q*U==reduced
assert all(mu[i,j]==0 for i in range(24) for j in range(i,24))
assert all(v>0 for v in ds) and (sp.eye(24)+mu)*sp.diag(*ds)*(sp.eye(24)+mu).T==reduced
bounds=record['coordinate_bounds'];inverse=reduced.inv()
assert len(bounds)==24 and all(type(v) is int and v>=0 for v in bounds)
assert all((v+1)**2>12*inverse[i,i] for i,v in enumerate(bounds))
vector=[0]*24;nodes=0;leaves=[]
mu=[[F(v) for v in row] for row in mu.tolist()];ds=list(map(F,ds))
def visit(i,remaining):
 global nodes
 nodes+=1
 if nodes>5000:raise RuntimeError('Bounded verification node cap exceeded')
 if i<0:
  assert (sp.Matrix(vector).T*reduced*sp.Matrix(vector))[0]==12-remaining
  leaves.append((vector[:],12-remaining));return
 shift=sum(mu[j][i]*vector[j] for j in range(i+1,24))
 for v in range(-bounds[i],bounds[i]+1):
  vector[i]=v;cost=ds[i]*(F(v)+shift)**2
  if cost<=remaining:visit(i-1,remaining-cost)
visit(23,F(12))
assert nodes==347 and leaves==[([0]*24,F(0))]
print(json.dumps({'status':'A_FULL_ENDOMORPHISM_TRACE12_OBSTRUCTION_PASS','dimension':24,'CM_type_stabilizer':[1,29],'nodes':nodes,'nonzero_closed_ball_vectors':0,'isomorphism_trace_required':12,'scope':'Exact full-lattice finite obstruction; actual geometric/Hodge identification remains the named human-proof input.'},sort_keys=True))
