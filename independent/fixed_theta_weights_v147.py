"""Separate quotient/quadratic-form count; no primary weight function or nullbasis."""
import hashlib, itertools, json, resource, signal, time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(20,20)); signal.alarm(30)
start=time.monotonic()
root=Path(__file__).resolve().parents[1]/'experiments/bc_jacobian_v147'
def rref(A):
 A=[[x%3 for x in r] for r in A]; piv=[]; row=0
 for j in range(len(A[0])):
  p=next((i for i in range(row,len(A)) if A[i][j]),None)
  if p is None: continue
  A[row],A[p]=A[p],A[row]
  c=pow(A[row][j],-1,3); A[row]=[c*x%3 for x in A[row]]
  for i in range(len(A)):
   if i!=row:
    c=A[i][j]; A[i]=[(a-c*b)%3 for a,b in zip(A[i],A[row])]
  piv.append(j); row+=1
  if row==len(A): break
 return piv

def determinant(A):
 A=[[x%3 for x in r] for r in A]; d=1
 for j in range(len(A)):
  p=next((i for i in range(j,len(A)) if A[i][j]),None)
  if p is None:return 0
  if p!=j:A[j],A[p]=A[p],A[j];d=-d
  a=A[j][j];d=d*a%3
  for i in range(j+1,len(A)):
   c=A[i][j]*pow(a,-1,3)%3
   A[i]=[(u-c*v)%3 for u,v in zip(A[i],A[j])]
 return d%3

def apply(F,x):return [sum(a*b for a,b in zip(row,x)) for row in F]
def pairing(x,y):return sum(x[i]*y[i+6]-x[i+6]*y[i] for i in range(6))
setup=json.loads((root/'period_setup.json').read_text())
expected=json.loads((root/'fixed_weights_exact.json').read_text())
answer={}
for name,m in setup['models'].items():
 F=[[int(x) for x in row] for row in m['delta_in_symplectic_basis']]
 unit=[[int(i==j) for i in range(12)] for j in range(12)]
 for e in unit:
  fe=apply(F,e); ffe=apply(F,fe)
  assert all(a+b+c==0 for a,b,c in zip(e,fe,ffe))
 for e in unit:
  for f in unit:assert pairing(apply(F,e),apply(F,f))==pairing(e,f)
 assert [pairing(e,apply(F,e))%2 for e in unit]==m['invariant_quadratic_refinement_chi']
 # x=-(F+2I)k/3 gives (F-I)x=k. Select independent image columns.
 V=[[-F[i][j]-2*int(i==j) for j in range(12)] for i in range(12)]
 selected=rref(V);assert len(selected)==6
 Q=[[2*(pairing(unit[i],apply(F,unit[j]))+pairing(unit[j],apply(F,unit[i])))%3 for j in selected] for i in selected]
 counts=[0]*3; points=set(); pairs=set()
 for a in itertools.product(range(3),repeat=6):
  k=[0]*12
  for j,c in zip(selected,a):k[j]=c
  v=tuple(x%3 for x in apply(V,k));assert v not in points;points.add(v)
  b=pairing(k,apply(F,k))%3;counts[b]+=1
  assert all(x%3==0 for x in [u-vv for u,vv in zip(apply(F,v),v)])
  if b==0 and any(v):pairs.add(min(v,tuple((-x)%3 for x in v)))
 assert sorted(pairs)==[tuple(v) for v in expected['models'][name]['weight_zero_nonzero_pairs']]
 assert counts==[261,234,234]
 answer[name]={'coset_coordinate_indices_zero_based':selected,'quadratic_matrix_mod3':Q,'rank':len(rref(Q)),'det_mod3':determinant(Q),'counts_b_mod3':counts,'eligible_pairs':len(pairs),'exact_pair_list_match':True}
print(json.dumps({'status':'PASS_FIXED_WEIGHTS_ONLY_NO_THETA_VALUES','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ['period_setup.json','fixed_weights_exact.json']},'method':'quotient k modulo (F-I)Lambda; b(k)=E(k,Fk) mod3','models':answer,'elapsed_seconds':time.monotonic()-start},indent=2))
