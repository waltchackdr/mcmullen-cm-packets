"""Small exact finite-field symmetry experiment, NOT a Galois certificate."""
from pathlib import Path
import hashlib,json,math,resource,signal,time
ROOT=Path(__file__).resolve().parents[2]
PINS={'experiments/s_order_models_v161/model.json':'1803891819609ee4c8ecb84222a898a2574207f06bfdf08c48fca95ef27f8184',
'experiments/s_order_models_v161/catalog_candidate.json':'e6a6d431107ce1ceb9cee31ace70558876800d56daf3869a2174dfbea9eb97e3'}
ONE=(1,0,0,0,0,0);T=(0,1,0,0,0,0)
def mul(a,b):
 c=[0]*11
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%3
 for k in range(10,5,-1):
  v=c[k]
  for j in range(6):c[k-6+j]=(c[k-6+j]-v)%3
 return tuple(c[:6])
def power(a,n):
 b=ONE
 while n:
  if n%2:b=mul(b,a)
  a=mul(a,a);n//=2
 return b
def transpose(a):return [list(x) for x in zip(*a)]
def mm(a,b):return [[sum(x*y for x,y in zip(row,col))%3 for col in zip(*b)] for row in a]
def mv(a,v):return tuple(sum(x*y for x,y in zip(row,v))%3 for row in a)
def identity(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def inverse(a):
 n=len(a);rows=[list(row)+eye for row,eye in zip(a,identity(n))]
 for col in range(n):
  p=next(i for i in range(col,n) if rows[i][col]%3);rows[col],rows[p]=rows[p],rows[col]
  u=pow(rows[col][col],-1,3);rows[col]=[u*x%3 for x in rows[col]]
  for i in range(n):
   if i!=col:
    u=rows[i][col];rows[i]=[(x-u*y)%3 for x,y in zip(rows[i],rows[col])]
 assert [r[:n] for r in rows]==identity(n)
 return [r[n:] for r in rows]
def rref(rows):
 a=[[int(x)%3 for x in row] for row in rows];k=0
 for col in range(6):
  p=next((i for i in range(k,len(a)) if a[i][col]),None)
  if p is None:continue
  a[k],a[p]=a[p],a[k];u=pow(a[k][col],-1,3);a[k]=[u*x%3 for x in a[k]]
  for i in range(len(a)):
   if i!=k:
    u=a[i][col];a[i]=[(x-u*y)%3 for x,y in zip(a[i],a[k])]
  k+=1
  if k==len(a):break
 return tuple(tuple(x) for x in a[:k])
def multiply_matrix(a):return transpose([mul(a,power(T,j)) for j in range(6)])
def orbits(permutations):
 unseen=set(range(160));result=[]
 while unseen:
  orbit={min(unseen)};front=list(orbit)
  while front:
   i=front.pop()
   for perm in permutations:
    v=perm[i]
    if v not in orbit:orbit.add(v);front.append(v)
  assert orbit<=unseen;unseen-=orbit;result.append(sorted(orbit))
 return result

def main():
 if not __debug__:raise RuntimeError('Assertions required')
 resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30);begun=time.monotonic()
 data={}
 for path,h in PINS.items():
  raw=(ROOT/path).read_bytes();assert hashlib.sha256(raw).hexdigest()==h;data[path]=json.loads(raw)
 seed=data['experiments/s_order_models_v161/model.json'];catalog=data['experiments/s_order_models_v161/catalog_candidate.json']
 assert power(T,7)==ONE and all(power(T,i)!=ONE for i in range(1,7))
 # Phi7 irreducible over F3 because ord7(3)=6. Multiplicative group has728 elements.
 generator=None
 for code in range(2,729):
  a=tuple((code//3**j)%3 for j in range(6))
  if all(power(a,728//p)!=ONE for p in [2,7,13]):generator=a;break
 assert generator is not None and power(generator,728)==ONE
 root=[[int(x)%3 for x in row] for row in seed['order7_action']]
 q=[[int(x)%3 for x in row] for row in seed['quotient_projection']]
 e=tuple(row[0] for row in q);cols=[];v=e
 for _ in range(6):cols.append(v);v=mv(root,v)
 basis=transpose(cols);inv=inverse(basis)
 assert mm(mm(inv,root),basis)==multiply_matrix(T)
 omega=[[int(x)%3 for x in row] for row in seed['discriminant_gram']]
 form=mm(mm(transpose(basis),omega),basis)
 def similitude(matrix):
  pulled=mm(mm(transpose(matrix),form),matrix)
  return next((u for u in [1,2] if pulled==[[u*x%3 for x in row] for row in form]),None)
 scalar=multiply_matrix(power(generator,13));normone=multiply_matrix(power(generator,26))
 assert similitude(scalar)==2 and similitude(normone)==1
 frob=transpose([power(T,3*j) for j in range(6)])
 corrections=[]
 for exponent in range(26):
  corrected=mm(multiply_matrix(power(generator,exponent)),frob);factor=similitude(corrected)
  if factor is not None:corrections.append((exponent,factor,corrected))
 assert len(corrections)==2
 exponent,factor,corrected=corrections[0]
 # Every catalog root orbit is explicit; no new Grassmann enumeration occurs.
 lookup={tuple(tuple(row) for row in h):m['index'] for m in catalog['models'] for h in m['root_orbit']}
 assert len(lookup)==1120
 def perm(matrix):
  actual=mm(mm(basis,matrix),inv)
  result=[lookup[rref([mv(actual,row) for row in m['lagrangian_rows']])] for m in catalog['models']]
  assert sorted(result)==list(range(160))
  assert all(catalog['models'][i]['endomorphism_order']==catalog['models'][v]['endomorphism_order'] for i,v in enumerate(result))
  return result
 permutations={'norm_one_scalar':perm(normone),'full_similitude_scalar':perm(scalar),'corrected_frobenius':perm(corrected)}
 groups={label:orbits([permutations[k] for k in generators]) for label,generators in {
  'reflex_residue_candidate':['norm_one_scalar'],
  'norm_one_plus_corrected_frobenius_candidate':['norm_one_scalar','corrected_frobenius'],
  'full_semilinear_similitude_candidate':['full_similitude_scalar','corrected_frobenius']}.items()}
 phi=seed['actual_tangent_type'];units=[a for a in range(42) if math.gcd(a,42)==1]
 stabilizer=[a for a in units if sorted(a*x%42 for x in phi)==sorted(phi)]
 logs={pow(3,j,7):j for j in range(6)}
 sums=[sum(3**logs[a%7] for a in phi),sum(3**logs[pow(a,-1,42)%7] for a in phi)]
 assert stabilizer==[1] and sums==[598,286] and all(math.gcd(x,728)==26 for x in sums)
 summary={name:{'orbits':len(os),'R_sizes':[len(o) for o in os if catalog['models'][o[0]]['endomorphism_order']=='R'],'S_sizes':[len(o) for o in os if catalog['models'][o[0]]['endomorphism_order']=='S']} for name,os in groups.items()}
 # Exact closure on the160 already certified classes; at most24 permutations.
 generators=[permutations['full_similitude_scalar'],permutations['corrected_frobenius']]
 closure={tuple(range(160))};front=list(closure)
 while front:
  a=front.pop()
  for b in generators:
   c=tuple(b[x] for x in a)
   if c not in closure:closure.add(c);front.append(c)
  assert len(closure)<=24
 assert len(closure)==24
 assert summary['full_semilinear_similitude_candidate']['R_sizes']==[4]
 assert sorted(summary['full_semilinear_similitude_candidate']['S_sizes'])==[4,8,12,12,12,12,24,24,24,24]
 summary['full_semilinear_similitude_candidate']['permutation_group_order']=len(closure)
 summary['full_semilinear_similitude_candidate']['S_representatives']=[o[0] for o in groups['full_semilinear_similitude_candidate'] if catalog['models'][o[0]]['endomorphism_order']=='S']
 return {'status':'PROVISIONAL_FINITE_SYMMETRY_NOT_ACTUAL_GALOIS_ACTION','pins':PINS,'field_polynomial':[1]*7,'generator':generator,'field_to_actual_basis':basis,'field_form':form,'norm_one_scalar':normone,'full_similitude_scalar':scalar,'corrected_frobenius':corrected,'frobenius_correction_exponent':exponent,'frobenius_similitude':factor,'all_corrections_exponent_factor':[[x,y] for x,y,_ in corrections],'type_stabilizer':stabilizer,'type_and_reflex_power_sums':sums,'residue_power_image_order':28,'permutations':permutations,'orbits':groups,'summary':summary,'wall_seconds':time.monotonic()-begun,'scope':'Exact finite-field/catalog actions only. Actual CM reciprocity image, polarized source descent/normalization and Galois propagation require proof and hostile review. No theta evaluations or exclusions.'}
if __name__=='__main__':
 result=main();print(json.dumps(result,indent=2,sort_keys=True))
