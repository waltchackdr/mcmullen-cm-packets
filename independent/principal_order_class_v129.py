"""Independent finite Picard quotient inputs: trial divisors and full field enumeration.
No SymPy or primary-certificate imports. Human classnumber/ideal/unit proof separate.
"""
import argparse,hashlib,itertools,json,resource,signal,sys
from fractions import Fraction
from pathlib import Path
PINS={'principal_bc_unpolarized_data.json':'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e','principal_bc_hom_data.json':'9f815f1592c8e22a2a14b15086d8d9cb63092f9ab84bb90a3d5ee4935b3fe22c'}
def check(directory):
 data={}
 for n,h in PINS.items():
  raw=(directory/n).read_bytes();assert hashlib.sha256(raw).hexdigest()==h;data[n]=json.loads(raw)
 def trim(p):
  p=[int(x)%3 for x in p]
  while p and not p[-1]:p.pop()
  return p
 def remainder(p,q):
  p=trim(p);q=trim(q)
  while len(p)>=len(q):
   c=p[-1]*pow(q[-1],-1,3)%3;k=len(p)-len(q)
   for i,v in enumerate(q):p[i+k]=(p[i+k]-c*v)%3
   p=trim(p)
  return p
 h=[1,2,1,2,1,2,1]
 # Every reducible monic degree6 polynomial has a monic factor ofdegree<=3.
 trials=0
 for degree in (1,2,3):
  for a in itertools.product(range(3),repeat=degree):
   assert remainder(h,list(a)+[1]);trials+=1
 assert trials==39
 zero=(0,)*6;one=(1,0,0,0,0,0);z=(0,1,0,0,0,0)
 def canon(p):return tuple((remainder(p,h)+[0]*6)[:6])
 def add(a,b):return canon([x+y for x,y in zip(a,b)])
 def neg(a):return canon([-x for x in a])
 def mul(a,b):
  c=[0]*11
  for i,x in enumerate(a):
   for j,y in enumerate(b):c[i+j]+=x*y
  return canon(c)
 def power(a,n):
  b=one
  while n:
   if n%2:b=mul(b,a)
   a=mul(a,a);n//=2
  return b
 def ratio(a):return mul(a,power(power(a,27),727))
 def order(a):
  b=one
  for i in range(1,729):
   b=mul(b,a)
   if b==one:return i
  raise AssertionError('not a unit')
 assert order(z)==14
 s=zero
 for e in (6,12,24):s=add(s,power(z,e))
 for e in (18,30,36):s=add(s,neg(power(z,e)))
 alpha=mul((2,0,0,0,0,0),add(one,s));gamma=add(power(z,6),neg(power(z,8)))
 assert mul(s,s)==neg(one) and ratio(alpha)==s
 assert ratio(gamma)==neg(one) and order(ratio(z))==7
 assert order(ratio(alpha))==4 and order(ratio(gamma))==2
 image={mul(power(ratio(z),i),power(ratio(gamma),j)) for i in range(7) for j in range(2)}
 assert len(image)==14 and ratio(alpha) not in image
 allfield=list(itertools.product(range(3),repeat=6));fixed={v for v in allfield if power(v,27)==v};normone={v for v in allfield if power(v,28)==one};ratios={ratio(v) for v in allfield if v!=zero}
 assert len(fixed)==27 and len(normone)==len(ratios)==28 and ratios==normone
 assert all(power(q,14)==one for q in image)
 # Rational coefficient arithmetic provides separate global norm checks.
 phi=[1,1,0,-1,-1,0,1,0,-1,-1,0,1,1]
 def zmul(a,b):
  c=[0]*(len(a)+len(b)-1)
  for i,v in enumerate(a):
   for j,w in enumerate(b):c[i+j]+=v*w
  return c
 def binomial_product(degrees):
  a=[1]
  for n in degrees:a=zmul(a,[-1]+[0]*(n-1)+[1])
  return a
 # Mobius cyclotomic identity, checked over Z rather than assumed coefficients.
 assert zmul(phi,binomial_product((21,14,6,1)))==binomial_product((42,7,3,2))
 assert [v%3 for v in phi]==[v%3 for v in zmul(h,h)]
 def qm(a,b):
  c=[Fraction(0)]*23
  for i,v in enumerate(a):
   for j,w in enumerate(b):c[i+j]+=v*w
  for i in range(22,11,-1):
   for j in range(12):c[i-12+j]-=c[i]*phi[j]
  return tuple(c[:12])
 qone=(Fraction(1),)+(Fraction(0),)*11;qz=(Fraction(0),Fraction(1))+(Fraction(0),)*10
 def qp(a,n):
  b=qone
  while n:
   if n%2:b=qm(b,a)
   a=qm(a,a);n//=2
  return b
 def determinant(A):
  A=[list(row) for row in A];out=Fraction(1)
  for k in range(12):
   pivot=next((i for i in range(k,12) if A[i][k]),None)
   if pivot is None:return Fraction(0)
   if pivot!=k:A[k],A[pivot]=A[pivot],A[k];out=-out
   d=A[k][k];out*=d
   for i in range(k+1,12):
    c=A[i][k]/d
    for j in range(k+1,12):A[i][j]-=c*A[k][j]
  return out
 def norm(a):
  cols=[qm(a,tuple(Fraction(i==j) for i in range(12))) for j in range(12)]
  return determinant(list(zip(*cols)))
 assert qp(qz,42)==qone and all(qp(qz,n)!=qone for n in (21,14,6))
 qgamma=tuple(a-b for a,b in zip(qp(qz,6),qp(qz,8)))
 qalpha=tuple(Fraction(v) for v in data['principal_bc_unpolarized_data.json']['alpha'])
 qpi=tuple(a-b for a,b in zip(qone,qp(qz,14)))
 assert norm(qgamma)==1 and norm(qalpha)==64 and norm(qpi)==729
 assert qgamma==tuple(Fraction(v) for v in data['principal_bc_hom_data.json']['square_generator'])
 def reduce_coeffs(coeffs):
  return canon([int(Fraction(v).numerator)*pow(int(Fraction(v).denominator),-1,3)%3 for v in coeffs])
 assert alpha==reduce_coeffs(data['principal_bc_unpolarized_data.json']['alpha'])
 assert gamma==reduce_coeffs(data['principal_bc_hom_data.json']['square_generator'])
 return {'status':'INDEPENDENT_PRINCIPAL_ORDER_PICARD_FINITE_INPUTS_PASS','irreducibility_trial_divisors':trials,'field_size':729,'fixed_field_size':27,'quotient_order':28,'z_unit_ratio_order':7,'gamma_unit_ratio_order':2,'generated_unit_image_order':14,'alpha_ratio_order':4,'alpha_outside_unit_image':True,'conditional_Picard_order':2,'global_norms':{'gamma':1,'alpha':64,'conductor_generator':729},'input_pins':PINS,'human_requirements':['Pic(OM)=1 exact cited theorem','All global CM unit ratios are roots in mu42, reductionorder divides14','Classification of invertible R ideals by normalized nonzero residue lines'],'scope':'Finite arithmetic only; full Picard theorem also uses stated human proof'}
def main():
 if sys.flags.optimize:raise RuntimeError('Optimized execution forbidden')
 resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
 ap=argparse.ArgumentParser();ap.add_argument('--certificates',type=Path,default=Path(__file__).resolve().parents[1]/'certificates');a=ap.parse_args();print(json.dumps(check(a.certificates),indent=2,sort_keys=True))
if __name__=='__main__':main()
