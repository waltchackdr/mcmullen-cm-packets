#!/usr/bin/env python3
"""Exact v57 common-level stabilizer, curve orbit, and CM-order certificate.

Completeness and family descent: docs/PRINCIPAL_DESCENT_V57.md.
Depends on primary principal_homology_v56; independent path uses its own library.
The finite-domain computation is on 28*8 residue line pairs, not a heuristic.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import principal_homology_v56 as base
import sympy as sp
from itertools import product
from collections import deque

t=base.t; f=base.f
# Reuse the primary v56 field identities, but derive the new ideal identity.
a14,b14=base.zpower(14)
assert abs(base.multiplication(b14).det())==27
G=sp.Poly(t**3-t**2+t+1,t,modulus=3)
assert sp.Poly(b14,t,modulus=3).rem(G).is_zero
assert G.is_irreducible
print('z14 coefficient b =',b14,'; N(b) =',base.multiplication(b14).det(),'; (b)=p3')

# Small prime-field polynomial arithmetic for the projective monodromy orbit.
# Each tuple stores constant-first coefficients of F27=F3[t]/(t^3-t^2+t+1).
zero=(0,0,0);one=(1,0,0)
def add(x,y):return tuple((a+b)%3 for a,b in zip(x,y))
def mul(x,y):
 c=[0]*5
 for i,a in enumerate(x):
  for j,b in enumerate(y):c[i+j]=(c[i+j]+a*b)%3
 for k in (4,3):
  # t^3=t^2-t-1
  for j,coef in enumerate((-1,-1,1)):c[k-3+j]=(c[k-3+j]+c[k]*coef)%3
 return tuple(c[:3])
def power(x,n):
 r=one
 for _ in range(n):r=mul(r,x)
 return r
field=list(product(range(3),repeat=3))
assert all(power(x,26)==one for x in field if x!=zero)
inverse={x:power(x,25) for x in field if x!=zero}
def reduce3(x):
 p=sp.Poly(base.red(x),t,modulus=3).rem(G)
 return tuple(int(p.nth(i))%3 for i in range(3))

def action3(M,m):
 x,y=(zero,one) if m is None else (one,m)
 a=add(mul(M[0][0],x),mul(M[0][1],y))
 b=add(mul(M[1][0],x),mul(M[1][1],y))
 return None if a==zero else mul(b,inverse[a])
def action7(M,m):
 x,y=(0,1) if m is None else (1,m)
 a=(M[0][0]*x+M[0][1]*y)%7;b=(M[1][0]*x+M[1][1]*y)%7
 return None if a==0 else b*pow(a,-1,7)%7
P3=field+[None];P7=list(range(7))+[None]
# Structural transitivity: lower translations and the quarter-turn generate
# all residue lines. Exhaustive SL2(F7) count is a bounded independent check.
sl7=[((a,b),(c,d)) for a,b,c,d in product(range(7),repeat=4) if (a*d-b*c)%7==1]
assert len(sl7)==336
assert len({action7(m,5) for m in sl7})==8
assert sum(action7(m,5)==5 for m in sl7)==42
trans3=[((one,zero),(x,one)) for x in field]
quarter3=((zero,tuple(-v%3 for v in one)),(one,zero))
assert {action3(m,zero) for m in trans3+[quarter3]}==set(P3)
print('ambient index = (27+1)*(7+1) =',28*8)

# McMullen section4 has ABCD=I; the fourth printed generator must agree.
D=sp.Matrix([[-t**5+2*t**4+7*t**3-6*t**2-13*t-1,2*t**5+t**4-9*t**3-5*t**2+6*t],
             [-t**5+5*t**3+2*t**2-5*t+1,t**5-2*t**4-7*t**3+6*t**2+14*t+1]])
D=base.matred(base.normalizer*D*base.normalizer_inv)
ABC=sp.eye(2)
for name in 'ABC':ABC=base.matred(ABC*base.vertices[name])
assert base.matred(ABC*D)==sp.eye(2)

# rho(g)=g^{-T}, the actual transition action (not the fixed-point action g^T).
assert all(base.red(v.det())==1 for v in base.vertices.values())
gens={name:base.matred(v.adjugate().T) for name,v in base.vertices.items()}
local={name:([[reduce3(m[i,j]) for j in range(2)] for i in range(2)],
             [[base.local(m[i,j])[0] for j in range(2)] for i in range(2)])
       for name,m in gens.items()}
space=list(product(P3,P7))
def act(name,p):
 m3,m7=local[name];return action3(m3,p[0]),action7(m7,p[1])
for name in local:assert len({act(name,p) for p in space})==224
start=(zero,5); seen={start};queue=deque([start])
while queue:
 p=queue.popleft()
 for name in local:
  q=act(name,p)
  if q not in seen:seen.add(q);queue.append(q)
assert seen==set(product(P3,[5]))
assert len(seen)==28
print('exact Delta0 orbit under rho(A/B/C):',len(seen))
print('local projections:',len({p[0] for p in seen}),len({p[1] for p in seen}))
# Power stabilizers at the common subgroup (identity actions included).
for name,order in (('A',14),('B',42),('C',42)):
 p=start;good=[]
 for n in range(order):
  if p==start:good.append(n)
  p=act(name,p)
 assert p==start
 assert good==([0,7] if name=='A' else [0,7,14,21,28,35])
 print(name,'common-level power stabilizer:',good)
print('PRINCIPAL DESCENT V57: PASS; ambient index224; curve subgroup index28; CM-order index27')
