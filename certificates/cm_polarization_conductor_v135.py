"""Exact bounded conductor-generator witness using only rational arithmetic.

Finite input certificate only; the classification remains a separate human proof.
Coefficient vectors are in increasing powers of zeta_42.
"""
from fractions import Fraction as Q
import json
import resource
import signal

resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
signal.alarm(30)
if not __debug__:
    raise RuntimeError("Enabled assertions are required for this certificate")
phi = [1,1,0,-1,-1,0,1,0,-1,-1,0,1,1]
one = [Q(1)] + [Q(0)]*11
z = [Q(0),Q(1)] + [Q(0)]*10

def add(a,b): return [x+y for x,y in zip(a,b)]
def scale(c,a): return [c*x for x in a]
def mul(a,b):
    p = [Q(0)]*23
    for i,x in enumerate(a):
        for j,y in enumerate(b): p[i+j] += x*y
    for k in range(22,11,-1):
        for j in range(12): p[k-12+j] -= p[k]*phi[j]
    return p[:12]

def power(a,n):
    b=one
    while n:
        if n%2: b=mul(b,a)
        a=mul(a,a);n//=2
    return b

def det(a):
    a=[r[:] for r in a];v=Q(1)
    for i in range(len(a)):
        j=next(j for j in range(i,len(a)) if a[j][i])
        if i!=j:a[i],a[j]=a[j],a[i];v=-v
        p=a[i][i];v*=p
        for j in range(i+1,len(a)):
            q=a[j][i]/p
            for k in range(i,len(a)):a[j][k]-=q*a[i][k]
    return v

assert power(z,42)==one
t=add(z,power(z,41))
beta=[Q(0)]*12
for i,c in enumerate([-2,-8,-2,6,1,-1]): beta=add(beta,scale(c,power(t,i)))
q=add(one,scale(-1,power(z,14)))
qbar=add(one,scale(-1,power(z,28)))
assert mul(q,qbar)==scale(3,one)
a=scale(Q(1,3),mul(beta,qbar))
assert mul(a,q)==beta and all(x.denominator==1 for x in a)
columns=[mul(a,power(z,j)) for j in range(12)]
norm=det([list(r) for r in zip(*columns)])
assert abs(norm)==1
print(json.dumps({'status':'CM_POLARIZATION_CONDUCTOR_INPUT_PASS',
    'coefficient_order':'1,z,...,z^11; z=zeta42',
    'phi42':phi,'beta':[int(x) for x in beta],
    'q':[int(x) for x in q], 'beta_over_q':[int(x) for x in a],
    'norm_beta_over_q':int(norm), 'identity_beta_equals_q_times_unit':True,
    'conclusion':'beta O = (1-z^14) O'},indent=2))
