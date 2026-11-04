#!/usr/bin/env python3
"""Restricted free-line incompatibility in Mobius coordinates at A,B,C.

Transposed homology and common nonfree H are verified in principal_homology_v56.py.
This script does not exclude general maximal isotropic subgroups."""
# This script is a focused exact check extracted from the vertex certificates.
import sympy as sp, itertools

t,u=sp.symbols("t u")
f=sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)
def P(x): return sp.Poly(sp.expand(x),t,domain=sp.QQ)
def red(x): return P(x).rem(f).as_expr()
ti=sp.invert(P(t),f).as_expr()
A=sp.Matrix([[t**5-4*t**3-t**2+3*t,-t**5+t**4+6*t**3-3*t**2-10*t],
             [t**3-1,-t**5+5*t**3+t**2-6*t]]).applyfunc(red)
B=sp.Matrix([[0,-1],[1,t]])
C=sp.Matrix([[2*t-4*t**2-4*t**3+t**4+t**5,-5*t-4*t**2+t**3+t**4],
             [-1+5*t+4*t**2-t**3-t**4,-t+4*t**2+4*t**3-t**4-t**5]]).applyfunc(red)
h=sp.diag(t,1);hi=sp.diag(ti,1)
At=(h*A*hi).applyfunc(red);Bt=(h*B*hi).applyfunc(red);Ct=(h*C*hi).applyfunc(red)
def mm(X,Y): return (X*Y).applyfunc(red)
def mpow(X,n):
    R=sp.eye(2);Q=X
    while n:
        if n&1:R=mm(R,Q)
        Q=mm(Q,Q);n//=2
    return R
def pmul(x,y):
    a,b=x;c,d=y
    return red(a*c-b*d),red(a*d+b*c+t*b*d)
def zpow(n):
    q=(sp.Integer(1),sp.Integer(0));z=(sp.Integer(0),sp.Integer(1))
    for _ in range(n):q=pmul(q,z)
    return q
lam=zpow(3);z14=zpow(14);s=(red(1+2*z14[0]),red(2*z14[1]))
y=red(s[1]*sp.invert(P(lam[1]),f).as_expr());x=red(s[0]-y*lam[0])
SA=(sp.eye(2)*x+At*y).applyfunc(red)
SB=(sp.eye(2)+2*mpow(Bt,14)).applyfunc(red)
SC=(sp.eye(2)+2*mpow(Ct,14)).applyfunc(red)
Ss={"A":SA,"B":SB,"C":SC}

# Work modulo p7^2 = (7,(t-1)^2).
N=2
def Ru(expr):
    e=sp.expand(P(expr).as_expr().subs(t,1+u))
    p=sp.Poly(e,u,modulus=7)
    out=[0]*N
    for (deg,),c in p.terms():
        if deg<N: out[deg]=int(c)%7
    return tuple(out)
def add(a,b):return tuple((a[i]+b[i])%7 for i in range(N))
def sub(a,b):return tuple((a[i]-b[i])%7 for i in range(N))
def mul(a,b):
    z=[0]*N
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<N:z[i+j]=(z[i+j]+x*y)%7
    return tuple(z)
def inv_eq(M,m):
    a,b=M[0];c,d=M[1]
    return sub(add(c,mul(d,m)),mul(m,add(a,mul(b,m))))

local={name:[[Ru(M[i,j]) for j in range(2)] for i in range(2)] for name,M in Ss.items()}
lift={}
for name,M in local.items():
    vals=[]
    for a in range(7):
        m=(4,a)
        if inv_eq(M,m)==(0,0):vals.append(a)
    assert len(vals)==1
    lift[name]=vals[0]
assert lift=={"A":4,"B":1,"C":5}
print("unique invariant lifts of common residue line (1,4):")
for name,a in lift.items():
    print(f"  {name}: (1, 4 + {a}(t-1)) mod p7^2")
assert len(set(lift.values()))==3
print("No pair has a common invariant FREE rank-one summand; nonfree levels are not excluded.")
print("THREE-VERTEX RESTRICTED FREE-LINE ASSERTIONS PASSED")
