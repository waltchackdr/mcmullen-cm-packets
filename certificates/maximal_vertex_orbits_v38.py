#!/usr/bin/env python3
"""Exact certificate for six maximal-order Hecke traces lying in B/C vertex SL2(O_L)-orbits.

For T with characteristic polynomial X^2-cX+pi, the rational matrix P0=[e1,Te1]
satisfies P0 A_c P0^{-1}=T, where A_c=[[0,-pi],[1,c]]. Its determinant
is T_21.  P0 need not be integral.  Nevertheless every other L-conjugator differs
from P0 by an element of the centralizer M^x, whose L-determinant is a relative
norm and is therefore totally positive.  Thus the determinant sign vector is an
invariant of the conjugacy orientation.  For maximal-order traces, the manuscript's
GL2(O_L)-collapse theorem supplies an integral conjugator Q.  If the sign vector
is all positive, det(Q) is a totally positive unit; the exact unit-norm theorem then
adjusts Q through the integral centralizer to determinant one.
"""
import sympy as sp

t=sp.symbols('t')
f=sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)

def red(x):
    return sp.Poly(sp.expand(x),t,domain=sp.QQ).rem(f).as_expr()

def mm(X,Y): return (X*Y).applyfunc(red)
def mpow(X,n):
    R=sp.eye(2); Q=X
    while n:
        if n&1: R=mm(R,Q)
        Q=mm(Q,Q); n//=2
    return R

tinv=sp.invert(sp.Poly(t,t,domain=sp.QQ),f).as_expr()
h=sp.diag(t,1); hi=sp.diag(tinv,1)
B=sp.Matrix([[0,-1],[1,t]])
C=sp.Matrix([[2*t-4*t**2-4*t**3+t**4+t**5,-5*t-4*t**2+t**3+t**4],
             [-1+5*t+4*t**2-t**3-t**4,-t+4*t**2+4*t**3-t**4-t**5]]).applyfunc(red)
Bt=mm(mm(h,B),hi); Ct=mm(mm(h,C),hi)
I=sp.eye(2)

def alpha0(X): return (I+X+mpow(X,4)).applyfunc(red)

def alpha_k(X,k): return mm(mpow(X,k),alpha0(X))

# Exact rational root intervals for the six real embeddings.  Their exact
# left endpoints give the decreasing order because consecutive isolating
# intervals are certified disjoint.
ivs=[ab for ab,m in sp.polys.polytools.intervals(f,eps=sp.Rational(1,10)**30)]
ivs=sorted(ivs,key=lambda ab:sp.Rational(ab[0]),reverse=True)
assert len(ivs)==6
assert all(sp.Rational(ivs[i][0])>sp.Rational(ivs[i+1][1]) for i in range(5))

def interval_eval(poly,ab):
    p=sp.Poly(poly,t,domain=sp.QQ)
    lo,hi=map(sp.Rational,ab)
    a=b=p.all_coeffs()[0]
    for c in p.all_coeffs()[1:]:
        vals=(a*lo,a*hi,b*lo,b*hi)
        a=min(vals)+c; b=max(vals)+c
    return a,b

def signs(poly):
    out=[]
    for ab in ivs:
        lo,hi=interval_eval(poly,ab)
        assert lo>0 or hi<0, (poly,ab,lo,hi)
        out.append(1 if lo>0 else -1)
    return tuple(out)

def norm(poly):
    return sp.Integer(sp.resultant(f.as_expr(),sp.Poly(poly,t,domain=sp.QQ).as_expr(),t))

# Chosen holomorphic representatives. True means replace T by cI-T.
cases=[('B',Bt,11,False),('B',Bt,32,True),
       ('C',Ct,16,False),('C',Ct,19,False),
       ('C',Ct,37,True),('C',Ct,40,True)]

for name,X,k,use_conj in cases:
    T=alpha_k(X,k)
    c=red(sp.trace(T)); pi=red(T.det())
    if use_conj:
        T=(I*c-T).applyfunc(red)
    # T still has trace c and determinant pi.
    assert red(sp.trace(T)-c)==0
    assert red(T.det()-pi)==0
    P=sp.Matrix([[1,T[0,0]],[0,T[1,0]]]).applyfunc(red)
    d=red(P.det())
    assert signs(d)==(1,1,1,1,1,1), (name,k,signs(d))
    assert norm(d)!=0
    Ac=sp.Matrix([[0,-pi],[1,c]])
    # Verify P A_c = T P, avoiding denominators.
    assert all(red(x)==0 for x in (P*Ac-T*P))
    print(f'k={k:2d} vertex={name} conjugate={use_conj}: det(P0) totally positive; P0 A_c=T P0 PASS; Norm(det P0)={norm(d)}')

print('RESOLVED MAXIMAL-ORDER INDICES: 11,16,19,32,37,40')
print('REMAINING MAXIMAL-ORDER d=3 INDICES: 1,17,20,22,38,41')
print('ALL EXACT ASSERTIONS PASSED')
