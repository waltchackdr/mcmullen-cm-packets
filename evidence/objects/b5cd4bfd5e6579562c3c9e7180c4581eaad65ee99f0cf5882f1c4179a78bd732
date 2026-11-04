#!/usr/bin/env python3
"""
Exact saturation certificate for the three normalized elliptic vertices.

Main checks:
  * B_t and C_t satisfy X^2-tX+1 and generate maximal CM order O_M.
  * For all 42 alpha_k=z^k(1+z+z^4), the matrices obtained by z->B_t and z->C_t
    are integral, have trace c_k and determinant pi, hence are degree-43 Hecke representatives
    fixing the B and C vertices.
  * A_t yields the previously identified four d=7 trace representatives.
  * C-vertex CM type, K3 signature/discriminant, and Hodge-ring structure.
"""
import sympy as sp, itertools, math

t=sp.symbols("t")
f=sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)
def P(x): return sp.Poly(sp.expand(x),t,domain=sp.QQ)
def red(x): return P(x).rem(f).as_expr()
def invL(x): return sp.invert(P(x),f).as_expr()
ti=invL(t)

A=sp.Matrix([
 [t**5-4*t**3-t**2+3*t, -t**5+t**4+6*t**3-3*t**2-10*t],
 [t**3-1, -t**5+5*t**3+t**2-6*t]
]).applyfunc(red)
B=sp.Matrix([[0,-1],[1,t]])
C=sp.Matrix([
 [2*t-4*t**2-4*t**3+t**4+t**5, -5*t-4*t**2+t**3+t**4],
 [-1+5*t+4*t**2-t**3-t**4, -t+4*t**2+4*t**3-t**4-t**5]
]).applyfunc(red)
h=sp.diag(t,1); hi=sp.diag(ti,1)
At=(h*A*hi).applyfunc(red)
Bt=(h*B*hi).applyfunc(red)
Ct=(h*C*hi).applyfunc(red)
for M in (Bt,Ct):
    assert red(M.trace()-t)==0 and red(M.det()-1)==0

# M=L[z]/(z^2-tz+1), pair a+bz
def pmul(x,y):
    a,b=x;c,d=y
    return red(a*c-b*d),red(a*d+b*c+t*b*d)
def zpow(n):
    q=(sp.Integer(1),sp.Integer(0));z=(sp.Integer(0),sp.Integer(1))
    for _ in range(n):q=pmul(q,z)
    return q
alpha0=(2-t**2,t**3-2*t+1)
pi=t**4+t**3-4*t**2-2*t+5
def alpha(k): return pmul(zpow(k),alpha0)
def eval_pair(pair,M):
    a,b=pair
    return (sp.eye(2)*a+M*b).applyfunc(red)

for k in range(42):
    ak=alpha(k)
    ck=red(2*ak[0]+t*ak[1])
    for name,M in (("B",Bt),("C",Ct)):
        T=eval_pair(ak,M)
        for q in T: sp.Poly(q,t,domain=sp.ZZ)
        assert red(T.trace()-ck)==0
        assert red(T.det()-pi)==0
print("ALL 42 trace classes have integral degree-43 representatives at B and C")

# Four A hits.
for k in (2,14,23,35):
    ak=alpha(k)
    # express z^3 basis lambda = A0+B0 z; solve ak=x+y lambda
    lam=zpow(3)
    y=red(ak[1]*invL(lam[1]))
    x=red(ak[0]-y*lam[0])
    T=(sp.eye(2)*x+At*y).applyfunc(red)
    for q in T: sp.Poly(q,t,domain=sp.ZZ)
    assert red(T.det()-pi)==0
print("k=2,14,23,35 also have integral A-vertex representatives")

# C CM type from exact sign intervals.  The six intervals are disjoint, so
# exact rational left endpoints order them as decreasing real embeddings
# t_r=2*cos(r*pi/21), r=1,5,11,13,17,19.
ivs=[iv for iv,m in sp.polys.polytools.intervals(f,eps=sp.Rational(1,10)**22)]
ivs=sorted(ivs,key=lambda ab:sp.Rational(ab[0]),reverse=True)
assert len(ivs)==6
assert all(sp.Rational(ivs[i][0])>sp.Rational(ivs[i+1][1]) for i in range(5))
ks=[1,5,11,13,17,19]
def qint(expr,iv):
    pp=sp.Poly(expr,t,domain=sp.QQ);lo,hi=map(sp.Rational,iv)
    a=b=pp.all_coeffs()[0]
    for c in pp.all_coeffs()[1:]:
        vals=[a*lo,a*hi,b*lo,b*hi]; a=min(vals)+c;b=max(vals)+c
    return a,b
def sign(expr,iv):
    a,b=qint(expr,iv); assert a>0 or b<0
    return 1 if a>0 else -1
cs=[sign(Ct[1,0],iv) for iv in ivs]
PhiC=[k if cs[i]>0 else (-k)%42 for i,k in enumerate(ks)]  # sin(k*pi/21)>0 for all listed k
assert PhiC==[1,5,11,13,17,19]
print("Phi_C =",PhiC)

# K3 operator and Hermitian determinant.
def mmred(X,Y): return (X*Y).applyfunc(red)
def mpow(X,n):
    R=sp.eye(2);Q=X
    while n:
        if n&1:R=mmred(R,Q)
        Q=mmred(Q,Q);n//=2
    return R
SC=(sp.eye(2)+2*mpow(Ct,14)).applyfunc(red)
assert all(red(q)==0 for q in SC*SC+3*sp.eye(2))
# Actual homology uses the transpose. In the K-basis (0,t^j),
# Gram=-Tr(eta*t^(i+j)); even rank six leaves this determinant unchanged.
# The actual 12-dimensional identity and basis are checked by v56.
eta=sp.Poly(SC[1,0],t,domain=sp.QQ)
discL=int(sp.discriminant(f.as_expr(),t))
detH=discL*int(sp.resultant(f.as_expr(),eta.as_expr(),t))
assert detH==-(2**6)*(3**6)*(7**5)*127
assert 127==10**2+3*3**2
print("det H_C =",detH,"=> K3 discriminant class -1")

# Hodge ring.
G=[a for a in range(1,42) if math.gcd(a,42)==1];PH=set(PhiC)
def isH(S,p):
    S=set(S)
    return all(len({g*s%42 for s in S}&PH)==p for g in G)
H={}
for p in range(1,7):
    H[p]={frozenset(S) for S in itertools.combinations(G,2*p) if isH(S,p)}
assert [len(H[p]) for p in range(1,7)]==[6,15,22,15,6,1]
pairs={frozenset({a,(-a)%42}) for a in G}
def divprod(S): return all(frozenset({a,(-a)%42}).issubset(S) for a in S)
exc=[S for S in H[3] if not divprod(S)]
Kplus=frozenset(a for a in G if a%3==1)
Kminus=frozenset(a for a in G if a%3==2)
assert set(exc)=={Kplus,Kminus}
print("C Hodge counts = 6,15,22,15,6,1; only exceptions are the two K3 Weil lines")
print("ALL VERTEX SATURATION ASSERTIONS PASSED")
