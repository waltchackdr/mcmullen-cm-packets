#!/usr/bin/env python3
# Exact certificate for the normalized A-vertex.
import sympy as sp
import itertools, math

t = sp.symbols("t")
f = sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)
def red(x):
    return sp.Poly(sp.cancel(x),t,domain=sp.QQ).rem(f).as_expr()
def fred(x):
    x=sp.cancel(x)
    n,d=sp.fraction(x)
    pn=sp.Poly(n,t,domain=sp.QQ).rem(f)
    pd=sp.Poly(d,t,domain=sp.QQ).rem(f)
    return (pn*sp.invert(pd,f)).rem(f).as_expr()

# McMullen's quadrilateral generator A, with A^7=-I.
A = sp.Matrix([
 [t**5-4*t**3-t**2+3*t, -t**5+t**4+6*t**3-3*t**2-10*t],
 [t**3-1, -t**5+5*t**3+t**2-6*t]
]).applyfunc(red)
assert red(A.det()-1)==0
assert red(A.trace()-(t**3-3*t))==0

# h_t=diag(t,1), t is a unit.
tinverse = sp.invert(sp.Poly(t,t,domain=sp.QQ),f).as_expr()
At=sp.Matrix([[A[0,0],red(t*A[0,1])],[red(tinverse*A[1,0]),A[1,1]]])
assert red(At.det()-1)==0

# Exact real-root intervals and lower-left signs.  The six intervals are
# disjoint, so their exact rational left endpoints determine their order.
# This is the decreasing order t_r=2*cos(r*pi/21), r=1,5,11,13,17,19.
ivs=[x for x,m in sp.polys.polytools.intervals(f,eps=sp.Rational(1,10)**24)]
ivs=sorted(ivs,key=lambda ab:sp.Rational(ab[0]),reverse=True)
assert len(ivs)==6
assert all(sp.Rational(ivs[i][0])>sp.Rational(ivs[i+1][1]) for i in range(5))
def qinterval(poly,ab):
    p=sp.Poly(poly,t,domain=sp.QQ)
    lo,hi=map(sp.Rational,ab)
    a=b=p.all_coeffs()[0]
    for c in p.all_coeffs()[1:]:
        vals=[a*lo,a*hi,b*lo,b*hi]
        a=min(vals)+c; b=max(vals)+c
    return a,b
def sign(poly,ab):
    a,b=qinterval(poly,ab)
    assert a>0 or b<0
    return 1 if a>0 else -1
ks=[1,5,11,13,17,19]
cs=[sign(At[1,0],iv) for iv in ivs]
sin_sign={1:1,5:1,11:-1,13:-1,17:1,19:1}
Phi=[k if cs[i]==sin_sign[k] else (-k)%42 for i,k in enumerate(ks)]
assert Phi==[1,5,31,29,17,19]
print("Phi_A =",Phi)

res7={1,2,4}
sig7=sum((a%7) in res7 for a in Phi)
sig3=sum(a%3==1 for a in Phi)
assert (sig7,6-sig7)==(2,4)
assert (sig3,6-sig3)==(3,3)
print("K=-7 signature = (2,4); K=-3 signature = (3,3)")

# Polynomial arithmetic in M=L[z]/(z^2-tz+1), represented A+Bz.
def pmul(x,y):
    a,b=x; c,d=y
    return red(a*c-b*d), red(a*d+b*c+t*b*d)
def zpow(n):
    q=(sp.Integer(1),sp.Integer(0))
    z=(sp.Integer(0),sp.Integer(1))
    for _ in range(n): q=pmul(q,z)
    return q

# sqrt(-3)=1+2 z^14 = x+y z^3.
lam=zpow(3); s14=zpow(14)
s3=(red(1+2*s14[0]),red(2*s14[1]))
y=fred(s3[1]/lam[1]); x=red(s3[0]-y*lam[0])
SA=(sp.eye(2)*x+At*y).applyfunc(red)
assert all(red(q)==0 for q in SA*SA+3*sp.eye(2))
for q in SA: sp.Poly(q,t,domain=sp.ZZ)

# Actual homology uses the transpose. In the K-basis (0,t^j),
# Gram=-Tr(eta*t^(i+j)); even rank six leaves this determinant unchanged.
# The actual 12-dimensional identity and basis are checked by v56.
eta=sp.Poly(SA[1,0],t,domain=sp.QQ)
discL=int(sp.discriminant(f.as_expr(),t))
Neta=int(sp.resultant(f.as_expr(),eta.as_expr(),t))
detH=discL*Neta
assert detH == -(2**6)*(3**6)*(7**5)*127
assert 127 == 10**2+3*3**2
print("det(H_A) =",detH,"=> discriminant class -1")

# Hodge character enumeration and full divisor-algebra closure.
G=[a for a in range(1,42) if math.gcd(a,42)==1]
P=set(Phi)
def hodge(S,p):
    S=set(S)
    return all(len({g*s%42 for s in S}&P)==p for g in G)
H={}
for p in range(1,7):
    H[p]={frozenset(S) for S in itertools.combinations(G,2*p) if hodge(S,p)}
assert [len(H[p]) for p in range(1,7)] == [12,51,88,51,12,1]
D=H[1]
generated=[]
for p in range(1,7):
    U=set()
    for comb in itertools.combinations(D,p):
        S=frozenset().union(*comb)
        if len(S)==2*p: U.add(S)
    assert U==H[p]
    generated.append(len(U))
print("Hodge counts =",[len(H[p]) for p in range(1,7)])
print("divisor-product counts =",generated)

stab=[g for g in G if {g*a%42 for a in P}==P]
assert stab==[1,29]
print("stabilizer(Phi_A) =",stab,"; fixed CM field is Q(zeta_7)")
print("ALL A-VERTEX ASSERTIONS PASSED")
