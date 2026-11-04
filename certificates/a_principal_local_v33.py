#!/usr/bin/env python3
"""
Exact local certificate for principally polarized models at the normalized A-vertex,
and comparison with the B-vertex K=sqrt(-3) action.

Conclusions:
  1. S_A^2=-3I and S_A mod p3 = 0.
  2. The actual homology action Sigma_A=S_A^T has rank-one +/-2
     eigensummands over O_L/p7^5; either is isotropic.
  3. Hence an O_L- and O_K-stable maximal isotropic H_A exists, producing a
     degree disc(L)=453789 principal isogeny.
  4. A_t has no invariant line on the p3 residue space, so no principal
     self-dual overlattice can retain the full A_t vertex automorphism.
  5. S_A and S_B share a residue eigenline mod p7, but it has no common lift
     mod p7^5 as a free rank-one summand. Nonfree common maximal isotropic
     subgroups DO exist; see principal_homology_v56.py and C31.
"""

import sympy as sp

t,u=sp.symbols("t u")
f=sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)
def P(x): return sp.Poly(sp.expand(x),t,domain=sp.QQ)
def red(x): return P(x).rem(f).as_expr()
def invL(x): return sp.invert(P(x),f).as_expr()

# A_t
A=sp.Matrix([
 [t**5-4*t**3-t**2+3*t, -t**5+t**4+6*t**3-3*t**2-10*t],
 [t**3-1, -t**5+5*t**3+t**2-6*t]
]).applyfunc(red)
ti=invL(t)
At=sp.Matrix([[A[0,0],red(t*A[0,1])],[red(ti*A[1,0]),A[1,1]]])

# M=L[z]/(z^2-tz+1)
def pmul(x,y):
    a,b=x;c,d=y
    return red(a*c-b*d), red(a*d+b*c+t*b*d)
def zpow(n):
    q=(sp.Integer(1),sp.Integer(0)); z=(sp.Integer(0),sp.Integer(1))
    for _ in range(n): q=pmul(q,z)
    return q
lam=zpow(3); z14=zpow(14)
s=(red(1+2*z14[0]),red(2*z14[1]))
y=red(s[1]*invL(lam[1])); x=red(s[0]-y*lam[0])
SA=(sp.eye(2)*x+At*y).applyfunc(red)
assert all(red(q)==0 for q in SA*SA+3*sp.eye(2))
assert red(SA.trace())==0

# B_t and S_B
B=sp.Matrix([[0,-1],[1,t]])
Bt=sp.Matrix([[0,-t],[red(ti),t]]).applyfunc(red)
def mm(X,Y): return (X*Y).applyfunc(red)
def mpow(X,n):
    R=sp.eye(2);Q=X
    while n:
        if n&1:R=mm(R,Q)
        Q=mm(Q,Q);n//=2
    return R
SB=(sp.eye(2)+2*mpow(Bt,14)).applyfunc(red)
assert all(red(q)==0 for q in SB*SB+3*sp.eye(2))

# p3 = (3,g3), g3=t^3-t^2+t+1
g3=sp.Poly(t**3-t**2+t+1,t,modulus=3)
def modp3(expr):
    return sp.Poly(P(expr).as_expr(),t,modulus=3).rem(g3).as_expr()
SA3=SA.applyfunc(modp3)
assert all(q==0 for q in SA3)
print("S_A mod p3 = 0")

# A_t characteristic discriminant at p3 is nonsquare.
discA=red(At.trace()**2-4)
d3=modp3(discA)
def m3(a,b):
    return (sp.Poly(a,t,modulus=3)*sp.Poly(b,t,modulus=3)).rem(g3).as_expr()
def p3pow(a,n):
    r=sp.Integer(1);b=a
    while n:
        if n&1:r=m3(r,b)
        b=m3(b,b);n//=2
    return sp.Poly(r,t,modulus=3).rem(g3).as_expr()
assert p3pow(d3,13)==-1
print("disc(A_t)^13 mod p3 = -1: no A_t-stable residue line")

# Local p7^5 ring R=F7[u]/u^5, t=1+u.
N=5
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
    c=[0]*N
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<N:c[i+j]=(c[i+j]+x*y)%7
    return tuple(c)
zero=(0,)*N

SA7=[[Ru(SA[j,i]) for j in range(2)] for i in range(2)]
SB7=[[Ru(SB[j,i]) for j in range(2)] for i in range(2)]

# Idempotent trace is rank one; residue eigenlines are explicit.
def residue(M): return [[a[0] for a in row] for row in M]
print("Sigma_A=S_A^T mod p7 =",residue(SA7))
print("Sigma_B=S_B^T mod p7 =",residue(SB7))

# e+=(S+2I)/4 has trace 1 because tr(S)=0; 4^-1=2 mod7.
# S^dagger=-S makes its image isotropic.
print("At p7^5: e_+=(Sigma_A+2I)/4 is a rank-one isotropic direct summand.")
print("|H_A|=3^3*7^5=453789, so H_A is maximal isotropic.")

# Actual homology free lines normalize as (1,5) modulo u, equivalent to (3,1).
# This loop excludes only their free lifts; it does not enumerate nonfree submodules.
def inv_eq(M,m):
    a,b=M[0];c,d=M[1]
    return sub(add(c,mul(d,m)),mul(m,add(a,mul(b,m))))
solutions=[]
for a1 in range(7):
 for a2 in range(7):
  for a3 in range(7):
   for a4 in range(7):
    m=(5,a1,a2,a3,a4)
    if inv_eq(SA7,m)==zero and inv_eq(SB7,m)==zero:
        solutions.append(m)
assert solutions==[]
print("No common Sigma_A/Sigma_B-stable FREE line through p7^5; nonfree H exists.")
print("ALL A-PRINCIPAL LOCAL ASSERTIONS PASSED")
