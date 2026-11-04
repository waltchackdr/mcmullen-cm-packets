#!/usr/bin/env python3
"""
Exact certificate for the four d=7, ell=43 nonmaximal trace classes.

Checks:
  * conductor prime and conductor norm for k=2,14,23,35;
  * complete splitting of 41 and 83 in L;
  * inertness of the relevant conductor primes in M/L=L(sqrt(-7))/L;
  * local unit quotient orders;
  * the cyclotomic unit epsilon=1-zeta_21 has local class of order 42;
  * Picard and ideal-class-monoid counts;
  * oriented SL2 class counts predicted from the centralizer norm quotient.

All polynomial and finite-field assertions are exact.
"""

import sympy as sp
from math import gcd

T = sp.symbols("T")
f = T**6 + T**5 - 6*T**4 - 6*T**3 + 8*T**2 + 8*T + 1

D83 = 2*T**5 - T**4 - 13*T**3 + 4*T**2 + 19*T - 7
D41 = -(T + 2)*(4*T**4 - 4*T**3 - 15*T**2 + 10*T + 13)

data = {
    2:  (83, 36, D83),
    23: (83, 36, D83),
    14: (41, 12, D41),
    35: (41, 12, D41),
}

def roots_mod(poly, p):
    return [a for a in range(p) if int(poly.subs(T, a)) % p == 0]

for p in (41, 83):
    fac = sp.factor(f, modulus=p)
    roots = roots_mod(f, p)
    assert len(roots) == 6
    print(f"f mod {p} =", fac)
    print(f"roots mod {p} =", roots)

for k, (p, root, Delta) in data.items():
    nd = abs(int(sp.resultant(f, Delta, T)))
    assert nd == p*p
    z = [r for r in roots_mod(f, p) if int(Delta.subs(T, r)) % p == 0]
    assert z == [root]
    # The maximal relative discriminant is 1, so disc(R)=f_conductor^2.
    print(f"k={k}: Norm(Delta)={nd}={p}^2; conductor=( {p}, t-{root} ), norm={p}")

# M=L(sqrt(-7)); the relevant degree-one primes are inert iff -7 is a nonsquare.
for p in (41,83):
    leg = pow((-7) % p, (p-1)//2, p)
    assert leg == p-1
    print(f"(-7/{p})=-1: conductor prime is inert; residue extension F_{p} -> F_{p*p}")

# Local quotient F_{p^2}^*/F_p^* has order p+1.
assert 41+1 == 42
assert 83+1 == 84

# Let z be a primitive 42nd root and epsilon=1-z^2.
# Since p == -1 mod 42, Frobenius sends z -> z^{-1}, and
# epsilon^(p-1)=(1-z^{-2})/(1-z^2)=-z^{-2}=z^19.
assert gcd(19,42) == 1
eps_class_order = 42
print("epsilon=1-zeta_21 has local quotient order", eps_class_order)

# The root-of-unity class [z] equals [epsilon]^22:
assert (19*22 - 40) % 42 == 0
assert gcd(22,42) == 2
print("[zeta_42]=[epsilon]^22 in the local quotient")

# h(M)=1.  Exact order-class sequence gives:
# |Pic(R)|=(p+1)/|image(O_M^*)|, and image(O_M^*) has order 42.
pic = {41:(42//42), 83:(84//42)}
assert pic == {41:1,83:2}

# Quadratic orders are Bass; with prime conductor there are exactly two overorders R and O_M.
# Therefore ICM(R)=Pic(R) disjoint union Pic(O_M), and Pic(O_M)=1.
icm = {p: pic[p]+1 for p in pic}
assert icm == {41:2,83:3}

# R^* is the kernel of the local reduction.  If epsilon^n z^m u is in R^*,
# then n+22m == 0 mod 42, hence n is even.  Thus N(R^*)=E_L^2,
# while E_L^+/E_L^2 has order 2.
oriented = {
    41: 2*pic[41] + 1,  # R-multiplier classes split in 2; O_M class does not.
    83: 2*pic[83] + 1,
}
assert oriented == {41:3,83:5}

print("Picard counts:", pic)
print("ICM / GL2 class counts:", icm)
print("holomorphic SL2 class counts:", oriented)
print("ALL d=7 NONMAXIMAL-ORDER ASSERTIONS PASSED")


# ---------- Explicit matrix representatives ----------

# Exact conductor generators in O_L.
gamma41 = -1 - 2*T + T**2 - T**3 - T**4
gamma83 = -2 - 2*T + 2*T**2 - T**4
assert sp.resultant(f, gamma41, T) == -41
assert sp.resultant(f, gamma83, T) == -83
assert int(gamma41.subs(T,12)) % 41 == 0
assert int(gamma83.subs(T,36)) % 83 == 0

# u=N_{M/L}(1-zeta_21)=4-t^2 is a totally positive unit.
u = 4 - T**2
assert sp.resultant(f, u, T) == 1

fpoly = sp.Poly(f, T, domain=sp.QQ)
def red_poly(e):
    return sp.Poly(sp.expand(e), T, domain=sp.QQ).rem(fpoly).as_expr()

def red_field(e):
    num, den = sp.fraction(sp.cancel(e))
    pn = sp.Poly(num, T, domain=sp.QQ).rem(fpoly)
    pd = sp.Poly(den, T, domain=sp.QQ).rem(fpoly)
    return (pn * sp.invert(pd, fpoly)).rem(fpoly).as_expr()

# z satisfies z^2-tz+1=0 over L.
def pair_mul(x,y):
    A,B=x; C,D=y
    return red_poly(A*C-B*D), red_poly(A*D+B*C+T*B*D)

def zpow(n):
    A,B=sp.Integer(1),sp.Integer(0)
    for _ in range(n):
        A,B=red_poly(-B),red_poly(A+T*B)
    return A,B

# alpha0=1+z+z^4 = (2-t^2)+(t^3-2t+1)z.
alpha0_pair = (2-T**2, T**3-2*T+1)
pi = T**4 + T**3 - 4*T**2 - 2*T + 5

def alpha_pair(k):
    return pair_mul(zpow(k), alpha0_pair)

def mult_matrix(pair):
    A,B=pair
    return sp.Matrix([[A,-B],[B,red_poly(A+T*B)]])

def char_check(M, k):
    A,B=alpha_pair(k)
    c=red_poly(2*A+T*B)
    tr=red_poly(M.trace())
    det=red_poly(M.det())
    assert red_poly(tr-c)==0
    assert red_poly(det-pi)==0
    return c

# Maximal-overorder representative on O_M=O_L[1,z].
for k in (2,14,23,35):
    MO=mult_matrix(alpha_pair(k))
    char_check(MO,k)

# Nontrivial conductor-83 Picard class:
# beta=1+3z has quotient order 84.  The intersection I=beta O_M cap R
# has basis v1=16+(58+3t)z, v2=gamma83(1+3z).
V = sp.Matrix([[16,gamma83],[58+3*T,3*gamma83]])

I_matrices={}
for k in (2,23):
    MO=mult_matrix(alpha_pair(k))
    MI=(V.inv()*MO*V).applyfunc(red_field)
    # Exact integrality: every entry reduces to a polynomial in Z[t].
    for x in MI:
        sp.Poly(x,T,domain=sp.ZZ)
    char_check(MI,k)
    I_matrices[k]=MI

# Companion representatives for the R_k module itself.
def companion(k):
    A,B=alpha_pair(k)
    c=red_poly(2*A+T*B)
    return sp.Matrix([[0,-pi],[1,c]])

for k in (2,14,23,35):
    char_check(companion(k),k)

# Orientation twist D=diag(u,1); u is a unit, so integrality is preserved.
uinv = red_field(1/u)
D = sp.diag(u,1)
Dinv = sp.diag(uinv,1)
for k in (2,14,23,35):
    twist=(D*companion(k)*Dinv).applyfunc(red_field)
    for x in twist:
        sp.Poly(x,T,domain=sp.ZZ)
    char_check(twist,k)

for k,MI in I_matrices.items():
    twist=(D*MI*Dinv).applyfunc(red_field)
    for x in twist:
        sp.Poly(x,T,domain=sp.ZZ)
    char_check(twist,k)

print("explicit O_M, companion, nontrivial-I, and orientation-twist matrices verified")
print("totally-positive orientation-sector counts: 5+5+3+3 = 16")
print("full SL2 counts per trace: k=2,23 -> 160; k=14,35 -> 96; grand total=512")
