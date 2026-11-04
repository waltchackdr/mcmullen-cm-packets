#!/usr/bin/env python3
"""Exact certificate for Theorem packetHodge in hecke_correction_note_v39.tex."""
import sympy as sp, itertools, math
z,t=sp.symbols('z t')
Phi=sp.Poly(sp.cyclotomic_poly(42,z),z,domain=sp.QQ)
f=sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)
def redz(expr):
    n,d=sp.fraction(sp.cancel(expr))
    return (sp.Poly(n,z,domain=sp.QQ)*sp.invert(sp.Poly(d,z,domain=sp.QQ),Phi)).rem(Phi)
def conjz(p): return redz(p.as_expr().subs(z,z**41))
tz=redz(z+z**41)
bas=[redz(tz.as_expr()**j) for j in range(6)]
B=sp.Matrix([[bas[j].nth(i) for j in range(6)] for i in range(12)])
def to_t(p):
    v=sp.Matrix([p.nth(i) for i in range(12)])
    s=next(iter(sp.linsolve((B,v))))
    return sp.Poly(sum(s[j]*t**j for j in range(6)),t,domain=sp.QQ).rem(f)
G=[a for a in range(1,42) if math.gcd(a,42)==1]
base=[1,5,11,13,17,19]
a0=sp.Poly(1+z+z**4,z,domain=sp.QQ)
sqrt3=redz(1+2*z**14)
discL=int(sp.discriminant(f.as_expr(),t))
indices=[1,11,16,17,19,20,22,32,37,38,40,41]

def redt(expr):
    return sp.Poly(sp.expand(expr),t,domain=sp.QQ).rem(f).as_expr()

# q_n(t)=(z^n-z^-n)/(z-z^-1), reconstructed in L by the exact recurrence
# q_0=0, q_1=1, q_{n+1}=t*q_n-q_{n-1}.
q=[sp.Integer(0),sp.Integer(1)]
for _ in range(1,46):
    q.append(redt(t*q[-1]-q[-2]))

# Exact rational isolating intervals, in the embedding order
# t_r=2*cos(r*pi/21), r=1,5,11,13,17,19.  Sorting rational endpoints is
# exact; no binary float is used to label the embeddings.
root_intervals=[ab for ab,_ in sp.polys.polytools.intervals(
    f,eps=sp.Rational(1,10)**30)]
root_intervals=sorted(root_intervals,key=lambda ab:ab[0],reverse=True)
assert len(root_intervals)==6 and all(0<r<21 for r in base)

def interval_eval(poly,ab):
    p=sp.Poly(poly,t,domain=sp.QQ)
    lo,hi=map(sp.Rational,ab)
    lower=upper=p.all_coeffs()[0]
    for coefficient in p.all_coeffs()[1:]:
        products=(lower*lo,lower*hi,upper*lo,upper*hi)
        lower=min(products)+coefficient
        upper=max(products)+coefficient
    return lower,upper

def cm_type(k):
    # Im(sigma_r(alpha_k)) = sin(r*pi/21)*sigma_r(y_k), and the sine
    # factor is strictly positive for every r in base.  The signs below are
    # therefore certified entirely by rational interval arithmetic in L.
    yk=redt(q[k]+q[k+1]+q[k+4])
    out=[]
    for r,ab in zip(base,root_intervals):
        lo,hi=interval_eval(yk,ab)
        assert lo>0 or hi<0,(k,r,lo,hi)
        out.append(r if lo>0 else (-r)%42)
    return frozenset(out)

def hodge_sets(PH,p):
    return {frozenset(S) for S in itertools.combinations(G,2*p)
            if all(len({g*s%42 for s in S}&PH)==p for g in G)}

types={}
for k in indices:
    ak=(sp.Poly(z**k,z,domain=sp.QQ)*a0).rem(Phi)
    types[k]=cm_type(k)

# Exact Galois-orbit partition of the 12 types.
remaining=set(indices); orbits=[]
while remaining:
    k=min(remaining); S=types[k]
    orb={frozenset((g*a)%42 for a in S) for g in G}
    block=sorted(j for j in indices if types[j] in orb)
    orbits.append(block); remaining-=set(block)
assert orbits==[[1,20,22,41],[11,16,17,19,32,37,38,40]]
print('CM-type orbits:',orbits)

# Stabilizers and Hodge profiles.
for block in orbits:
    S=types[block[0]]
    stab=[g for g in G if frozenset((g*a)%42 for a in S)==S]
    counts=[len(hodge_sets(S,p)) for p in range(1,7)]
    if block[0]==1:
        assert stab==[1,29]
        assert counts==[12,51,88,51,12,1]
    else:
        assert stab==[1]
        assert counts==[6,15,22,15,6,1]
        # Standard opposite-pair divisor products; precisely two middle lines remain.
        pairs={frozenset({a,(-a)%42}) for a in G}
        div3={frozenset().union(*cc) for cc in itertools.combinations(pairs,3)}
        assert len(hodge_sets(S,3)-div3)==2
    print('block',block,'stabilizer',stab,'Hodge counts',counts)

# Actual companion homology is (x_k I+y_k A_c)^T. In the K-basis
# (0,t^j), Gram=-Tr(y_k*t^(i+j)); since rank is six its determinant
# equals disc(L)*N(y_k). See the v56 human proof for the convention.
for k in indices:
    ak=(sp.Poly(z**k,z,domain=sp.QQ)*a0).rem(Phi)
    den=ak-conjz(ak)
    y_z=redz((sqrt3-conjz(sqrt3)).as_expr()*sp.invert(den,Phi).as_expr())
    y_t=to_t(y_z)
    Ny=int(sp.resultant(f.as_expr(),y_t.as_expr(),t))
    detH=discL*Ny
    assert Ny==-1728
    assert detH==-(2**6)*(3**6)*(7**5)
print('All 12: N_L/Q(y_k)=-1728 and det(H_k)=-2^6*3^6*7^5, hence discriminant class -1.')
print('PACKET HODGE v39: PASS')
