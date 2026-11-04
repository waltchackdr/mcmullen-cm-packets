#!/usr/bin/env python3
"""Exact v56 special-fiber homology and common nonfree level certificate.

Human completeness and trace-pairing bridge: docs/PRINCIPAL_HOMOLOGY_REPAIR_V56.md.
No global family, curve monodromy, or exact surviving order is asserted.
This SymPy path develops the September audit witness; the separate standard-
library implementation reconstructs the cyclotomic algebra independently.
"""
import sympy as sp

t, u = sp.symbols("t u")
f = sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1, t)

def red(x):
    return sp.rem(sp.Poly(sp.expand(x), t, domain=sp.QQ), f).as_expr()

def matred(m):
    return m.applyfunc(red)

def mpow(m, n):
    r = sp.eye(2)
    for _ in range(n):
        r = matred(r*m)
    return r

def zpower(n):
    a, b = sp.Integer(1), sp.Integer(0)
    for _ in range(n):
        a, b = red(-b), red(a+t*b)
    return a, b

a = sp.Matrix([[t**5-4*t**3-t**2+3*t,
                -t**5+t**4+6*t**3-3*t**2-10*t],
               [t**3-1, -t**5+5*t**3+t**2-6*t]])
b = sp.Matrix([[0,-1],[1,t]])
c = sp.Matrix([[2*t-4*t**2-4*t**3+t**4+t**5,
                -5*t-4*t**2+t**3+t**4],
               [-1+5*t+4*t**2-t**3-t**4,
                -t+4*t**2+4*t**3-t**4-t**5]])
normalizer = sp.diag(t,1)
normalizer_inv = sp.diag(sp.invert(t,f.as_expr()),1)
a,b,c = [matred(normalizer*m*normalizer_inv) for m in (a,b,c)]
z3a,z3b = zpower(3)
z14a,z14b = zpower(14)
coef = red(2*z14b*sp.invert(z3b,f.as_expr()))
const = red(1+2*z14a-coef*z3a)
operators = {
    "A": matred(const*sp.eye(2)+coef*a),
    "B": matred(sp.eye(2)+2*mpow(b,14)),
    "C": matred(sp.eye(2)+2*mpow(c,14)),
}

def local(x):
    p = sp.Poly(sp.expand(x.subs(t,1+u)),u,domain=sp.QQ)
    return tuple(int(p.nth(i).p)*pow(int(p.nth(i).q),-1,7)%7 for i in range(5))

def add(x,y):
    return tuple((a+b)%7 for a,b in zip(x,y))

def mul(x,y):
    return tuple(sum(x[j]*y[i-j] for j in range(i+1))%7 for i in range(5))

def sub(x,y):
    return tuple((a-b)%7 for a,b in zip(x,y))

def rank(rows):
    rows = [list(r) for r in rows]
    pivot = 0
    for col in range(10):
        hit = next((i for i in range(pivot,len(rows)) if rows[i][col]%7),None)
        if hit is None:
            continue
        rows[pivot],rows[hit] = rows[hit],rows[pivot]
        inv = pow(rows[pivot][col],-1,7)
        rows[pivot] = [(x*inv)%7 for x in rows[pivot]]
        for i in range(len(rows)):
            if i != pivot:
                q = rows[i][col]
                rows[i] = [(x-q*y)%7 for x,y in zip(rows[i],rows[pivot])]
        pivot += 1
    return pivot

def shifted(v,n):
    return tuple(tuple(v[j] if i==n else 0 for i in range(5)) for j in range(2))

def action(m,v):
    return tuple(add(mul(m[i][0],v[0]),mul(m[i][1],v[1])) for i in range(2))

def flat(v):
    return v[0]+v[1]

expected = {"A":[[1,1],[3,6]],"B":[[6,5],[2,1]],"C":[[5,0],[5,2]]}
for name,s in operators.items():
    assert matred(s*s+3*sp.eye(2)) == sp.zeros(2)
    assert s.trace() == 0
    assert all(sp.rem(sp.Poly(x,t,modulus=3),
                      sp.Poly(t**3-t**2+t+1,t,modulus=3)).is_zero for x in s)
    residue = [[local(s[i,j])[0] for j in range(2)] for i in range(2)]
    assert residue == expected[name]
    print(name,"residue:",residue)

for transpose,v,w in ((False,(1,4),(0,1)),(True,(3,1),(1,0))):
    # H = u^2 R v + u^3 R^2 has this five-dimensional F7 basis.
    basis = [shifted(v,n) for n in (2,3,4)] + [shifted(w,n) for n in (3,4)]
    rows = [flat(q) for q in basis]
    assert rank(rows) == 5
    for x in basis:
        for y in basis:
            assert sub(mul(x[0],y[1]),mul(x[1],y[0])) == (0,)*5
    for name,s in operators.items():
        m = s.T if transpose else s
        m = [[local(m[i,j]) for j in range(2)] for i in range(2)]
        for q in basis:
            assert rank(rows+[flat(action(m,q))]) == 5
    print("transpose" if transpose else "matrix","convention:",
          "dimension 5; order",7**5,"; isotropic; stable under A/B/C")

sb = sp.Matrix(expected["B"])
wrong_line = sp.Matrix([1,2])
sv = list((sb*wrong_line).applyfunc(lambda x:x%7))
stv = list((sb.T*wrong_line).applyfunc(lambda x:x%7))
assert sv == [2,4] and stv == [3,0]
assert int(sp.Matrix.hstack(wrong_line,sb.T*wrong_line).det())%7 == 1
print("Displayed B eigenline:",sv,"under S;",stv,"under S^T; not S^T-stable")

# Rational homology, with columns (a,b) representing az+b.  Use the literal
# alternating form from the manuscript, not its former negative block matrix.
def multiplication(x):
    return sp.Matrix(6, 6, lambda i,j: sp.Poly(red(x*t**j),t).nth(i))

trace_basis = [multiplication(t**j).trace() for j in range(6)]
def trace(x):
    q = sp.Poly(red(x),t)
    return sum(q.nth(j)*trace_basis[j] for j in range(6))

T = sp.Matrix(6,6,lambda i,j:trace(t**(i+j)))
E = sp.BlockMatrix([[sp.zeros(6),-T],[T,sp.zeros(6)]]).as_explicit()
assert T.det() == 453789 == sp.discriminant(f.as_expr(),t)
F = sp.eye(12)[:,6:]  # f_j=(0,t^j), j=0,...,5
expected_det = {"A":-2**6*3**6*7**5*127,
                "B":-2**6*3**6*7**5,
                "C":-2**6*3**6*7**5*127}
vertices = dict(zip("ABC",(a,b,c)))
for name,s in operators.items():
    sigma = s.T
    U = matred((sigma-sp.eye(2))/2)
    assert all(sp.Poly(q,t).domain == sp.ZZ or
               all(k.q == 1 for k in sp.Poly(q,t).all_coeffs()) for q in U)
    assert matred(U*U+U+sp.eye(2)) == sp.zeros(2)
    Q = sp.BlockMatrix([[multiplication(sigma[i,j]) for j in range(2)]
                       for i in range(2)]).as_explicit()
    assert Q.T*E+E*Q == sp.zeros(12)
    assert Q*Q == -3*sp.eye(12)
    assert F.row_join(Q*F).det() != 0  # actual K-basis, not an assumption
    assert F.T*E*F == sp.zeros(6)
    gram = F.T*Q.T*E*F
    assert gram == sp.Matrix(6,6,lambda i,j:-trace(s[1,0]*t**(i+j)))
    assert gram.det() == expected_det[name]
    vertex = vertices[name].T
    assert matred(sigma*vertex-vertex*sigma) == sp.zeros(2)
    assert sp.rem(sp.Poly((vertices[name].trace()**2-4),t,modulus=3),
                  sp.Poly(t**3-t**2+t+1,t,modulus=3)).pow(13).rem(
                      sp.Poly(t**3-t**2+t+1,t,modulus=3)) == sp.Poly(-1,t,modulus=3)
    print(name,"integral O_K generator; Rosati; actual K-basis; det(H) =",gram.det())

    # Full local projectors for the corrected *homology* action.
    for sign in (1,-1):
        e = matred((2*sp.eye(2)+sign*sigma)/4)
        local_e = [[local(e[i,j]) for j in range(2)] for i in range(2)]
        generated = [action(local_e,shifted(v,n))
                     for v in ((1,0),(0,1)) for n in range(5)]
        assert rank([flat(v) for v in generated]) == 5
        for v in generated:
            assert action(local_e,v) == v
            for w in generated:
                assert sub(mul(v[0],w[1]),mul(v[1],w[0])) == (0,)*5

# Scope of the old obstruction: every pair has exactly one common residue
# line, but its free lifts disagree already modulo u^2. Enumerate all P^1(F7).
residue_lines = [(1,m) for m in range(7)]+[(0,1)]
def invariant_mod_u2(s,v):
    M = [[local(s[i,j]) for j in range(2)] for i in range(2)]
    image = action(M,v)
    return sub(mul(v[0],image[1]),mul(v[1],image[0]))[:2] == (0,0)

from itertools import combinations
for X,Y in combinations("ABC",2):
    common = []
    for v in residue_lines:
        def stable(name):
            m = operators[name].T
            w = [sum(local(m[i,j])[0]*v[j] for j in range(2))%7 for i in range(2)]
            return (v[0]*w[1]-v[1]*w[0])%7 == 0
        if stable(X) and stable(Y): common.append(v)
    assert common == [(1,5)]  # same line as (3,1)
lifts = {}
for name,s in operators.items():
    # Normalize as (3+a*u,1), convenient for the common vector (3,1).
    lifts[name] = [k for k in range(7)
                   if invariant_mod_u2(s.T,((3,k,0,0,0),(1,0,0,0,0)))]
assert lifts == {"A":[3],"B":[6],"C":[2]}
print("No common free rank-one summand; transposed first-order lifts:",lifts)
print("COMMON NONFREE H: order 3^3*7^5 = 453789; see human trace-pairing proof")
print("PRINCIPAL HOMOLOGY V56: PASS (special fibers only)")
