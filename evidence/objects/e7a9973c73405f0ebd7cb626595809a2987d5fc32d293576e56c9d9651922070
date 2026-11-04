#!/usr/bin/env python3
import sympy as sp

t=sp.symbols('t')
f=sp.Poly(t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1,t,domain=sp.QQ)
def red(x): return sp.Poly(sp.expand(x),t,domain=sp.QQ).rem(f).as_expr()
def mm(X,Y): return (X*Y).applyfunc(red)
def mpow(X,n):
    R=sp.eye(2); Q=X
    while n:
        if n&1:R=mm(R,Q)
        Q=mm(Q,Q); n//=2
    return R

tinv=sp.invert(sp.Poly(t,t,domain=sp.QQ),f).as_expr()
h=sp.diag(t,1); hi=sp.diag(tinv,1)
B=sp.Matrix([[0,-1],[1,t]])
C=sp.Matrix([[t**5+t**4-4*t**3-4*t**2+2*t, t**4+t**3-4*t**2-5*t],
             [-t**4-t**3+4*t**2+5*t-1, -t**5-t**4+4*t**3+4*t**2-t]]).applyfunc(red)
D=sp.Matrix([[-t**5+2*t**4+7*t**3-6*t**2-13*t-1, 2*t**5+t**4-9*t**3-5*t**2+6*t],
             [-t**5+5*t**3+2*t**2-5*t+1, t**5-2*t**4-7*t**3+6*t**2+14*t+1]]).applyfunc(red)
for name,X in [('B',B),('C',C),('D',D)]:
    assert red(X.det()-1)==0
    print(name,'trace=',red(sp.trace(X)))
Bt,Ct,Dt=[mm(mm(h,X),hi) for X in (B,C,D)]
I=sp.eye(2)
def alpha0(X): return (I+X+mpow(X,4)).applyfunc(red)
def alpha_k(X,k): return mm(mpow(X,k),alpha0(X))
ivs=[ab for ab,m in sp.polys.polytools.intervals(f,eps=sp.Rational(1,10)**30)]
ivs=sorted(ivs,key=lambda ab:sp.Rational(ab[0]),reverse=True)
assert len(ivs)==6
assert all(sp.Rational(ivs[i][0])>sp.Rational(ivs[i+1][1]) for i in range(5))
def interval_eval(poly,ab):
    p=sp.Poly(poly,t,domain=sp.QQ); lo,hi=map(sp.Rational,ab); a=b=p.all_coeffs()[0]
    for c in p.all_coeffs()[1:]:
        vals=(a*lo,a*hi,b*lo,b*hi); a=min(vals)+c; b=max(vals)+c
    return a,b
def signs(poly):
    out=[]
    for ab in ivs:
        lo,hi=interval_eval(poly,ab); assert lo>0 or hi<0
        out.append(1 if lo>0 else -1)
    return tuple(out)
def norm(poly): return sp.Integer(sp.resultant(f.as_expr(),sp.Poly(poly,t,domain=sp.QQ).as_expr(),t))
maxd3=[1,11,16,17,19,20,22,32,37,38,40,41]
for name,X in [('B',Bt),('C',Ct),('D',Dt)]:
    print('\n',name)
    for k in maxd3:
        T=alpha_k(X,k); c=red(sp.trace(T)); pi=red(T.det())
        vals=[]
        for conj in (False,True):
            U=(I*c-T).applyfunc(red) if conj else T
            d=red(U[1,0])
            vals.append((conj,signs(d),norm(d)))
        hits=[v for v in vals if v[1]==(1,1,1,1,1,1)]
        if hits: print('HIT',k,hits)
        elif k in [1,17,20,22,38,41]: print('residual',k,vals)

# Structural B-D equivalence check.  B is the companion matrix for X^2-tX+1.
# P_D=[e1,D_t e1] conjugates B to D_t; h conjugates B to B_t.
# Hence P_D h^{-1} conjugates B_t to D_t.  Its determinant sign is the
# quotient of the two sign vectors and is totally positive.
dD=red(Dt[1,0])
assert signs(dD)==(1,1,-1,-1,-1,-1)
assert signs(t)==(1,1,-1,-1,-1,-1)
ratio=red(dD*tinv)
assert signs(ratio)==(1,1,1,1,1,1)
print('B-D normalized orientation class: TOTALLY POSITIVE; Norm(det P_D)=',norm(dD))
print('FOUR-VERTEX AUDIT PASSED')
