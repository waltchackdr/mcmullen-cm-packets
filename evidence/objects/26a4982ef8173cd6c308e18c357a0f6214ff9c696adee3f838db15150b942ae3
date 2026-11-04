#!/usr/bin/env python3
"""Exact all-42-trace embedded-order and affine-action identities.

Supports the human coverage lemmas, not completeness of proper module
classes or curve membership. Every operation is in QQ[t]/(f); all claimed
integral units and their inverses are checked in the integral power basis.
"""
from hashlib import sha256
import json
import sympy as sp

if not __debug__:
    raise RuntimeError("Assertions must be enabled")

T = sp.Symbol("T")
F = sp.Poly(T**6+T**5-6*T**4-6*T**3+8*T**2+8*T+1,T,domain=sp.QQ)
PI = T**4+T**3-4*T**2-2*T+5
ORDER_GROUPS = (
    ((1,0,0),(1,7,11,13,16,17,19,20,22,28,32,34,37,38,40,41)),
    ((7,1,1),(0,3,6,9,15,18,21,24,27,30,36,39)),
    ((7,1,2),(12,33)),
    ((41,12,1),(14,35)),
    ((41,22,1),(5,26)),
    ((43,20,1),(4,25)),
    ((43,31,1),(8,29)),
    ((83,36,1),(2,23)),
    ((127,94,1),(10,31)),
)


def poly(value):
    return sp.Poly(value,T,domain=sp.QQ).rem(F)


def reduce(value):
    return poly(value).as_expr()


def inverse(value):
    return sp.invert(poly(value),F).as_expr()


def integral(value):
    return all(c.q==1 for c in poly(value).all_coeffs())


def norm(value):
    return sp.Rational(sp.resultant(F,poly(value)))


def coefficients(value):
    p=poly(value)
    assert integral(value)
    return [int(p.nth(i)) for i in range(6)]


def multiply(left,right):
    a,b=left; c,d=right
    return reduce(a*c-b*d),reduce(a*d+b*c+T*b*d)


def main():
    powers=[(sp.Integer(1),sp.Integer(0))]
    for _ in range(45):
        powers.append(multiply(powers[-1],(0,1)))
    assert powers[42]==(1,0)
    alpha=[tuple(reduce(sum(powers[j][i] for j in (k,k+1,k+4)))
                 for i in range(2)) for k in range(42)]
    assert abs(norm(PI))==43 and abs(norm(T*T-4))==1
    assert sorted(k for _,indices in ORDER_GROUPS for k in indices)==list(range(42))
    rows=[]; traces=[]
    for k,(a,b) in enumerate(alpha):
        assert integral(a) and integral(b) and b!=0
        trace=reduce(2*a+T*b); traces.append(tuple(coefficients(trace)))
        assert reduce(a*a+T*a*b+b*b-PI)==0
        assert reduce(trace*trace-4*PI-b*b*(T*T-4))==0
        assert tuple(reduce(-v) for v in alpha[k])==alpha[(k+21)%42]
    assert len(set(traces))==42
    conjugates=[(reduce(a+T*b),reduce(-b)) for a,b in alpha]
    elements={tuple(tuple(coefficients(v)) for v in pair) for pair in alpha+conjugates}
    assert len(elements)==84
    assert abs(norm(alpha[0][1]))==7
    groups=[]
    for key,indices in ORDER_GROUPS:
        prime,root,exponent=key
        representative=indices[0]
        ref_a,ref_b=alpha[representative]
        assert abs(norm(ref_b))==prime**exponent
        if prime!=1:
            assert sp.isprime(prime)
            assert int(F.as_expr().subs(T,root))%prime==0
            assert int(ref_b.subs(T,root))%prime==0
            if exponent==1:
                # Integral principal ideal (b) of norm p lies in (p,t-r),
                # also of norm p, hence the two ideals are equal.
                assert abs(norm(ref_b))==prime
            else:
                assert key==(7,1,2) and abs(norm(T-1))==7
                ratio=reduce(ref_b*inverse((T-1)**2))
                assert integral(ratio) and integral(inverse(ratio))
                assert abs(norm(ratio))==1
        members=[]
        for k in indices:
            a,b=alpha[k]
            ratio=reduce(b*inverse(ref_b)); inv=inverse(ratio)
            assert integral(ratio) and integral(inv) and abs(norm(ratio))==1
            assert reduce(ratio*inv)==1
            shift=reduce(a-ratio*ref_a)
            assert integral(shift) and reduce(b-ratio*ref_b)==0
            rows.append(dict(trace_index=k,trace_coefficients=list(traces[k]),
                             alpha_a=coefficients(a),alpha_b=coefficients(b),
                             conductor_key=list(key),representative=representative,
                             unit_ratio=coefficients(ratio),inverse_unit=coefficients(inv),
                             affine_shift=coefficients(shift)))
            members.append(k)
        groups.append(dict(conductor_key=list(key),trace_indices=members))
    # Distinct roots over the same rational prime define distinct ideals;
    # the p7-square differs from p7 by its norm. No rational-norm shortcut.
    assert len({key for key,_ in ORDER_GROUPS})==9
    assert int(alpha[5][1].subs(T,12))%41==6
    assert int(alpha[14][1].subs(T,22))%41==17
    assert not integral(reduce(alpha[5][1]*inverse(alpha[14][1])))
    payload=dict(status="PASS_EXACT_PACKET_ORDER_COVERAGE",trace_count=42,norm_fiber_element_count=84,
                 maximal_order_trace_count=16,proper_embedded_order_count=8,
                 groups=groups,rows=sorted(rows,key=lambda r:r["trace_index"]),
                 scope="Embedded orders and affine-unit action identities; no proper-module count or new geometry decision")
    encoded=json.dumps(payload,sort_keys=True,separators=(",",":"))
    print(encoded)
    print("SHA256",sha256(encoded.encode()).hexdigest())


if __name__=="__main__":
    main()
