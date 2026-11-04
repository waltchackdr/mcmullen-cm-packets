"""Exact short witness for the positive-unit curve-image comparison.

Checks field identities, both directions of principal lattice preservation,
and an explicit Gamma_std transport. Geometric equivariance is HUMAN-PROOF
in experiments/curve_twist_v138/IMAGE_PROPOSAL.md; no orbit exhaustiveness
or full theorem follows from this executable alone. Pending review/adoption.
"""
import resource,signal,json
resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
from pathlib import Path
import sympy as sp
if not __debug__:
 raise RuntimeError("Enabled assertions are required")
t=sp.Symbol('t');f=t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1
red=lambda x:sp.rem(sp.Poly(sp.expand(x),t,domain=sp.QQ),sp.Poly(f,t,domain=sp.QQ)).as_expr()
inv=lambda x:sp.invert(x,f,t)
mr=lambda m:m.applyfunc(red)
A=sp.Matrix([[t**5-4*t**3-t**2+3*t,-t**5+t**4+6*t**3-3*t**2-10*t],[t**3-1,-t**5+5*t**3+t**2-6*t]])
B=sp.Matrix([[0,-1],[1,t]])
C=sp.Matrix([[2*t-4*t**2-4*t**3+t**4+t**5,-5*t-4*t**2+t**3+t**4],[-1+5*t+4*t**2-t**3-t**4,-t+4*t**2+4*t**3-t**4-t**5]])
N=sp.diag(t,1);Ni=sp.diag(inv(t),1)
gens={k:mr(N*v*Ni) for k,v in zip('ABC',(A,B,C))};C=gens['C']
q=4-t*t
p0,p1=sp.Integer(2),t
for i in range(2,11):p0,p1=p1,red(t*p1-p0)
s=p1
assert red(s*s)==2-t
assert sp.resultant(f,s,t) in (1,-1)
assert sp.resultant(f,q,t)==1
h=mr(s*(sp.eye(2)+C))
assert red(h.det())==q
rho=lambda g:mr(inv(red(g.det()))*g.adjugate().T)
word='BAC'
g=sp.eye(2)
for letter in word:g=mr(g*gens[letter])
k=mr(g*h)
assert red(k.det())==q
# Exact global lattice, not only residue agreement.
beta=-t**5+t**4+6*t**3-2*t*t-8*t-2;u=t-1
a=inv(red(beta*u*u));b=inv(red(u**3));delta=red(a*b)
T=sp.Matrix([[a,3*b],[0,b]]);Ti=sp.Matrix([[red(beta*u*u),-3*red(beta*u*u)],[0,red(u**3)]])
K=mr(Ti*rho(k)*T)
Ki=mr(inv(red(K.det()))*K.adjugate())
assert red(K.det())==inv(q)
is_int=lambda m:all(c.q==1 for x in m for c in sp.Poly(x,t).all_coeffs())
assert is_int(K) and is_int(Ki)
W=sp.Matrix([[red(a*a),0],[3,1]]);Wi=sp.Matrix([[inv(red(a*a)),0],[-3*inv(red(a*a)),1]])
kstd=mr(W*k*Wi)
Dq=sp.diag(inv(q),1)
ambient=mr(Dq*kstd)
standard=mr(sp.diag(inv(delta),1)*ambient*sp.diag(delta,1))
assert red(standard.det())==1 and is_int(standard)
def poly(x):return [str(sp.Poly(red(x),t).nth(i)) for i in range(6)]
def matrix(m):return [[poly(m[i,j]) for j in range(2)] for i in range(2)]
report={'status':'CURVE_TWIST_MATRIX_WITNESS_PASS','field_polynomial':str(f),'s':poly(s),'q':poly(q),'compensating_word_left_to_right':word,'odd_projective_representative':matrix(h),'k_equals_word_times_h':matrix(k),'lattice_transition_Tinv_rho_k_T':matrix(K),'ambient_SL2O_matrix':matrix(standard),'identities':['s^2=2-t','det(h)=q','rho(word*h) P=P','Dq W(word*h)W^-1 in Gamma_std'],'scope':'Exact matrix identities only; geometric proof and independent review required.'}
print(json.dumps(report,indent=2,sort_keys=True))
