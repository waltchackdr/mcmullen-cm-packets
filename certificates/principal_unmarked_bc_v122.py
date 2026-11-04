#!/usr/bin/env python3
"""Self-contained exact B/C principal-polarization separation candidate.

Actual geometric types and common lattice are established input hypotheses.
This reconstructs their CM modules and polarizations in coefficient arithmetic,
then proves the necessary trace-norm shell empty by finite rational recursion.
No repository arithmetic module is imported. Not promoted until review.
"""
import argparse, hashlib, json, math, sys
from fractions import Fraction as F
from pathlib import Path
import sympy as s

DATA=Path(__file__).with_name('principal_unmarked_bc_data.json')
DATA_SHA='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
MOD=(1,1,0,-1,-1,0,1,0,-1,-1,0,1,1)
N=12

def need(ok,message):
    if not ok:raise ValueError(message)

class K:
    """Q[z]/Phi_42(z); exact coefficient operations and rational inversion."""
    def __init__(self,value=0):
        if isinstance(value,K):self.c=value.c;return
        c=[F(v) for v in value] if isinstance(value,(tuple,list)) else [F(value)]
        c.extend([F(0)]*max(0,N-len(c)))
        for k in range(len(c)-1,N-1,-1):
            if c[k]:
                for j,a in enumerate(MOD[:-1]):
                    if a:c[k-N+j]-=c[k]*a
        self.c=tuple(c[:N])
    def __add__(self,o):return K([a+b for a,b in zip(self.c,K(o).c)])
    __radd__=__add__
    def __neg__(self):return K([-a for a in self.c])
    def __sub__(self,o):return self+-K(o)
    def __rsub__(self,o):return K(o)+-self
    def __mul__(self,o):
        c=[F(0)]*(2*N-1)
        for i,a in enumerate(self.c):
            if a:
                for j,b in enumerate(K(o).c):
                    if b:c[i+j]+=a*b
        return K(c)
    __rmul__=__mul__
    def __pow__(self,n):
        n=int(n);need(n>=0,'negative field exponent');r=K(1);a=self
        while n:
            if n&1:r=r*a
            a=a*a;n//=2
        return r
    def __eq__(self,o):return self.c==K(o).c
    def inverse(self):
        z=K([0,1]);columns=[(self*z**j).c for j in range(N)]
        a=[[columns[j][i] for j in range(N)]+[F(i==0)] for i in range(N)]
        for j in range(N):
            k=next((i for i in range(j,N) if a[i][j]),None);need(k is not None,'zero field inverse')
            a[j],a[k]=a[k],a[j];v=a[j][j];a[j]=[x/v for x in a[j]]
            for i in range(N):
                if i!=j and a[i][j]:
                    v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
        return K([row[-1] for row in a])

def verify_lattice_bridge(record):
    q=s.Matrix(record['data']['positive_definite_trace_Gram_Q']);r=record['lattice']
    u=s.Matrix(r['U_columns_old_coordinates']);h=s.Matrix(r['reduced_Gram_H'])
    mu=s.Matrix(r['gram_schmidt_mu']);d=list(map(s.Rational,r['gram_schmidt_positive_diagonal']))
    need(q.shape==u.shape==h.shape==mu.shape==(N,N),'lattice matrix dimensions')
    need(all(x.q==1 for x in u) and abs(u.det())==1,'non-unimodular transform')
    need(u.T*q*u==h,'Gram congruence')
    need(all(mu[i,j]==0 for i in range(N) for j in range(i,N)),'nontriangular LDL')
    need(all(v>0 for v in d),'nonpositive LDL')
    lower=s.eye(N)+mu
    need(lower*s.diag(*d)*lower.T==h,'LDL identity')
    bounds=r['reduced_coordinate_bounds'];inv=h.inv()
    need(len(bounds)==N and all(type(v) is int and v>=0 for v in bounds),'integer bounds')
    need(all((v+1)**2>12*inv[i,i] for i,v in enumerate(bounds)),'incomplete coordinate bounds')
    return h,[[F(v) for v in row] for row in mu.tolist()],list(map(F,d)),bounds

def reconstruct(record):
    z=K([0,1]);t=z+z**41
    need(z**42==1 and t**6+t**5-6*t**4-6*t**3+8*t**2+8*t+1==0,'field identities')
    powers=[z**i for i in range(N)]
    trace_powers=[sum((p*z**j).c[j] for j in range(N)) for p in powers]
    def trace(x):return sum(a*b for a,b in zip(x.c,trace_powers))
    def auto(x,g):return sum((a*z**((g*i)%42) for i,a in enumerate(x.c)),K())
    def bar(x):return auto(x,41)
    types=record['actual_tangent_types'];units=[i for i in range(42) if math.gcd(i,42)==1]
    need(types=={'B':[1,5,23,25,29,31],'C':[1,5,11,13,17,19]},'actual input types')
    twists=[g for g in units if {g*i%42 for i in types['C']}==set(types['B'])]
    need(twists==[31] and record['galois_exponent_B_to_C']==31,'Hodge twist direction')
    u=t-1;beta=-t**5+t**4+6*t**3-2*t**2-8*t-2;a=(beta*u**2).inverse();b=(u**3).inverse()
    delta=a*b;c={'B':t.inverse(),'C':(-1+5*t+4*t**2-t**3-t**4)*t.inverse()}
    d={'B':t,'C':-t+4*t**2+4*t**3-t**4-t**5};ideals={};xi={};data=record['data']
    target=s.Matrix(N,N,lambda i,j:0 if (i<6)==(j<6) else (-1 if i<6 else 1)*trace(delta*t**int(i%6+j%6))/2)
    need(target.det()==1,'principal polarization determinant')
    for name in ('B','C'):
        w=(z-d[name])*c[name].inverse();q=a*w;p=b*(1+3*w)
        ideals[name]=[t**j*q for j in range(6)]+[t**j*p for j in range(6)]
        xi[name]=-c[name]*(z-z**41).inverse()
        need([v.c for v in ideals[name]]==[tuple(map(F,v)) for v in data['I_'+name+'_basis']],'actual lattice '+name)
        need(xi[name].c==tuple(map(F,data['xi_'+name])),'polarization scalar '+name)
        gram=s.Matrix(N,N,lambda i,j:trace(xi[name]*ideals[name][i]*bar(ideals[name][j])))
        need(gram==target==s.Matrix(data['principal_alternating_matrix']),'actual alternating form '+name)
    r=auto(xi['B'],31)*xi['C'].inverse()
    need(r==bar(r) and r.c==tuple(map(F,data['relative_norm_target_r'])),'norm target')
    denominator=auto(ideals['B'][0],31).inverse();j=[v*denominator for v in ideals['C']]
    need([v.c for v in j]==[tuple(map(F,v)) for v in data['search_overlattice_J_basis']],'search over-lattice')
    ri=r.inverse();q=s.Matrix(N,N,lambda a,b:trace(j[a]*bar(j[b])*ri))
    need(q==s.Matrix(data['positive_definite_trace_Gram_Q']),'reconstructed trace Gram')

def finite_ball(record):
    h,mu,d,bounds=verify_lattice_bridge(record);x=[0]*N;leaves=[];nodes=0
    def visit(i,remaining):
        nonlocal nodes
        nodes+=1
        if i<0:
            need((s.Matrix(x).T*h*s.Matrix(x))[0]==12-remaining,'leaf norm')
            leaves.append((x[:],12-remaining));return
        shift=sum(mu[j][i]*x[j] for j in range(i+1,N))
        for v in range(-bounds[i],bounds[i]+1):
            x[i]=v;cost=d[i]*(F(v)+shift)**2
            if cost<=remaining:visit(i-1,remaining-cost)
    visit(N-1,F(12))
    need(leaves==[([0]*N,F(0))],'nonzero trace ball survivor')
    need(nodes==record['expected']['nodes_total']==59,'deterministic traversal')
    need(all(v.q==1 for v in h) and all(h[i,i]%2==0 for i in range(N)),'even integral Gram')
    witness=s.Matrix(record['expected']['minimum14_witness'])
    need(witness.shape==(N,1) and (witness.T*h*witness)[0]==14,'minimum14 witness')
    return {'status':'PRINCIPAL_UNMARKED_BC_EXACT_NORM_OBSTRUCTION_PASS','trace_norm_target':12,
            'lattice_minimum':14,'ball_nodes':nodes,'norm12_vectors':0,
            'scope':'Specified C56 principal B/C; polarization retained, O_L marking forgotten; geometric inputs are inherited.'}

def check(record):
    need(record['schema']=='PRINCIPAL_UNMARKED_BC_FINITE_NORM_V122_1','data schema')
    reconstruct(record);return finite_ball(record)

def main():
    need(not sys.flags.optimize,'optimized Python forbidden')
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=DATA);a=p.parse_args();raw=a.data.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==DATA_SHA,'frozen mathematical data hash')
    print(json.dumps(check(json.loads(raw)),sort_keys=True))
if __name__=='__main__':main()
