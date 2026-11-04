#!/usr/bin/env python3
"""Exact scalar-polarized B/C isogenies; uses the accepted C58 primary field path.

Completeness is a finite rational suffix recursion, followed by the full norm,
integral map and polarization tests. The separate checker uses companion
matrices and Schur prefix bounds. Human geometric hypotheses are inherited.
"""
import argparse,hashlib,importlib.util,json,math,resource,signal,sys
from fractions import Fraction as F
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
DATA=HERE/'principal_bc_isogeny_data.json'
DATA_SHA='4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee'
BASE_SHA='a73582098061839d7405a1fe00576b8bc50b703d0d830bb7d3977ecad446fdd3'
def need(ok,message):
    if not ok: raise ValueError(message)
def load_base():
    p=HERE/'principal_unmarked_bc_v122.py'
    need(hashlib.sha256(p.read_bytes()).hexdigest()==BASE_SHA,'C58 primary source pin')
    spec=importlib.util.spec_from_file_location('accepted_c58_primary',p)
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    return base

def check(record):
    need(record['schema']=='PRINCIPAL_BC_MULTIPLIER2_V123_1','schema')
    base=load_base();raw=base.DATA.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==base.DATA_SHA==record['C58_data_sha256'],'C58 data pin')
    old=json.loads(raw);need(base.check(old)['lattice_minimum']==14,'multiplier1 obstruction')
    K=base.K;z=K([0,1]);N=12
    powers={g:[z**((g*i)%42) for i in range(N)] for g in (31,41)}
    def auto(v,g):return sum((c*w for c,w in zip(v.c,powers[g])),K())
    def mat(vs):return s.Matrix.hstack(*(s.Matrix(v.c) for v in vs))
    d=old['data'];lat=old['lattice']
    IB=[K(v) for v in d['I_B_basis']];IC=[K(v) for v in d['I_C_basis']]
    BB=mat(IB);BC=mat(IC);BI=BB.inv();CI=BC.inv();sb=[auto(v,31) for v in IB]
    J=[K(v) for v in d['search_overlattice_J_basis']];r=K(d['relative_norm_target_r'])
    E=s.Matrix(d['principal_alternating_matrix']);U=s.Matrix(lat['U_columns_old_coordinates'])
    H,mu,diag,_=base.verify_lattice_bridge(old);inverse=H.inv()
    bounds=[math.isqrt(int(s.floor(24*inverse[i,i]))) for i in range(N)]
    need(bounds==[1,1,1,1,2,1,1,1,1,1,1,1],'complete multiplier2 coordinate bounds')
    vector=[0]*N;nodes=0;shell=[];found=[];bad_norm=bad_integral=0
    def visit(i,remaining):
        nonlocal nodes,bad_norm,bad_integral
        nodes+=1
        if i<0:
            need((s.Matrix(vector).T*H*s.Matrix(vector))[0]==24-remaining,'leaf trace norm')
            if remaining:return
            coords=U*s.Matrix(vector)
            m=sum((int(c)*w for c,w in zip(coords,J)),K())
            shell.append(tuple(vector))
            if m*auto(m,41)!=2*r:bad_norm+=1;return
            T=CI*mat([m*w for w in sb])
            if any(v.q!=1 for v in T):bad_integral+=1;return
            need(T.T*E*T==2*E and T.det()==64,'direct polarization and degree')
            need(all(v.q==1 for v in 2*T.inv()),'kernel killed by2')
            found.append((m,T));return
        shift=sum(mu[j][i]*vector[j] for j in range(i+1,N))
        for value in range(-bounds[i],bounds[i]+1):
            vector[i]=value;cost=diag[i]*(F(value)+shift)**2
            if cost<=remaining:visit(i-1,remaining-cost)
    visit(N-1,F(24))
    expected=record['expected']
    need((nodes,len(shell),bad_norm,bad_integral,len(found))==(1033,84,72,6,6),'complete deterministic counts')
    need(expected=={'multiplier':2,'degree':64,'trace_shell':84,'norm_rejected':72,'integrality_rejected':6,'isogenies':6,'nodes':1033,'automorphism_order':6},'expected count contract')
    actual={m.c for m,_ in found}
    need(actual=={K(v).c for v in record['expected_scalar_set']},'complete scalar set')
    m0=K(record['witness_m']);chosen=[T for m,T in found if m==m0]
    need(len(chosen)==1 and chosen[0]==s.Matrix(record['witness_matrix']),'displayed witness')
    unit=z**7
    need({(m0*unit**j).c for j in range(6)}==actual,'single explicit target unit orbit')
    GB=BI*mat([unit*v for v in IB]);GC=CI*mat([unit*v for v in IC])
    for G in (GB,GC):
        need(all(v.q==1 for v in G) and G.det()==1 and G.T*E*G==E,'integral polarized unit')
        need(G**6==s.eye(N) and all(G**j!=s.eye(N) for j in (1,2,3)),'exact order6')
    need(all(T*GB==GC*T for _,T in found),'unit intertwining')
    return {'status':'PRINCIPAL_BC_SCALAR_POLARIZED_ISOGENY_PASS','minimum_positive_scalar_multiplier':2,
            'degree':64,'exact_number_of_multiplier2_isogenies':6,'kernel':'(Z/2Z)^6',
            'explicit_target_automorphism_orbits':1,'polarized_automorphism_orders':[6,6],
            'proof_bridge':'Complete six-map list and faithful pre/postcomposition prove full origin-preserving polarized automorphism groups are the displayed mu6. All maps share their kernel.',
            'scope':'Specified C56 principal B/C. Does not minimize arbitrary unpolarized isogeny degree.'}

def main():
    need(not sys.flags.optimize,'optimized Python forbidden')
    resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=DATA);a=p.parse_args();raw=a.data.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==DATA_SHA,'witness data pin')
    print(json.dumps(check(json.loads(raw)),sort_keys=True))
if __name__=='__main__':main()
