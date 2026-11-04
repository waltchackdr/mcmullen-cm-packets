"""Exact degree7 witness and small-prime norm obstruction for specified B/C.
Human saturation argument reduces arbitrary degrees to integer field norms;
adopted C60 excludes degree1. This is not a census of unpolarized maps.
"""
import argparse,hashlib,importlib.util,json,resource,signal,sys
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;DATA=HERE/'principal_bc_min_degree_data.json'
DATA_SHA='1e6b331975adac3175523d6a0a1de7152e4dc230be440e7e96a71d2243e95a6a'
C58_SOURCE_SHA='a73582098061839d7405a1fe00576b8bc50b703d0d830bb7d3977ecad446fdd3'
def need(ok,msg):
    if not ok:raise ValueError(msg)
def check(record):
    need(record['schema']=='PRINCIPAL_BC_MIN_DEGREE_V126_1','schema');p=HERE/'principal_unmarked_bc_v122.py';need(hashlib.sha256(p.read_bytes()).hexdigest()==C58_SOURCE_SHA,'C58 source pin');spec=importlib.util.spec_from_file_location('c58_primary',p);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    datasets={}
    for label,name,pin in [('C58','principal_unmarked_bc_data.json','b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'),('C59','principal_bc_isogeny_data.json','4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee'),('C60','principal_bc_unpolarized_data.json','f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e')]:
        raw=(HERE/name).read_bytes();need(hashlib.sha256(raw).hexdigest()==pin==record[label+'_data_sha256'],label+' data pin');datasets[label]=json.loads(raw)
    K=base.K;z=K([0,1]);zp=[z**i for i in range(12)];mat=lambda vs:s.Matrix.hstack(*(s.Matrix(v.c) for v in vs))
    def auto(v,g):return sum((c*z**((g*i)%42) for i,c in enumerate(v.c)),K())
    old=datasets['C58'];types=old['actual_tangent_types'];need(types=={'B':[1,5,23,25,29,31],'C':[1,5,11,13,17,19]},'actual types');need([g for g in (1,5,11,13,17,19,23,25,29,31,37,41) if {g*j%42 for j in types['C']}==set(types['B'])]==[31],'unique Hodge twist')
    d=old['data'];IB=[K(v) for v in d['I_B_basis']];IC=[K(v) for v in d['I_C_basis']];BB=mat(IB);BC=mat(IC)
    alpha=K(datasets['C60']['alpha']);m0=K(datasets['C59']['witness_m']);beta=1+z-z**3
    need(record['root_exponent']==5 and K(record['beta_coefficients'])==beta,'witness formula');m7=z**5*beta*alpha.inverse()*m0
    need(m7==K(record['m7']) and mat([beta*v for v in zp]).det()==7,'explicit coefficient and norm7 identities')
    T=BC.inv()*mat([m7*auto(v,31) for v in IB]);need(all(v.q==1 for v in T) and T==s.Matrix(record['matrix']),'full integral isogeny matrix');need(T.det()==7 and all(v.q==1 for v in 7*T.inv()),'degree7 and cyclic kernel')
    # Recheck exact saturation identities used by the human norm reduction;
    # maximal-order saturation equality itself is the accepted C60 dependency.
    SB=s.Matrix(datasets['C60']['S_B_basis']);SC=s.Matrix(datasets['C60']['S_C_basis']);basis=lambda A:[K(list(A[:,j])) for j in range(12)]
    for S,I in [(SB,BB),(SC,BC)]:
        U=S.inv()*I;need(all(v.q==1 for v in U) and abs(U.det())==27,'equal saturation indices')
    U=mat([alpha*v for v in basis(SC)]).inv()*mat([m0*auto(v,31) for v in basis(SB)]);need(all(v.q==1 for v in U) and abs(U.det())==1,'F0 saturation bridge')
    x=s.Symbol('x');f=s.cyclotomic_poly(42,x);gcds=[]
    for p,exponent in [(2,4),(3,3),(5,5)]:
        g=s.gcd(s.Poly(f,x,modulus=p),s.Poly(x**exponent-x,x,modulus=p));need(g.degree()==0,'no residue degree yielding prime norm<=6');gcds.append([p,exponent])
    need(record['expected_degree']==7 and record['kernel']=='Z/7Z','expected contract')
    return {'status':'PRINCIPAL_BC_MINIMUM_UNPOLARIZED_DEGREE7_INPUTS_PASS','minimum_degree':7,'kernel':'Z/7Z','finite_gcd_checks':gcds,'scalar_polarized':False,'proof_bridge':'C60 excludes1; equal saturation indices turn any degree into norm of integral eta, and finite polynomial gcds exclude prime ideal factors of norms2,3,4,5.','scope':'Specified C56 principal B/C, all unpolarized isogenies; no complete map count and no new curve point.'}
def main():
    need(not sys.flags.optimize,'optimized Python forbidden');resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30);ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=DATA);a=ap.parse_args();raw=a.data.read_bytes();need(hashlib.sha256(raw).hexdigest()==DATA_SHA,'data pin');print(json.dumps(check(json.loads(raw)),sort_keys=True))
if __name__=='__main__':main()
