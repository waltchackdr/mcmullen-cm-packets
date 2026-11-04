"""Exact lattice inputs for the structural unpolarized B/C obstruction.

Human proof: a hypothetical isomorphism forces a CM unit ratio whose residue
has order4, whereas all roots of unity reduce to order dividing14. This file
certifies actual full multiplier order, saturation and all finite identities.
"""
import argparse,hashlib,importlib.util,json,resource,signal,sys
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form
HERE=Path(__file__).resolve().parent
DATA=HERE/'principal_bc_unpolarized_data.json'
DATA_SHA='f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e'
BASE_SHA='a73582098061839d7405a1fe00576b8bc50b703d0d830bb7d3977ecad446fdd3'
C59_SHA='4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee'
def need(ok,msg):
    if not ok:raise ValueError(msg)
def check(record):
    need(record['schema']=='PRINCIPAL_BC_UNPOLARIZED_V125_1','schema')
    p=HERE/'principal_unmarked_bc_v122.py';need(hashlib.sha256(p.read_bytes()).hexdigest()==BASE_SHA,'C58 source pin')
    spec=importlib.util.spec_from_file_location('accepted_c58',p);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    raw=base.DATA.read_bytes();need(hashlib.sha256(raw).hexdigest()==base.DATA_SHA==record['C58_data_sha256'],'C58 data pin');old=json.loads(raw);base.check(old)
    raw=(HERE/'principal_bc_isogeny_data.json').read_bytes();need(hashlib.sha256(raw).hexdigest()==C59_SHA==record['C59_data_sha256'],'C59 data pin');c59=json.loads(raw)
    K=base.K;z=K([0,1]);N=12;zp=[z**i for i in range(N)]
    def auto(v,g):return sum((c*z**((g*i)%42) for i,c in enumerate(v.c)),K())
    def mat(vs):return s.Matrix.hstack(*(s.Matrix(v.c) for v in vs))
    def basis(A):return [K(list(A[:,j])) for j in range(N)]
    def integer(A):return all(q.q==1 for q in A)
    def same(A,B):
        T=A.inv()*B;return integer(T) and abs(T.det())==1
    def hnf(A):
        den=s.ilcm(*[q.q for q in A]);return hermite_normal_form((den*A).applyfunc(int))/den
    def saturation(vs):return hnf(mat([q*v for q in zp for v in vs]))
    d=old['data'];IB=[K(v) for v in d['I_B_basis']];IC=[K(v) for v in d['I_C_basis']];BB=mat(IB);BC=mat(IC);CI=BC.inv()
    SB=saturation(IB);SC=saturation(IC)
    need(same(SB,s.Matrix(record['S_B_basis'])) and same(SC,s.Matrix(record['S_C_basis'])),'displayed maximal-order saturations')
    need(integer(SB.inv()*BB) and integer(SC.inv()*BC) and abs((SB.inv()*BB).det())==abs((SC.inv()*BC).det())==27,'saturation indices')
    # Full multiplier order: x is a multiplier iff A*x is integral.
    A=s.Matrix.hstack(*(s.Matrix(list(CI*mat([q*v for v in IC]))) for q in zp))
    den=s.ilcm(*[q.q for q in A]);H=hermite_normal_form((den*A).applyfunc(int).T);R=den*H.T.inv()
    need(same(R,s.Matrix(record['R_C_basis'])),'full multiplier order, not a suborder')
    need(integer(R) and abs(R.det())==27 and integer(A*R),'integral order index27')
    s7=sum((z**(6*i) for i in (1,2,4)),K())-sum((z**(6*i) for i in (3,5,6)),K());alpha=(1+s7)*K([s.Rational(1,2)]);baralpha=auto(alpha,41)
    need(s7*s7==K([-7]) and alpha*baralpha==K([2]),'quadratic norm2 element')
    need(s7==K(record['sqrt_minus7']) and alpha==K(record['alpha']),'explicit quadratic witness')
    m0=K(c59['witness_m']);T=CI*mat([m0*auto(v,31) for v in IB]);E=s.Matrix(d['principal_alternating_matrix'])
    need(T==s.Matrix(c59['witness_matrix']) and integer(T) and T.det()==64 and T.T*E*T==2*E,'actual C59 map is a local isomorphism at3')
    need(same(mat([alpha*v for v in basis(SC)]),mat([m0*auto(v,31) for v in basis(SB)])),'F0(SB)=alpha SC')
    pi=1-z**14;P=mat([pi*q for q in zp]);PI=P.inv()
    need(abs(P.det())==729 and integer(PI*s.Matrix(K([3]).c)),'prime3 modulus and characteristic3')
    need(all(integer(PI*s.Matrix((v-auto(v,41)).c)) for v in basis(R)),'every multiplier residue fixed by conjugation')
    need(integer(PI*mat([auto(v,41) for v in basis(P)])),'conjugation stable modulus')
    # OM/P = F3[zeta7], irreducible Phi7 since ord_7(3)=6.
    need([pow(3,k,7) for k in range(1,7)]==[3,2,6,4,5,1],'prime residue degree6')
    need(integer(PI*s.Matrix((z**14-1).c)),'all mu42 reductions have order dividing14')
    def zero_mod_P_at3(v):return all(q.q%3 for q in PI*s.Matrix(v.c))
    ratio=alpha*baralpha.inverse()
    need(zero_mod_P_at3(ratio-s7) and zero_mod_P_at3(ratio**2+1),'quadratic residue ratio of order4')
    need(not zero_mod_P_at3(K([2])) and zero_mod_P_at3(ratio**14+1) and not zero_mod_P_at3(ratio**14-1),'order4 incompatible with dividing14')
    return {'status':'PRINCIPAL_BC_UNPOLARIZED_OBSTRUCTION_INPUTS_PASS','saturation_indices':[27,27],'full_multiplier_order_index':27,'modulus_index':729,'saturation_identity':'F0(SB)=alpha SC','residue_ratio_order':4,'CM_unit_ratio_reduction_order_divides':14,'proof_scope':'Human local-unit argument excludes unpolarized B/C isomorphism for specified C56 models; unique rational Hodge twist and actual geometric identifications inherited from C58.','not_claimed':'Minimum degree among arbitrary unpolarized isogenies; new curve points; Mostaed image equality.'}
def main():
    need(not sys.flags.optimize,'optimized Python forbidden');resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=DATA);a=ap.parse_args();raw=a.data.read_bytes();need(hashlib.sha256(raw).hexdigest()==DATA_SHA,'witness data pin');print(json.dumps(check(json.loads(raw)),sort_keys=True))
if __name__=='__main__':main()
