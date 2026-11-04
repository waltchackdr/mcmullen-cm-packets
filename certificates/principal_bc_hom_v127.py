"""Exact complete normalized Hom ideal and polarization identities for B/C.
All eta satisfy finite linear integral-map constraints; no search truncation.
The nonprincipal idealclass and degree formula use the stated human bridges.
"""
import argparse,hashlib,importlib.util,json,resource,signal,sys
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form,smith_normal_form
HERE=Path(__file__).resolve().parent;DATA=HERE/'principal_bc_hom_data.json';DATA_SHA='9f815f1592c8e22a2a14b15086d8d9cb63092f9ab84bb90a3d5ee4935b3fe22c';BASE_SHA='a73582098061839d7405a1fe00576b8bc50b703d0d830bb7d3977ecad446fdd3'
PINS={'principal_unmarked_bc_data.json':'b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f','principal_bc_isogeny_data.json':'4796d9f3b2ddfe3fefbf80e6fff1b7842c7e7d40cd35b3795deec0974ab49cee','principal_bc_unpolarized_data.json':'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e','principal_bc_min_degree_data.json':'1e6b331975adac3175523d6a0a1de7152e4dc230be440e7e96a71d2243e95a6a'}
def need(ok,msg):
    if not ok:raise ValueError(msg)
def check(record):
    need(record['schema']=='PRINCIPAL_BC_HOM_IDEAL_V127_1' and record['input_pins']==PINS,'schema and inherited pins');p=HERE/'principal_unmarked_bc_v122.py';need(hashlib.sha256(p.read_bytes()).hexdigest()==BASE_SHA,'C58 arithmetic source');spec=importlib.util.spec_from_file_location('c58',p);base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base);values=[]
    for n,h in PINS.items():
        raw=(HERE/n).read_bytes();need(hashlib.sha256(raw).hexdigest()==h,'inherited data');values.append(json.loads(raw))
    c58,c59,c60,c61=values;K=base.K;z=K([0,1]);zp=[z**i for i in range(12)];mat=lambda vs:s.Matrix.hstack(*(s.Matrix(v.c) for v in vs));basis=lambda A:[K(list(A[:,j])) for j in range(12)]
    def auto(v,g):return sum((c*z**((g*i)%42) for i,c in enumerate(v.c)),K())
    def integral(A):return all(q.q==1 for q in A)
    def same(A,B):
        U=A.inv()*B;return integral(U) and abs(U.det())==1
    def lattice(vs):
        A=mat(vs);den=s.ilcm(*[q.q for q in A]);return hermite_normal_form((den*A).applyfunc(int))/den
    d=c58['data'];IB=[K(v) for v in d['I_B_basis']];IC=[K(v) for v in d['I_C_basis']];CI=mat(IC).inv();m0=K(c59['witness_m']);alpha=K(c60['alpha']);factor=m0*alpha.inverse();sb=[auto(v,31) for v in IB]
    A=s.Matrix.hstack(*(s.Matrix(list(CI*mat([factor*q*v for v in sb]))) for q in zp));den=s.ilcm(*[q.q for q in A]);H=hermite_normal_form((den*A).applyfunc(int).T);fullHom=den*H.T.inv()
    t=z+z**41;pi=1-z**14;OL=[t**i for i in range(6)];P=[pi*q for q in zp];R=s.Matrix(c60['R_C_basis']);need(same(R,lattice(OL+P)),'R=OL+P');J=lattice([alpha*v for v in OL]+P)
    need(same(J,s.Matrix(record['J_basis'])) and same(J,fullHom),'all integral maps equal normalized J, not just inclusion');need(integral(J) and abs(J.det())==abs(R.det())==27 and integral(A*J),'integral order and ideal indices')
    Jv=basis(J);Rv=basis(R);need(integral(J.inv()*mat([a*v for a in Rv for v in Jv])),'R module closure')
    need(same(lattice([a*auto(v,41) for a in Jv for v in Jv]),R),'J barJ=R explicit inverse')
    gamma=z**6*(1-z**2);need(gamma==K(record['square_generator']) and mat([gamma*v for v in zp]).det()==1,'unit square generator');need(same(lattice([a*v for a in Jv for v in Jv]),mat([gamma*v for v in Rv])),'J squared=gamma R')
    r=K(d['relative_norm_target_r']);need(factor*auto(factor,41)==r,'target relative norm normalization');need(record['source_pullback_twist']==19 and (19*31)%42==1,'source inverse twist')
    xiB=K(d['xi_B']);xiC=K(d['xi_C']);need(auto(xiC*r,19)==xiB,'all source pullbacks use sigma19 of eta bar_eta by trace invariance')
    # Human degree bridge uses exact maximal-order saturations inherited from C60.
    SB=s.Matrix(c60['S_B_basis']);SC=s.Matrix(c60['S_C_basis']);need(same(mat([factor*auto(v,31) for v in basis(SB)]),SC),'normalized rational map identifies saturations')
    for S,I in [(SB,mat(IB)),(SC,mat(IC))]:need(integral(S.inv()*I) and abs((S.inv()*I).det())==27,'equal saturation indices for all degrees')
    # Explicit degree7 pullback polarization type, supplementing C61.
    T7=s.Matrix(c61['matrix']);m7=K(c61['m7']);need(T7==CI*mat([m7*v for v in sb]) and T7.det()==7,'C61 actual map');E=s.Matrix(d['principal_alternating_matrix']);D=smith_normal_form(T7.T*E*T7,domain=s.ZZ);need([abs(int(D[i,i])) for i in range(12)]==[1]*10+[7]*2,'degree7 pullback type')
    need(record['degree7_polarization_type']==[1]*5+[7],'type contract')
    return {'status':'PRINCIPAL_BC_COMPLETE_HOM_IDEAL_INPUTS_PASS','J_index_OM':27,'full_Hom':'eta in J=alpha OL+P; F_eta=(m0/alpha)*eta*sigma31','degree':'N_M/Q(eta) for eta nonzero','actual_source_pullback_multiplier':'sigma19(eta*bar(eta))','inverse_ideal':'bar(J)','square_ideal':'z^6*(1-z^2)*R','ideal_class_order':2,'class_order_bridge':'J invertible by JbarJ=R; J^2 principal; if J principal index ratio forces norm1 and a forbidden C60 isomorphism.','degree7_pullback_type':[1,1,1,1,1,7],'scope':'Specified C56 models; no total Picard group order or surjectivity onto all totally positive multipliers asserted.'}
def main():
    need(not sys.flags.optimize,'optimized Python forbidden');resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30);ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=DATA);a=ap.parse_args();raw=a.data.read_bytes();need(hashlib.sha256(raw).hexdigest()==DATA_SHA,'witness data pin');print(json.dumps(check(json.loads(raw)),sort_keys=True))
if __name__=='__main__':main()
