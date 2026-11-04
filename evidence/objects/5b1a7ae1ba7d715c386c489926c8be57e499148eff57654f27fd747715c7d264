"""Independent exact serialized-model audit. No producer imports/HNF/LLL.

Local CLI admits only the pinned actual model000 smoke. Remote CLI admits one
fixed ten-record original task stdout, with its inventory and byte bindings.
This verifies mathematics of those records, not provider/runtime provenance
or the authority to execute this new checker. No theta/period computation.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import re
import resource
import signal
import sys
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[2]
COMMON_PINS={
 'experiments/s_order_models_v161/model.json':'1803891819609ee4c8ecb84222a898a2574207f06bfdf08c48fca95ef27f8184',
 'experiments/s_order_models_v161/adoption.json':'bd08dfa1f408a3dc17385c6d6e76de27c35f21336f9411898cbdb68ceab1156b',
 'experiments/s_order_models_v161/catalog_candidate.json':'e6a6d431107ce1ceb9cee31ace70558876800d56daf3869a2174dfbea9eb97e3',
 'experiments/s_order_census_v162/review.json':'b8f0e1ef527b05b45601a41a6b0a4c25ef2e3b536144d85f485b3e898c68c049',
 'build/s-order-v162/models01/compact_inventory.json':'a8c4d55907e1be64106916d8b23a528f924da5bbe41ce9b04e3a8009792a2c3a',
}
EXPECTED_INPUT_PINS={
 'experiments/s_order_models_v161/model_candidate.py':'6bb0f3d1797a779ab8665ed47c3cf1d38e34b9205f7a5c9e19ab89615d54cf5d',
 'experiments/s_order_models_v161/model.json':COMMON_PINS['experiments/s_order_models_v161/model.json'],
 'experiments/s_order_models_v161/catalog_candidate.json':COMMON_PINS['experiments/s_order_models_v161/catalog_candidate.json'],
 'certificates/principal_unmarked_bc_data.json':'b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f',
 'certificates/principal_bc_unpolarized_data.json':'f4481b8162c2a8aac0db7dcc9afdbb0e8d6f56eb252e1ede92670034d6deb86e',
 'certificates/actual_bc_polarizations_data.json':'76c0520fbbae970dd2abfbe1b70d9d71df59534a8c7f4b1f20433047c83f311c',
}
AUTHORITY_SHA='ef4bd15b6401b6df390c1aefc85e9b14f850ca47c9a634b7059b23ab713e0476'
SMOKE_PATH='build/s-order-v162/models01/model000.json'
SMOKE_SHA='807940d7a8c65a58a86e367976feee0c30a406ebd1cbcdc42e63c2c620979a15'
RECORD_KEYS={'actual_tangent_type','catalog_index','delta_in_symplectic_basis','eligible_nonzero_pairs',
 'endomorphism_order','fiber_levels','fixed_torsion_numerator_basis','input_pins','invariant_refinement_chi',
 'lagrangian_rows','runtime','short_basis_change','short_basis_form','status','symplectic_change',
 'symplectic_power_basis','target_form','target_lattice_basis','theta_evaluations'}


def require(ok,message):
    if not ok:raise ValueError(message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def pairs_object(pairs):
    out={}
    for key,value in pairs:
        require(key not in out,'Duplicate JSON key');out[key]=value
    return out
def decode(raw):
    return json.loads(raw,object_pairs_hook=pairs_object,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON')))
def checked(path,digest,maximum=2*1024**2):
    path=Path(path)
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<=maximum,'Input file bound')
    raw=path.read_bytes();require(len(raw)<=maximum and sha(raw)==digest,'Input hash differs');return raw
def exact_matrix(raw,rows,cols,integers=False):
    require(type(raw) is list and len(raw)==rows,'Matrix row count')
    require(all(type(row) is list and len(row)==cols for row in raw),'Matrix column count')
    values=[]
    for row in raw:
        for x in row:
            require(type(x) is str and len(x)<=256 and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?',x),'Rational matrix syntax')
            value=s.Rational(x);require(not integers or value.q==1,'Nonintegral matrix');values.append(value)
    return s.Matrix(rows,cols,values)
def integral(m):return all(x.q==1 for x in m)
def zero_mod(m,p):return integral(m) and all(int(x)%p==0 for x in m)
def vector(raw,n,p):
    require(type(raw) is list and len(raw)==n and all(type(x) is int and 0<=x<p for x in raw),'Residue vector')
    return raw
def matvec(m,x):return [sum(a*b for a,b in zip(row,x)) for row in m]
def pairing(x,y):return sum(x[i]*y[i+6]-x[i+6]*y[i] for i in range(6))


def context():
    data={p:decode(checked(ROOT/p,h)) for p,h in COMMON_PINS.items()}
    seed=data['experiments/s_order_models_v161/model.json'];adoption=data['experiments/s_order_models_v161/adoption.json']
    catalog=data['experiments/s_order_models_v161/catalog_candidate.json']
    review=data['experiments/s_order_census_v162/review.json'];inventory=data['build/s-order-v162/models01/compact_inventory.json']
    require(inventory['status']=='COMPACT_RETAINED_MODEL_PROVENANCE_NOT_MATHEMATICAL_ADOPTION' and inventory['parent_inspection']=={'path':'build/s-order-v162/models01/inspection.json','sha256':'9f2aed79f0eed6ca8e80c3161a843dcfc59cccc8c22ae7606611e67cfc151604'},'Compact inventory parent binding')
    require(adoption['pins']['experiments/s_order_models_v161/model.json']==COMMON_PINS['experiments/s_order_models_v161/model.json'],'Seed adoption')
    require(review['verdict']=='ACCEPT' and review['pins']['experiments/s_order_models_v161/catalog_candidate.json']==COMMON_PINS['experiments/s_order_models_v161/catalog_candidate.json'],'Catalog review binding')
    require(catalog['input_model_sha256']==COMMON_PINS['experiments/s_order_models_v161/model.json'],'Catalog seed')
    require([x['index'] for x in catalog['models']]==list(range(160)),'Catalog coverage')
    require(inventory['authority_sha256']==AUTHORITY_SHA and [x['index'] for x in inventory['models']]==list(range(160)),'Inventory authority/coverage')
    for i,item in enumerate(inventory['models']):
        require((item['job'],item['task'],item['line'])==(i//80,(i%80)//10,i%10+1),'Inventory task partition')
    e=exact_matrix(seed['ambient_form'],12,12);pi=exact_matrix(seed['source_lattice_basis'],12,12,True)
    lift=exact_matrix(seed['quotient_lifts'],12,6,True);q=exact_matrix(seed['quotient_projection'],6,12,True)
    omega=exact_matrix(seed['discriminant_gram'],6,6,True);real=exact_matrix(seed['real_subfield_action'],6,6,True)
    require(e.T==-e and abs(pi.det())==729 and integral(3*pi.inv()),'Seed lattice')
    require(zero_mod(q*pi,3) and zero_mod(q*lift-s.eye(6),3),'Seed quotient')
    require(zero_mod(omega-3*lift.T*e*lift,3) and int(omega.det())%3!=0,'Seed pairing')
    j=s.zeros(12);j[:6,6:]=s.eye(6);j[6:,:6]=-s.eye(6)
    return {'seed':seed,'catalog':catalog,'inventory':inventory,'e':e,'pi':pi,'lift':lift,'q':q,'omega':omega,'real':real,'j':j,'f':s.eye(12)-pi}


def audit_record(ctx,d,index):
    require(type(index) is int and index in range(160),'Fixed catalog domain')
    require(type(d) is dict and set(d)==RECORD_KEYS,'Record schema')
    require(type(d['catalog_index']) is int and d['catalog_index']==index,'Catalog index')
    require(d['input_pins']==EXPECTED_INPUT_PINS,'Record source inputs')
    require(d['status']=='EXACT_MODEL_CANDIDATE_NOT_THETA_EVALUATION' and type(d['theta_evaluations']) is int and d['theta_evaluations']==0,'Model-only status')
    require(d['actual_tangent_type']==[1,5,11,13,17,19] and d['actual_tangent_type']==ctx['seed']['actual_tangent_type'],'Actual type')
    entry=ctx['catalog']['models'][index]
    require(d['lagrangian_rows']==entry['lagrangian_rows'],'Exact catalog representative')
    rows=d['lagrangian_rows']
    for row in rows:vector(row,6,3)
    require(len(rows)==3,'Lagrangian dimension')
    hp={tuple(sum(a[k]*rows[k][i] for k in range(3))%3 for i in range(6)) for a in itertools.product(range(3),repeat=3)}
    require(len(hp)==27,'Independent Lagrangian generators')
    h=s.Matrix(rows).T
    require(zero_mod(h.T*ctx['omega']*h,3),'Actual isotropy')
    real=[[int(x)%3 for x in row] for row in ctx['real'].tolist()]
    stable=all(tuple(x%3 for x in matvec(real,v)) in hp for v in hp)
    require(d['endomorphism_order']==entry['endomorphism_order']==('R' if stable else 'S'),'Full order via accepted catalog/multiplier theorem')
    mat=lambda key:exact_matrix(d[key],12,12,True)
    lat=mat('target_lattice_basis');lf=mat('target_form')
    require(abs(lat.det())==27,'Target index in O')
    li=lat.inv()
    require(integral(li*ctx['pi']) and abs((li*ctx['pi']).det())==27,'P inclusion and quotient degree')
    require(integral(li*ctx['lift']*h),'Full chosen inverse-image lattice inclusion')
    # The generated P+VH lattice already has index27 in O because H has
    # order27 in O/P. Equal index proves equality without computing HNF.
    require(lf==lat.T*ctx['e']*lat and lf.T==-lf and lf.det()==1,'Principal form identity')
    change=mat('short_basis_change');short=mat('short_basis_form');symp=mat('symplectic_change');u=mat('symplectic_power_basis');act=mat('delta_in_symplectic_basis')
    require(abs(change.det())==1 and short==change.T*lf*change,'Short basis identity')
    require(abs(symp.det())==1 and u==lat*symp and symp.T*lf*symp==ctx['j'],'Symplectic lattice equality')
    require(u.T*ctx['e']*u==ctx['j'] and ctx['f']*u==u*act,'Power-basis form/action bridge')
    require(act.T*ctx['j']*act==ctx['j'] and act**2+act+s.eye(12)==s.zeros(12),'Integral polarized cubic action')
    chi=vector(d['invariant_refinement_chi'],12,2)
    # Equality of homogeneous quadratic polynomials modulo2 is determined
    # by 12 diagonal and 66 cross coefficients; no 4096-point replay.
    form=ctx['j']*act
    def coefficients(a):return [int(a[i,i])%2 for i in range(12)]+[int(a[i,j]+a[j,i])%2 for i in range(12) for j in range(i+1,12)]
    expected=chi+[int(j==i+6) for i in range(12) for j in range(i+1,12)]
    require(coefficients(form)==expected and coefficients(act.T*form*act)==expected,'Full invariant quadratic refinement')
    require(sum(chi[i]*chi[i+6] for i in range(6))%2==0,'Refinement parity')
    quadratic=lambda x:(sum(x[i]*x[i+6] for i in range(6))+sum(a*b for a,b in zip(chi,x)))%2
    nums=exact_matrix(d['fixed_torsion_numerator_basis'],12,6,True)
    require(zero_mod((act-s.eye(12))*nums,3),'Fixed numerator columns')
    nn=[[int(x)%3 for x in row] for row in nums.tolist()];aa=[[int(x) for x in row] for row in act.tolist()]
    points={tuple(sum(nn[i][k]*a[k] for k in range(6))%3 for i in range(12)) for a in itertools.product(range(3),repeat=6)}
    require(len(points)==729,'Complete fixed group')
    sectors={0:set(),2:set(),4:set()}
    for point in points:
        delta=[x-y for x,y in zip(matvec(aa,point),point)]
        require(all(x%3==0 for x in delta),'Fixed point lift')
        k=[x//3 for x in delta]
        exponent=(3*quadratic(k)+pairing(point,k))%6
        require(exponent in sectors,'Order-three normalized fiber character')
        sectors[exponent].add(point)
    counts=[len(sectors[x]) for x in (0,2,4)]
    require(d['fiber_levels']==counts==[261,234,234],'Actual fiber sectors')
    eligible=sectors[0]
    require(all(tuple(-x%3 for x in v) in eligible for v in eligible),'Eligible negation')
    pairs=sorted(v for v in eligible if any(v) and v<tuple(-x%3 for x in v))
    require(type(d['eligible_nonzero_pairs']) is list and len(d['eligible_nonzero_pairs'])==130,'Pair count')
    for v in d['eligible_nonzero_pairs']:vector(v,12,3)
    require(pairs==[tuple(v) for v in d['eligible_nonzero_pairs']],'Every exact eligible pair label')
    require(type(d['runtime']) is dict and set(d['runtime'])=={'python','sympy'} and all(type(v) is str for v in d['runtime'].values()),'Observed runtime metadata schema')
    return {'catalog_index':index,'status':'EXACT_SERIALIZED_MODEL_CHECK_PASS','endomorphism_order':d['endomorphism_order'],
            'lattice_and_form':'PASS','full_order_criterion':'PASS_USING_ACCEPTED_CATALOG_MULTIPLIER_THEOREM',
            'symplectic_action':'PASS','quadratic_refinement_all_coefficients':'PASS','fixed_points_checked':729,
            'fiber_counts':counts,'all_130_pairs_match':True,'runtime_metadata':d['runtime'],'theta_evaluations':0}


def audit_line(ctx,line,index,stdout_sha256):
    """Audit one original stdout line after caller hashes its complete chunk.

    Caller must bind the whole raw chunk to stdout_sha256 and enforce its
    fixed ten-record range. This API checks the reviewed inventory's expected
    stdout digest and the model's exact line+LF digest independently.
    """
    require(type(index) is int and index in range(160),'Fixed catalog domain')
    require(type(line) is bytes and len(line)<=2*1024**2 and b'\n' not in line and b'\r' not in line,'One raw model line')
    item=ctx['inventory']['models'][index]
    require(stdout_sha256==item['stdout_sha256'],'Original stdout digest')
    model_sha=sha(line+b'\n')
    require(model_sha==item['sha256'],'Collected model byte identity')
    begun=time.perf_counter();record=audit_record(ctx,decode(line),index)
    record.update(model_sha256=model_sha,elapsed_seconds=time.perf_counter()-begun,
                  raw_provenance={k:item[k] for k in ('job','task','stdout_sha256','line')})
    return record


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--input',type=Path);parser.add_argument('--input-sha256');parser.add_argument('--start',type=int)
    args=parser.parse_args()
    if args.smoke:
        require(args.input is None and args.input_sha256 is None and args.start is None,'Smoke has one fixed input')
        resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
        path=ROOT/SMOKE_PATH;digest=SMOKE_SHA;start=0;count=1
    else:
        require(platform.system()=='Linux' and bool(os.environ.get('RUNPOD_POD_ID')),'Remote chunks only')
        require(args.start in range(0,160,10) and args.input is not None and re.fullmatch('[0-9a-f]{64}',args.input_sha256 or ''),'Fixed chunk arguments')
        resource.setrlimit(resource.RLIMIT_CPU,(120,120));signal.alarm(150)
        resource.setrlimit(resource.RLIMIT_FSIZE,(2*1024**2,2*1024**2))
        path=args.input;digest=args.input_sha256;start=args.start;count=10
    begun=time.perf_counter();cpu=time.process_time();ctx=context()
    raw=checked(path,digest);lines=raw.splitlines()
    require(len(lines)==count and raw.endswith(b'\n'),'Exact complete task lines')
    if not args.smoke:require(all(ctx['inventory']['models'][i]['stdout_sha256']==digest for i in range(start,start+10)),'Original task stdout binding')
    records=[]
    for offset,line in enumerate(lines):
        index=start+offset
        records.append(audit_line(ctx,line,index,ctx['inventory']['models'][index]['stdout_sha256'] if args.smoke else digest))
    result={'status':'INDEPENDENT_EXACT_DATA_CHUNK_PASS_NOT_ADOPTION','mode':'ONE_ACTUAL_MODEL_SMOKE' if args.smoke else 'REMOTE_FIXED_TEN_MODEL_CHUNK',
            'checker_sha256':sha(Path(__file__).read_bytes()),'input_sha256':digest,'indices':list(range(start,start+count)),
            'common_pins':COMMON_PINS,'record_expected_input_pins':EXPECTED_INPUT_PINS,'execution_authority_sha256':AUTHORITY_SHA,
            'records':records,'wall_seconds':time.perf_counter()-begun,'cpu_seconds':time.process_time()-cpu,
            'checker_runtime':{'python':sys.version,'sympy':s.__version__},
            'method':'Exact serialized identities, accepted catalog order proof, quadratic coefficient identities, actual-numerator AH labels; no producer imports/HNF/LLL',
            'provenance_runtime_adoption':False,'mathematical_data_adoption':False,'theta_evaluations':0}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
