"""Short independent checks of the V141 four supplied homology matrices.

Standard-library integer arithmetic; no primary arithmetic imports. Rebuilds
q on the common O_L basis, checks four form equations and determinants,
enumerates all 4096 two-torsion vectors, and factors Phi42 over F2.
Does not independently reconstruct the complete integral Hom module.
"""
import hashlib
import json
import resource
import signal
import sys
import time
from pathlib import Path

resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(30)
if not __debug__:raise RuntimeError('Assertions required')
start=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
wb=(ROOT/'certificates/bc_polarization_square_data.json').read_bytes();w=json.loads(wb)
eb=(ROOT/'certificates/principal_unmarked_bc_data.json').read_bytes()
assert hashlib.sha256(eb).hexdigest()=='b9792b7ebe38ab8235d8b7b603a2596768d182b7863c6ca975eac56f9cb0437f'
E=[[int(x) for x in row] for row in json.loads(eb)['data']['principal_alternating_matrix']]
def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def transpose(A):return list(map(list,zip(*A)))
def determinant(A):
    A=[row[:] for row in A];sign=1;previous=1;n=len(A)
    for k in range(n-1):
        p=next(i for i in range(k,n) if A[i][k])
        if p!=k:A[k],A[p]=A[p],A[k];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*pivot-A[i][k]*A[k][j]
                assert numerator%previous==0
                A[i][j]=numerator//previous
        for i in range(k+1,n):A[i][k]=0
        previous=pivot
    return sign*A[-1][-1]
# Multiplication by q=4-t^2 on the ordered 1,t,...,t^5 basis.
f=[1,8,8,-6,-6,1]
columns=[]
for j in range(6):
    v=[0]*8;v[j]=4;v[j+2]=-1
    for k in range(7,5,-1):
        for i,a in enumerate(f):v[k-6+i]-=v[k]*a
    columns.append(v[:6])
small=transpose(columns)
Q=[[small[i%6][j%6] if i//6==j//6 else 0 for j in range(12)] for i in range(12)]
I=[[int(i==j) for j in range(12)] for i in range(12)]
assert mm(transpose(Q),E)==mm(E,Q) and determinant(Q)==1
kernels={}
for key,m in w['models'].items():
    e,f=map(int,key);T=[[int(x) for x in row] for row in m['matrix']]
    assert determinant(T)==64
    assert mm(mm(transpose(T),mm(E,Q if f else I)),T)==[[2*x for x in row] for row in mm(E,Q if e else I)]
    masks=[sum((a%2)<<j for j,a in enumerate(row)) for row in T]
    kernels[key]={v for v in range(4096) if all((v&r).bit_count()%2==0 for r in masks)}
    assert len(kernels[key])==64
assert kernels['00']==kernels['11']
assert kernels['01']==kernels['10']
assert kernels['00']&kernels['01']=={0}
assert len({x^y for x in kernels['00'] for y in kernels['01']})==4096
def rem(a,b):
    while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
    return a
def polymul(a,b):
    r=0
    while b:
        if b&1:r^=a
        b>>=1;a<<=1
    return r
phi=sum(1<<i for i in (0,1,3,4,6,8,9,11,12))
factors=[p for p in range(64,128) if rem(phi,p)==0]
assert len(factors)==2 and polymul(*factors)==phi
for p in factors:
    assert all(rem(p,d)!=0 for degree in (1,2,3) for d in range(1<<degree,1<<(degree+1)))
print(json.dumps({'status':'INDEPENDENT_FOUR_FORMS_AND_KERNELS_PASS','degrees':[64]*4,'kernel_sizes':{key:len(value) for key,value in kernels.items()},'two_complementary_kernel_classes':True,'two_torsion_vectors_exhausted':4096,'Phi42_mod2_irreducible_sextics_binary':factors,'elapsed_seconds':time.monotonic()-start,'python_version':sys.version,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'witness_sha256':hashlib.sha256(wb).hexdigest(),'adopted_principal_data_sha256':hashlib.sha256(eb).hexdigest(),'scope':'Four supplied form equations, degrees, finite kernel groups and mod2 factorization; norm-family completeness and full Hom membership retain their named inputs.'},indent=2,sort_keys=True))
