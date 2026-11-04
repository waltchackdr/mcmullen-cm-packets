#!/usr/bin/env python3
"""Independent CRT weights and integral pair-lattice gluing, without search.

Uses actual J7 differentials x^(j-1)dx/y, j=1,2,3, under x->zeta7*x.
Does not import a primary triple implementation or its output.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import json
import resource
import signal

ROOT=Path(__file__).resolve().parents[1]
CENSUS='experiments/fixed_field_cm_census_v98/independent.json'
CENSUS_SHA='8d2a5d643e3fca75f06270d9fa5d102d7d7711f1dfebd1dc7184c75445a1d4fe'
UNITS=tuple(a for a in range(42) if gcd(a,42)==1)
GROUP=tuple(product(range(2),range(6)))
CRT={h:next(a for a in UNITS if a%3==pow(2,h[0],3) and a%7==pow(3,h[1],7)) for h in GROUP}
LOG={a:h for h,a in CRT.items()}
LOG7={pow(3,j,7):j for j in range(6)}
R=(1,5,11,13,17,19)

def rref(matrix):
    a=[list(map(Fraction,row)) for row in matrix];piv=[];k=0
    for j in range(len(a[0])):
        q=next((i for i in range(k,len(a)) if a[i][j]),None)
        if q is None:continue
        a[k],a[q]=a[q],a[k];c=a[k][j];a[k]=[x/c for x in a[k]]
        for i in range(len(a)):
            if i!=k and a[i][j]:
                c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[k])]
        piv.append(j);k+=1
        if k==len(a):break
    return a,piv


def rank(columns):
    return len(rref(list(map(list,zip(*columns))))[1])


def integral(x):
    assert x.denominator==1
    return int(x)


def run():
    raw=(ROOT/CENSUS).read_bytes();assert sha256(raw).hexdigest()==CENSUS_SHA
    census=json.loads(raw)['orbits']
    types={f:set(census[i]['representative']) for f,i in [('B',0),('D',2)]}
    def field_column(a,cm):
        l=LOG[a]
        return tuple(1 if CRT[((l[0]-h[0])%2,(l[1]-h[1])%6)] in cm else -1 for h in GROUP)
    B=[field_column(a,types['B']) for a in R]
    D=[field_column(a,types['D']) for a in R]
    Y=[tuple(1 if pow(3,(LOG7[a]-h[1])%6,7) in {1,2,3} else -1 for h in GROUP) for a in (1,2,3)]
    rows=[list(row) for row in zip(*(B+Y))]
    reduced,piv=rref(rows)
    assert len(piv)==5 and all(j<6 for j in piv)
    H=[[0]*3 for _ in range(6)]
    for i,j in enumerate(piv):H[j]=[integral(reduced[i][6+k]) for k in range(3)]
    assert all(sum(B[j][i]*H[j][k] for j in range(6))==Y[k][i] for i in range(12) for k in range(3))
    free=next(j for j in range(6) if j not in piv)
    KB=[0]*6;KB[free]=1
    for i,j in enumerate(piv):KB[j]=-integral(reduced[i][free])
    assert all(sum(B[j][i]*KB[j] for j in range(6))==0 for i in range(12))
    # Universal gluing identities checked on the fifteen coordinate vectors;
    # individual pieces become kernel vectors whenever the input is a relation.
    for j in range(15):
        v=[int(i==j) for i in range(15)];hy=[sum(H[i][k]*v[12+k] for k in range(3)) for i in range(6)]
        bd=[v[i]+hy[i] for i in range(6)]+v[6:12]+[0]*3
        by=[-x for x in hy]+[0]*6+v[12:]
        assert [a+b for a,b in zip(bd,by)]==v
        assert all(sum((B+D+Y)[k][i]*by[k] for k in range(15))==0 for i in range(12))
    bd_reduced,bd_piv=rref(list(map(list,zip(*(B+D)))))
    bd_free=[j for j in range(12) if j not in bd_piv]
    bd_basis=[]
    for f in bd_free:
        z=[0]*12;z[f]=1
        for i,j in enumerate(bd_piv):z[j]=-integral(bd_reduced[i][f])
        bd_basis.append(z)
    joint_basis=[z+[0]*3 for z in bd_basis]
    joint_basis += [[-H[j][k] for j in range(6)]+[0]*6+[int(i==k)for i in range(3)]for k in range(3)]
    assert all(all(sum((B+D+Y)[j][i]*z[j] for j in range(15))==0 for i in range(12)) for z in joint_basis)
    result={'status':'EXACT_INDEPENDENT_CRT_INTEGRAL_LIFT_AND_GLUING_HUMAN_PROOF',
      'input_sha256':{CENSUS:CENSUS_SHA},'row_CRT_coordinates':GROUP,'row_units':list(CRT.values()),
      'column_order':[(f,a)for f in ('B','D')for a in R]+[('Y',a)for a in (1,2,3)],
      'types':{**{f:sorted(t)for f,t in types.items()},'Y_actual':[1,2,3]},
      'full_weight_columns':B+D+Y,'ranks':{'B':rank(B),'D':rank(D),'Y':rank(Y),'BD':rank(B+D),'BY':rank(B+Y),'BDY':rank(B+D+Y)},
      'B_pivot_columns':piv,'B_free_column':free,'B_integer_kernel_generator':KB,'Y_integral_lift_in_B':H,
      'BD_pivot_columns':bd_piv,'BD_free_columns':bd_free,'BD_integral_kernel_basis':bd_basis,
      'independent_joint_integral_basis':joint_basis,
      'joint_basis_argument':'The BD basis is the full integer kernel because its pivot map is integral and its free rows are identity. The three additional generators (-H e_k,0,e_k) split off arbitrary Y coordinates. The gluing formula therefore gives a unique integral expansion in these nine columns.',
      'exact_identity':'W_Y = W_B H, with H integral; independently checked in all twelve CRT rows.',
      'gluing':'Every joint relation(b,d,y) equals(b+Hy,d,0)+(-Hy,0,y). The first is in K(B,D), the second in K(B,Y), and both are integral. Thus the sum of the embedded pair kernels equals the full joint integer kernel, not just its rational span.',
      'intersection':'The intersection of the embedded pair kernels is precisely K(B), supported on B alone, of rank1. Pair ranks6 and4 therefore glue to joint rank9.',
      'actual_Y_mark':'Y is J(y^2=x^7-1) with zeta7 acting as x->zeta7*x. Its holomorphic differentials x^(j-1)dx/y have eigenvalues zeta7^j,j=1,2,3. Earlier type{1,4,5}=4*{1,2,3} relabels the same field action and permutes signed Y columns; no change of variety or independent Galois translation is used.',
      'scope':'Scientific arithmetic and conditional lattice gluing only. Algebraicity requires reviewed C46,C50 and the explicit contraction lemma with factor marks and rational descent. No new theorem promotion, geometry or paid compute.'}
    # Comparison occurs only AFTER the independent weights, integral lift and
    # integral bases have been constructed; no primary data guide them.
    primary_path=ROOT/'experiments/mixed_triple_v102/primary.json'
    if primary_path.exists():
        primary=json.loads(primary_path.read_bytes())
        row_perm=[primary['synchronized_galois_rows'].index(pow(CRT[h],-1,42)) for h in GROUP]
        primary_Y=[a%7 for f,a in primary['column_order'] if f=='Y']
        primary_to_actual=[2*a%7 for a in primary_Y]
        assert primary_to_actual==[2,3,1]
        col_perm=list(range(12))+[12+primary_to_actual.index(a)for a in (1,2,3)]
        assert all((B+D+Y)[j][i]==primary['full_weight_columns'][col_perm[j]][row_perm[i]] for j in range(15)for i in range(12))
        for z_primary in zip(*primary['integral_kernel_chart']):
            z=[z_primary[j]for j in col_perm]
            assert all(sum((B+D+Y)[j][i]*z[j]for j in range(15))==0 for i in range(12))
            hy=[sum(H[j][k]*z[12+k]for k in range(3))for j in range(6)]
            bd=[z[j]+hy[j]for j in range(6)]+z[6:12]
            coefficients=[bd[j]for j in bd_free]+z[12:]
            assert [sum(coefficients[k]*joint_basis[k][j]for k in range(9))for j in range(15)]==z
        result['primary_comparison']={'primary_sha256':sha256(primary_path.read_bytes()).hexdigest(),
          'independent_row_to_primary_row':row_perm,'independent_column_to_primary_column':col_perm,
          'primary_Y_to_actual_J7_labels':primary_to_actual,'all_180_weight_entries_equal':True,
          'all_nine_primary_chart_columns_integrally_reconstructed':True,
          'marking_argument':'The identity on the underlying J7 variety changes the marked generator from r7 to r7^4. H10 labels become4*{1,2,3}={1,4,5}; inverse4=2 sends primarylabels1,5,4 to actual2,3,1. This is a Ycolumn relabeling only. One common inverseCRT row permutation is used on B,D,Y simultaneously; no relative Galois twist occurs.'}
    return result

if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_CPU,(18,19));signal.alarm(30)
    print(json.dumps(run(),sort_keys=True,indent=2))
