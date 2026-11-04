const fs = require('fs');
const data = JSON.parse(fs.readFileSync('/mnt/data/qag_paper_a_v1/pp_modular_data.json','utf8'));
const G=data.Gpp, U=data.Upp, N=12;
const mod=(a,p)=>((a%p)+p)%p;
function rowMul(v,M,p){let o=Array(N).fill(0);for(let j=0;j<N;j++){let s=0;for(let i=0;i<N;i++)s+=v[i]*M[i][j];o[j]=mod(s,p);}return o;}
function rref(rows,p,ncols=N){let A=rows.filter(r=>r.some(x=>mod(x,p)!==0)).map(r=>r.map(x=>mod(x,p)));let rr=0,piv=[];for(let c=0;c<ncols;c++){let q=-1;for(let i=rr;i<A.length;i++)if(A[i][c]){q=i;break;}if(q<0)continue;[A[rr],A[q]]=[A[q],A[rr]];let inv=1; if(p===3 && A[rr][c]===2) inv=2;A[rr]=A[rr].map(x=>mod(inv*x,p));for(let i=0;i<A.length;i++)if(i!==rr&&A[i][c]){let f=A[i][c];for(let j=0;j<ncols;j++)A[i][j]=mod(A[i][j]-f*A[rr][j],p);}piv.push(c);rr++;if(rr===A.length)break;}return {rows:A.slice(0,rr),piv};}
function nullspace(rows,p,ncols=N){const z=rref(rows,p,ncols), free=[];for(let c=0;c<ncols;c++)if(!z.piv.includes(c))free.push(c);return free.map(f=>{let v=Array(ncols).fill(0);v[f]=1;z.piv.forEach((c,i)=>v[c]=mod(-z.rows[i][f],p));return v;});}
function rank(A,p,ncols){return rref(A,p,ncols).rows.length;}
// Verify U relation modulo 2 and 3 as an independent regression check.
for (const p of [2,3]) {
  for(let i=0;i<N;i++)for(let j=0;j<N;j++){
    let s=mod(U[i][j]+(i===j?1:0),p);for(let k=0;k<N;k++)s=mod(s+U[i][k]*U[k][j],p);
    if(s!==0)throw new Error(`U relation failed mod ${p}`);
  }
}
// index 3 left eigenspace
let A3=[];for(let i=0;i<N;i++){let row=[];for(let j=0;j<N;j++)row.push(mod(U[j][i]-(i===j?1:0),3));A3.push(row);}let E=nullspace(A3,3);if(E.length!==6)throw new Error('bad eigenspace');console.log('index-3 projective kernels:',(3**E.length-1)/2);
// index 4 invariant codim2 kernels
let spaces=new Map();for(let num=1;num<(1<<N);num++){let v=Array(N).fill(0).map((_,i)=>(num>>i)&1);let vu=rowMul(v,U,2);let R=rref([v,vu],2).rows;if(R.length!==2)continue;let key=JSON.stringify(R);spaces.set(key,R);}console.log('index-4 invariant kernels:',spaces.size);
let counts={'(1,1,1,1,2,2)':0,'(1,1,1,1,1,4)':0};for(const R of spaces.values()){let K=nullspace(R,2);let H=Array.from({length:10},()=>Array(10).fill(0));for(let a=0;a<10;a++)for(let b=0;b<10;b++){let s=0;for(let i=0;i<N;i++)for(let j=0;j<N;j++)s^=(K[a][i]&mod(G[i][j],2)&K[b][j]);H[a][b]=s;}let rk=rank(H,2,10);if(rk===8)counts['(1,1,1,1,2,2)']++;else if(rk===10)counts['(1,1,1,1,1,4)']++;else throw new Error('unexpected rank '+rk);}console.log('index-4 type counts:',counts);
if(spaces.size!==1365||counts['(1,1,1,1,2,2)']!==693||counts['(1,1,1,1,1,4)']!==672)throw new Error('count mismatch');
console.log('independent Node enumeration: PASS');
