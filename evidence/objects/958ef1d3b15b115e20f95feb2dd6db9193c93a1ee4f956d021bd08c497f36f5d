// Independent Node.js check of the A-vertex CM-character and divisor algebra.
const gcd=(a,b)=>b?gcd(b,a%b):Math.abs(a);
const G=[...Array(42).keys()].filter(a=>a>0 && gcd(a,42)===1);
const Phi=new Set([1,5,31,29,17,19]);
const key=S=>[...S].sort((a,b)=>a-b).join(",");
function combos(arr,k,start=0,prefix=[],out=[]){
  if(prefix.length===k){out.push([...prefix]); return out;}
  for(let i=start;i<=arr.length-(k-prefix.length);i++){
    prefix.push(arr[i]); combos(arr,k,i+1,prefix,out); prefix.pop();
  }
  return out;
}
function isHodge(S,p){
  for(const g of G){
    let n=0;
    for(const s of S) if(Phi.has((g*s)%42)) n++;
    if(n!==p) return false;
  }
  return true;
}
const H={};
for(let p=1;p<=6;p++){
  H[p]=new Map();
  for(const S of combos(G,2*p)) if(isHodge(S,p)) H[p].set(key(S),new Set(S));
}
const counts=[1,2,3,4,5,6].map(p=>H[p].size);
if(JSON.stringify(counts)!==JSON.stringify([12,51,88,51,12,1])) throw Error("Hodge counts");
const divisors=[...H[1].values()];
const gen=[];
for(let p=1;p<=6;p++){
  const U=new Set();
  for(const C of combos(divisors,p)){
    const s=new Set();
    for(const D of C) for(const x of D) s.add(x);
    if(s.size===2*p) U.add(key(s));
  }
  if(U.size!==H[p].size) throw Error(`generated count p=${p}`);
  for(const k of H[p].keys()) if(!U.has(k)) throw Error(`missing p=${p} ${k}`);
  gen.push(U.size);
}
const stab=G.filter(g=>{
  const S=new Set([...Phi].map(x=>(g*x)%42));
  return key(S)===key(Phi);
});
if(JSON.stringify(stab)!==JSON.stringify([1,29])) throw Error("stabilizer");
console.log("Phi_A =",key(Phi));
console.log("Hodge counts =",counts.join(","));
console.log("divisor-product counts =",gen.join(","));
console.log("stabilizer =",stab.join(","));
console.log("ALL NODE A-VERTEX HODGE ASSERTIONS PASSED");
