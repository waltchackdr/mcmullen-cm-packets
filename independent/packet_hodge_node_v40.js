#!/usr/bin/env node
// Independent combinatorial recheck of the twelve-class CM-type orbit/Hodge profile.
const gcd=(a,b)=>b?gcd(b,a%b):Math.abs(a);
const G=[...Array(41)].map((_,i)=>i+1).filter(a=>gcd(a,42)===1);
const base=[1,5,11,13,17,19];
const indices=[1,11,16,17,19,20,22,32,37,38,40,41];
const mod=(x,n=42)=>((x%n)+n)%n;
const key=S=>[...S].sort((a,b)=>a-b).join(',');
function cmType(k){
  const out=[];
  for(const r of base){
    const th=2*Math.PI*r/42;
    let im=0;
    for(const e of [k,k+1,k+4]) im+=Math.sin(e*th);
    if(Math.abs(im)<1e-10) throw Error('ambiguous sign');
    out.push(im>0?r:mod(-r));
  }
  return new Set(out);
}
const types=Object.fromEntries(indices.map(k=>[k,cmType(k)]));
function orbit(S){return new Set(G.map(g=>key(new Set([...S].map(a=>mod(g*a))))));}
let remaining=new Set(indices), blocks=[];
while(remaining.size){
  const k=Math.min(...remaining), o=orbit(types[k]);
  const b=indices.filter(j=>o.has(key(types[j])));
  blocks.push(b); for(const j of b) remaining.delete(j);
}
const expected=JSON.stringify([[1,20,22,41],[11,16,17,19,32,37,38,40]]);
if(JSON.stringify(blocks)!==expected) throw Error('orbit partition mismatch '+JSON.stringify(blocks));
function combinations(arr,k,start=0,prefix=[],out=[]){
 if(prefix.length===k){out.push(prefix.slice());return out;}
 for(let i=start;i<=arr.length-(k-prefix.length);i++){prefix.push(arr[i]);combinations(arr,k,i+1,prefix,out);prefix.pop();}
 return out;
}
function hodgeSets(PH,p){
 const P=new Set(PH), out=new Set();
 for(const S of combinations(G,2*p)){
   let ok=true;
   for(const g of G){
     let n=0; for(const s of S) if(P.has(mod(g*s))) n++;
     if(n!==p){ok=false;break;}
   }
   if(ok) out.add(key(new Set(S)));
 }
 return out;
}
for(const block of blocks){
 const S=types[block[0]];
 const stab=G.filter(g=>key(new Set([...S].map(a=>mod(g*a))))===key(S));
 const counts=[]; for(let p=1;p<=6;p++) counts.push(hodgeSets(S,p).size);
 const es=block[0]===1?[12,51,88,51,12,1]:[6,15,22,15,6,1];
 if(JSON.stringify(counts)!==JSON.stringify(es)) throw Error('count mismatch');
 if(block[0]===1 && JSON.stringify(stab)!==JSON.stringify([1,29])) throw Error('imprimitive stabilizer');
 if(block[0]!==1 && JSON.stringify(stab)!==JSON.stringify([1])) throw Error('primitive stabilizer');
 if(block[0]!==1){
   const pairs=[]; const seen=new Set();
   for(const a of G){const q=key(new Set([a,mod(-a)])); if(!seen.has(q)){seen.add(q);pairs.push(q.split(',').map(Number));}}
   const div3=new Set(combinations(pairs,3).map(cc=>key(new Set(cc.flat()))));
   const H=hodgeSets(S,3); const rem=[...H].filter(x=>!div3.has(x));
   if(rem.length!==2) throw Error('expected exactly two Weil lines');
 }
 console.log('block',block.join(','),'stabilizer',stab.join(','),'counts',counts.join(','));
}
console.log('PACKET HODGE NODE v40: PASS');
