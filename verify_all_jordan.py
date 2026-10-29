"""Exact full-determinant checks for every nonzero nilpotent type through N=10.

All entries and computations are Python integers. Complete tensor grids in
lambda (degree <=N) and t (degree <=g) prove the two-variable identity for
each listed Jordan type. This finite check complements the all-N graph proof.
"""
from pathlib import Path
import json
import time
from verify_higher_rank import det


def partitions(n,cap=None):
    if n==0:
        yield ()
        return
    if cap is None:cap=n
    for first in range(min(n,cap),0,-1):
        for tail in partitions(n-first,first):yield (first,)+tail


def matrices(parts):
    n=sum(parts)
    g=len(parts)
    a=[[0]*n for _ in range(n)]
    b=[[0]*n for _ in range(n)]
    heads=[]
    ends=[]
    cursor=0
    for p in parts:
        heads.append(cursor)
        ends.append(cursor+p-1)
        for j in range(p-1):a[cursor+j][cursor+j+1]=1
        cursor+=p
    for i in range(g):b[ends[(i+1)%g]][heads[i]]=1
    chosen=next(i for i,p in enumerate(parts) if p>=2)
    shortcut_row=ends[(chosen+1)%g]
    b[shortcut_row][heads[chosen]+1]=1
    return a,b,heads,ends,shortcut_row


def run():
    started=time.perf_counter()
    records=[]
    total=0
    for n in range(2,11):
        for parts in partitions(n):
            g=len(parts)
            if g==n:continue
            a,b,heads,ends,shortcut=matrices(parts)
            rank_minor=det([[b[i][j] for j in heads] for i in ends])
            assert abs(rank_minor)==1
            gram=[[sum(b[i][k]*b[j][k] for k in range(n))
                   for j in ends] for i in ends]
            for i in range(g):
                for j in range(g):
                    assert gram[i][j]==((2 if ends[i]==shortcut else 1) if i==j else 0)
            checked=0
            for lam in range(n+1):
                for t in range(g+1):
                    full=[[(lam if i==j else 0)-a[i][j]-t*b[i][j]
                           for j in range(n)] for i in range(n)]
                    actual=det(full)
                    expected=lam**n-t**g*(1+lam)
                    assert actual==expected,(parts,lam,t,actual,expected)
                    checked+=1
            total+=checked
            records.append(dict(N=n,jordan_type=parts,profile_rank=g,
                                rank_minor=rank_minor,
                                full_integer_determinants=checked,
                                gram_diagonal=[gram[i][i] for i in range(g)],
                                degree_bounds=dict(lam=n,t=g),
                                interpolation_status='complete tensor grid; exact'))
    result=dict(status='passed',jordan_types=len(records),
                full_integer_determinants=total,cases=records,
                range='Every nonzero nilpotent Jordan type in dimensions 2 through 10',
                arithmetic='exact integer full N-by-N determinants; no low-rank determinant reduction',
                conclusion='P(lambda,t)=lambda^N−t^g(1+lambda), rank B=g, norm B=sqrt(2)',
                scope='Finite exact algebraic checks complement the all-dimension proof; no numerical monodromy substituted for the group proof.',
                elapsed_seconds=time.perf_counter()-started)
    Path(__file__).with_name('all_jordan_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))


if __name__=='__main__':run()
