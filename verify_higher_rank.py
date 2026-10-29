"""Independent integer determinant checks of arbitrary-block fusion.

No symbolic-algebra or numerical package is required. The complete tensor
interpolation grids establish the displayed characteristic identity for each
listed (m,r), with arbitrary complex lambda,t,delta: both sides have degrees
at most (mr,r,1). They do not replace the proof for arbitrary m,r.
"""
from pathlib import Path
import json
import time


def det(a):
    """Fraction-free Bareiss elimination with exact integer pivoting."""
    a=[list(row) for row in a]
    n=len(a)
    if not n:return 1
    sign=1
    previous=1
    for k in range(n-1):
        pivot_row=next((i for i in range(k,n) if a[i][k]),None)
        if pivot_row is None:return 0
        if pivot_row!=k:
            a[k],a[pivot_row]=a[pivot_row],a[k]
            sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=a[i][j]*pivot-a[i][k]*a[k][j]
                assert value%previous==0
                a[i][j]=value//previous
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]


def matrices(m,r,delta):
    n=m*r
    a=[[0]*n for _ in range(n)]
    b=[[0]*n for _ in range(n)]
    for i in range(r):
        for j in range(m-1):a[i*m+j][i*m+j+1]=1
        b[i*m+m-1][i*m]=1
        if i<r-1:b[i*m+m-1][(i+1)*m]=1
    b[n-1][0]=-delta
    b[n-1][1]=1
    return a,b


def run():
    started=time.perf_counter()
    cases=[]
    total=0
    for m in (2,3,4):
        for r in (2,3,4,5):
            n=m*r
            checked=0
            rank_minors=[]
            for delta in (0,2):
                a,b=matrices(m,r,delta)
                rows=[i*m+m-1 for i in range(r)]
                cols=[i*m for i in range(r-1)]+[1]
                minor=det([[b[i][j] for j in cols] for i in rows])
                assert minor==1
                rank_minors.append(minor)
                for lam in range(n+1):
                    for t in range(r+1):
                        matrix=[[(lam if i==j else 0)-a[i][j]-t*b[i][j]
                                 for j in range(n)] for i in range(n)]
                        actual=det(matrix)
                        expected=(lam**m-t)**r-t**r*(lam-delta)
                        assert actual==expected,(m,r,delta,lam,t,actual,expected)
                        checked+=1
            total+=checked
            cases.append(dict(m=m,r=r,N=n,integer_determinants=checked,
                              rank_minor_determinants=rank_minors,
                              polynomial_degree_bounds=dict(lam=n,t=r,delta=1),
                              interpolation_status='complete tensor grid; exact'))
    result=dict(status='passed',cases=cases,total_integer_determinants=total,
                parameter_scope='arbitrary complex lambda,t,delta for each listed (m,r)',
                arithmetic='exact Python integers; full mr-by-mr determinant, no low-rank formula',
                elapsed_seconds=time.perf_counter()-started)
    dest=Path(__file__).with_name('higher_rank_verification.json')
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
