C=== This was 'sort()' in  gamfit's  mysort.f  [or sortdi() in sortdi.f ] :
C
C     a[] is real*8 because the caller reuses a double (sometimes v[] itself!)

      subroutine sort (v,ii,jj,n)
c
Cf2py intent(in,out) v(n)
Cf2py intent(in) ii
Cf2py intent(in) jj
Cf2py integer intent(hide),depend(v) :: n=shape(v, 0)
c
c     Puts into a the permutation vector which sorts v into
c     increasing order.  Only elements from ii to jj are considered.
c     Arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c
c     This is a modification of CACM algorithm #347 by R. C. Singleton,
c     which is a modified Hoare quicksort.
c
      integer ii,jj,n
      real*8 v(n)
c
      integer iu(20),il(20)
      integer m,i,j,ij,k,l
      real*8 vt, vtt

      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      vt=v(ij)
      if (v(i).le.vt) go to 30
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      v(k+1)=vt
      go to 100
      end

C=== This was 'sort()' in  gamfit's  mysort.f  [or sortdi() in sortdi.f ] :
C
C     a[] is real*8 because the caller reuses a double (sometimes v[] itself!)
      subroutine sort_i (v,a,ii,jj,n)
c
Cf2py intent(in) v(n)
Cf2py intent(in,out) a(n)
Cf2py intent(in) ii
Cf2py intent(in) jj
Cf2py integer intent(hide),depend(v) :: n=shape(v, 0)
c
c     Puts into a the permutation vector which sorts v into
c     increasing order.  Only elements from ii to jj are considered.
c     Arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c
c     This is a modification of CACM algorithm #347 by R. C. Singleton,
c     which is a modified Hoare quicksort.
c
      integer ii,jj,n,a(n)
      real*8 v(n)
c
      integer iu(20),il(20)
      integer t,tt, m,i,j,ij,k,l
      real*8 vt, vtt

      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=int(a(ij))
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=a(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      a(k+1)=t
      v(k+1)=vt
      go to 100
      end
