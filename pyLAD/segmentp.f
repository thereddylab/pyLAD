      subroutine bsegp(n, gendat, ostat, pval, ng, tol)
      integer n, ng
      real*8 gendat(n), ostat, pval, tol
c
Cf2py integer intent(in),depend(gendat) :: n=shape(gendat,0)
Cf2py intent(in) gendat
Cf2py intent(in,out) ostat
Cf2py intent(in,out) pval
Cf2py intent(in) ng
Cf2py intent(in) tol
c
      real*8 btmax, btailp
      external btmax, btailp

      ostat = btmax(n, gendat)
c      call dblepr("Max Stat",8,ostat,1)
      pval = btailp(ostat, n, ng, tol)
      if (pval .gt. 1) pval = 1.0d0

      return
      end

      real*8 function btmax(n, x)
      integer n
      real*8 x(n)
c
Cf2py integer intent(in),depend(x) :: n=shape(x,0)
Cf2py intent(in) x
Cf2py intent(out) btmax
c
      integer i
      real*8 sumxi, btmaxi, dn, di, ostat

      sumxi = x(1)
      ostat = 0.0
      dn = dfloat(n)
      di = 1.0
      do 20 i = 2,n-2
         di = di + 1.0
         sumxi = sumxi + x(i)
         btmaxi = dn*(sumxi**2)/(di*(dn-di))
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
c            ibseg = i
         endif
 20   continue
      btmax = sqrt(ostat)

      return
      end

c     pseudo confidence interval based on permutations
      subroutine bsegci(n, k, sumxk, x, px, sr, vfact, nperm, bsloc)
      integer n, k, sr(2), nperm, bsloc(nperm)
      real*8 sumxk, x(n), px(n), vfact(n)
c
Cf2py integer intent(in),depend(x) :: n=shape(x,0)
Cf2py intent(in) k
Cf2py intent(in) sumxk
Cf2py intent(in) x
Cf2py real*8 intent(in),depend(n) :: px
Cf2py intent(in) sr
Cf2py real*8 intent(in),depend(n) :: vfact
Cf2py integer intent(in),depend(bsloc) :: nperm=shape(bsloc,0)
Cf2py intent(in,out) bsloc
c
      integer k1, nk, np, ibseg

      call random_seed()
      k1 = k+1
      nk = n-k
      do 10 np = 1, nperm
         call xperm(k,x,px)
         call xperm(nk,x(k1),px(k1))
         call btmxci(n,k,sr,px,vfact,ibseg,sumxk)
         bsloc(np) = ibseg
 10   continue

      return
      end

      subroutine btmxci(n,k,sr,x,vfact,ibseg,sumxk)
      integer n,k,sr(2),ibseg
      real*8 x(n),vfact(n),sumxk
c
Cf2py integer intent(in),depend(x) :: n=shape(x,0)
Cf2py intent(in) k
Cf2py intent(in) sr
Cf2py intent(in) x
Cf2py real*8 intent(in),depend(n) :: vfact
Cf2py intent(in,out) ibseg
Cf2py intent(in) sumxk
c
      integer i
      real*8 sumxi, ostat, btmaxi

      ostat = vfact(k)*(sumxk**2)
      ibseg = k
      sumxi = sumxk
      do 10 i = k-1,sr(1),-1
         sumxi = sumxi - x(i+1)
         btmaxi = vfact(i)*(sumxi**2)
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            ibseg = i
         endif
 10   continue

      sumxi = sumxk
      do 20 i = k+1,sr(2)
         sumxi = sumxi + x(i)
         btmaxi = vfact(i)*(sumxi**2)
         if (ostat .lt. btmaxi) then
            ostat = btmaxi
            ibseg = i
         endif
 20   continue
      ostat = sqrt(ostat)

      return
      end
