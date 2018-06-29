c     tail probability of circular binary segmentation statistic
c     from Siegmund (1988) or Yao (1989) paper
      double precision function tailp(b, delta, m, ngrid, tol)
      double precision b, delta, tol
      integer m, ngrid
c
Cf2py intent(in) b
Cf2py intent(in) delta
Cf2py intent(in) m
Cf2py intent(in) ngrid
Cf2py intent(in) tol
Cf2py intent(out) tailp
c
c     it1tsq is the integral of 1/(t*(1-t))**2
      double precision nu, it1tsq
      external nu, it1tsq

      double precision t, tl, dincr, bsqrtm, x, nux
      integer i

      dincr = (0.5d0 - delta)/dfloat(ngrid)
      bsqrtm = b/sqrt(dfloat(m))

      tl = 0.5d0 - dincr
      t = 0.5d0 - 0.5d0*dincr
      tailp = 0.0d0
      do 10 i = 1,ngrid
         tl = tl + dincr
         t = t + dincr
         x = bsqrtm/sqrt(t*(1-t))
         nux = nu(x, tol)
         tailp = tailp + (nux**2)*it1tsq(tl, dincr)
 10   continue
      tailp = 9.973557d-2*(b**3)*exp(-b**2/2)*tailp
c     since test is two-sided need to multiply tailp by 2
      tailp = 2.0d0*tailp

      return
      end

c     integral of 1/(t*(1-t))**2 from x to x+a
      double precision function it1tsq(x, a)
      double precision x, a
c
Cf2py intent(in) x
Cf2py intent(in) a
Cf2py intent(out) it1tsq
c
      double precision y

      y = x + a - 0.5d0
      it1tsq = (8.0d0*y)/(1.0d0 - 4.0d0*y**2) + 
     1     2.0d0*log((1.0d0 + 2.0d0*y)/(1.0d0 - 2.0d0*y))
      y = x - 0.5d0
      it1tsq = it1tsq - (8.0d0*y)/(1.0d0 - 4.0d0*y**2) -
     1     2.0d0*log((1.0d0 + 2.0d0*y)/(1.0d0 - 2.0d0*y))

      return
      end

      double precision function nu(x, tol)
      double precision x, tol
c
Cf2py intent(in) x
Cf2py intent(in) tol
Cf2py intent(out) nu
c
      double precision alnorm
      external alnorm

      double precision lnu0, lnu1, dk, xk
      integer i, k

      if (x .gt. 0.01d0) then
         lnu1 = log(2.0d0) - 2*log(x)
         lnu0 = lnu1
         k = 2
         dk = 0
         do 10 i = 1, k
            dk = dk + 1
            xk = -x*sqrt(dk)/2.0d0
            lnu1 = lnu1 - 2.0d0*alnorm(xk, 0)/dk
 10      continue

         do 50 while (dabs((lnu1-lnu0)/lnu1) .gt. tol)
            lnu0 = lnu1
            do 20 i = 1,k
               dk = dk + 1
               xk = -x*sqrt(dk)/2.0d0
               lnu1 = lnu1 - 2.0d0*alnorm(xk, 0)/dk
 20         continue
            k = 2*k
 50      enddo
      else
         lnu1 = -0.583d0*x
      endif
      nu = exp(lnu1)

      return
      end

c     tail probability of binary segmentation statistic
c     from page 387 of Siegmund (1986) paper
      double precision function btailp(b, m, ng, tol)
      integer m, ng
      double precision b, tol
c
Cf2py intent(in) b
Cf2py intent(in) m
Cf2py intent(in) ng
Cf2py intent(in) tol
Cf2py intent(out) btailp
c
      double precision ll, ul, dincr, nulo, nuhi, x, x1, dm
      integer i, k

      double precision alnorm, nu
      external alnorm, nu

      dm = dfloat(m)
      k = 2
      ll = b*sqrt(1.0/dfloat(m-k) - 1.0/dfloat(m))
      ul = b*sqrt(1.0/dfloat(k) - 1.0/dfloat(m))
      dincr = (ul - ll)/dfloat(ng)

      btailp = 0.0
      x = ll
      x1 = x + (b**2)/(dm*x)
      nulo = nu(x1, tol)/x
      do 10 i = 1, ng
         x = x + dincr
         x1 = x + (b**2)/(dm*x)
         nuhi = nu(x1, tol)/x
         btailp = btailp + (nuhi + nulo)*dincr
         nulo = nuhi
 10   continue
      btailp = b*exp(-b**2/2)*btailp/2.506628275

      btailp =  btailp + 2*(1.0-alnorm(b, 0))

      return
      end
