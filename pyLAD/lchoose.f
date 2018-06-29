      function lchoose(n, k)
      real*8 n, k
c
Cf2py intent(in) n
Cf2py intent(in) k
Cf2py intent(out) l
c
      real*8 lchoose
      lchoose = lgamma(n + 1) - lgamma(n - k + 1) - lgamma(k + 1)
      return
      end