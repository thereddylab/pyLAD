
      logical function is_infinite(A)
      real*8 A
c
Cf2py intent(in) A
Cf2py intent(out) is_infinite
c
      if (A .eq. A + 1) then
        is_infinite = .TRUE.
      else
        is_infinite = .FALSE.
      endif
      return
      end

      logical function is_nonint(A)
      real*8 A
c
Cf2py intent(in) A
Cf2py intent(out) is_nonint
c
      if (abs(A) - floor(A + 0.5) .gt. 1e-7) then
        is_nonint = .TRUE.
      else
        is_nonint = .FALSE.
      endif
      return
      end

      logical function is_negnonint(A)
      real*8 A
c
Cf2py intent(in) A
Cf2py intent(out) is_negnonint
c
      logical is_nonint
      external is_nonint
      if ((A .lt. 0) .OR. (is_nonint(A))) then
        is_negnonint = .TRUE.
      else
        is_negnonint = .FALSE.
      endif
      return
      end

      real*8 function pdhyper (x, NR, NB, n, log_p)
      real*8 x, NR, NB, n
      logical log_p
c
Cf2py intent(in) x
Cf2py intent(in) NR
Cf2py intent(in) NB
Cf2py intent(in) n
Cf2py intent(in) log_p
Cf2py intent(out) pdhyper
c
      real*8 sum, term, x_

      sum = 0
      term = 1
      x_ = x
      do 10 while ((x_ .gt. 0) .AND. (term .ge. EPSILON(term) * sum))
        term = term * (x_ * (NB - n + x_) / (n + 1 - x) / (NR + 1 - x_))
        sum = sum + term
        x_ = x_ - 1
10    continue
      sum = sum + 1
      if (log_p) then
        pdhyper = log(sum)
      else
        pdhyper = sum
      endif
      return 
      end

      real*8 function phyper(x, NR, NB, n, lower_tail, log_p)
      real*8 x, NR, NB, n
      logical lower_tail, log_p
c
Cf2py intent(in) x
Cf2py intent(in) NR
Cf2py intent(in) NB
Cf2py intent(in) n
Cf2py intent(in) lower_tail
Cf2py intent(in) log_p
Cf2py intent(out) phyper
c
      real*8 x_i, NR_i, NB_i, n_i, old_NB, dh, pdh, dhyper, pdhyper,
     &  nan64
      logical is_infinite, lower_tail_i
      external dhyper, pdhyper, is_infinite

      nan64 = transfer(Z'7ff8000000000000', nan64)

      if (isnan(x) .OR. isnan(NR) .OR. isnan(NB) .OR. isnan(n)) then
        phyper = nan64
        go to 10
      endif
      
      x_i = FLOOR(x + 1e-7)
      NR_i = NINT(NR)
      NB_i = NINT(NB)
      n_i = NINT(n)
      lower_tail_i = lower_tail

      
      if ((NR .lt. 0) .OR. (NB .lt. 0) .OR. is_infinite(NR+NB) .OR.
     &    (n .lt. 0) .OR. (n_i .gt. NR_i+NB_i)) then
        phyper = nan64
        go to 10
      endif

      if (x_i*(NR_i+NB_i) .gt. n_i*NR_i) then
        old_NB = NB_i
        NB_i = NR_i
        NR_i = old_NB
        x_i = n_i - x_i - 1
        lower_tail_i = .NOT. lower_tail
      endif

      if (x .lt. 0) then
        phyper = 0
        go to 10
      endif
      if ((x_i .ge. NR_i) .OR. (x_i .ge. n_i)) then
        phyper = 1
      endif

      dh = dhyper(x_i, NR_i, NB_i, n_i, log_p)
      pdh = pdhyper(x_i, NR_i, NB_i, n_i, log_p)

      if (lower_tail_i) then
        if (log_p) then
          phyper = dh + pdh
        else
          phyper = dh * pdh
        endif
      else
        if (log_p) then
          phyper = log(1 - exp(dh + pdh))
        else
          phyper = 0.5 - (dh * pdh) + 0.5
        endif
      endif
10    return
      end

      real*8 function dhyper(x, r, b, n, give_log)
      real*8 x, r, b, n
      logical give_log
c
Cf2py intent(in) x
Cf2py intent(in) r
Cf2py intent(in) b
Cf2py intent(in) n
Cf2py intent(in) give_log
Cf2py intent(out) dhyper
c
      real*8 p, q, p1, p2, p3, x_i, r_i, b_i, n_i, dbinom_raw, nan64,
     &  infinity
      logical is_negnonint
      external is_negnonint, dbinom_raw

      infinity = transfer(Z'7ff0000000000000', infinity)
      nan64 = transfer(Z'7ff8000000000000', nan64)

      if (isnan(x) .OR. isnan(r) .OR. isnan(b) .OR. isnan(n)) then
        dhyper = nan64
        go to 10
      endif

      if (is_negnonint(r) .OR. is_negnonint(b) .OR. is_negnonint(n)
     &   .OR. (n .gt. r+b)) then
        dhyper = nan64
        go to 10
      endif
      if (x .lt. 0) then
        if (give_log) then
          dhyper = -infinity
        else
          dhyper = 0
        endif
        go to 10
      endif

      x_i = REAL(NINT(x), 8)
      r_i = REAL(NINT(r), 8)
      b_i = REAL(NINT(b), 8)
      n_i = REAL(NINT(n), 8)

      if ((n_i .lt. x_i) .OR. (r_i .lt. x_i) .OR. (n-x .gt. b)) then
        dhyper = 0
        go to 10
      endif
      if (n_i .eq. 0) then
        if (x .eq. 0) then
          if (give_log) then
            dhyper = 0.
          else
            dhyper = 1.
          endif
        else
          if (give_log) then
            dhyper = -infinity
          else
            dhyper = 0.
          endif
        endif
        go to 10
      endif

      p = n_i / (r_i+b_i)
      q = (r_i+b_i-n_i) / (r_i+b_i)

      p1 = dbinom_raw(x_i, r_i, p, q, give_log)
      p2 = dbinom_raw(n_i-x_i, b_i, p, q, give_log)
      p3 = dbinom_raw(n_i, r_i+b_i, p, q, give_log)

      if (give_log) then
        dhyper = p1 + p2 - p3
      else
        dhyper = p1 * p2 / p3
      endif
10    return
      end


      real*8 function dbinom_raw(x, n, p, q, give_log)
      real*8 x, n, p, q
      logical give_log
c
Cf2py intent(in) x
Cf2py intent(in) n
Cf2py intent(in) p
Cf2py intent(in) q
Cf2py intent(in) give_log
Cf2py intent(out) dbinom_raw
c
      real*8, PARAMETER :: LN_2PI = 1.837877066409345483560659472811
      real*8 lf, lc, stirlerr, bd0, infinity
      external stirlerr, bd0

      infinity = transfer(Z'7ff0000000000000', infinity)

      if (p .eq. 0) then
        if (x .eq. 0) then
          if (give_log) then
            dbinom_raw = 0.
          else
            dbinom_raw = 1.
          endif
        else
          if (give_log) then
            dbinom_raw = -infinity
          else
            dbinom_raw = 0.
          endif
        endif
        go to 10
      endif
      if (q .eq. 0) then
        if (x .eq. n) then
          if (give_log) then
            dbinom_raw = 0.
          else
            dbinom_raw = 1.
          endif
        else
          if (give_log) then
            dbinom_raw = -infinity
          else
            dbinom_raw = 0.
          endif
        endif
        go to 10
      endif

      if (x .eq. 0) then
        if (n .eq. 0) then
          if (give_log) then
            dbinom_raw = 0.
          else
            dbinom_raw = 1.
          endif
          go to 10
        endif
        if (p .lt. 0.1) then
          lc = -bd0(n, n*q) - n*p
        else
          lc = n*log(q)
        endif
        if (give_log) then
          dbinom_raw = lc
        else
          dbinom_raw = exp(lc)
        endif
        go to 10  
      endif

      if (x .eq. n) then
        if (q .lt. 0.1) then
          lc = -bd0(n, n*p) - n*q
        else
          lc = n*log(p)
        endif
        if (give_log) then
          dbinom_raw = lc
        else
          dbinom_raw = exp(lc)
        endif
        go to 10
      endif
      if ((x .lt. 0) .OR. (x .gt. n)) then
        if (give_log) then
          dbinom_raw = -infinity
        else
          dbinom_raw = 0.
        endif
        go to 10
      endif

c n*p or n*q can underflow to zero if n and p or q are small.  This
c used to occur in dbeta, and gives NaN as from R 2.3.0.
      lc = stirlerr(n)-stirlerr(x)-stirlerr(n-x)-bd0(x,n*p)-bd0(n-x,n*q)

c f = (M_2PI*x*(n-x))/n; could overflow or underflow
c Upto R 2.7.1:
c lf = log(M_2PI) + log(x) + log(n-x) - log(n);
c -- following is much better for  x << n :
      lf = LN_2PI + log(x) + log(1 - x/n)
      if (give_log) then
        dbinom_raw = lc - 0.5 * lf
      else
        dbinom_raw = exp(lc - 0.5 * lf)
      endif
10    return
      end

      real*8 function stirlerr(n)
      real*8 n
c
Cf2py intent(in) n
Cf2py intent(out) stirlerr
c
      real*8 nn
c 1/12
      real*8, PARAMETER :: S0 = 0.083333333333333333333
c 1/360
      real*8, PARAMETER :: S1 = 0.00277777777777777777778
c 1/1260
      real*8, PARAMETER :: S2 = 0.00079365079365079365079365
c 1/1680
      real*8, PARAMETER :: S3 = 0.000595238095238095238095238
c 1/1188
      real*8, PARAMETER :: S4 = 0.0008417508417508417508417508
      real*8, PARAMETER :: LN_SQRT_2PI = 
     &  0.918938533204672741780329736406
c error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
      real*8, dimension(0:30) :: sferr_halves
      sferr_halves(0:30) = (/ 0.0,
     &  0.1534264097200273452913848,
     &  0.0810614667953272582196702,
     &  0.0548141210519176538961390,
     &  0.0413406959554092940938221,
     &  0.03316287351993628748511048,
     &  0.02767792568499833914878929,
     &  0.02374616365629749597132920,
     &  0.02079067210376509311152277,
     &  0.01848845053267318523077934,
     &  0.01664469118982119216319487,
     &  0.01513497322191737887351255,
     &  0.01387612882307074799874573,
     &  0.01281046524292022692424986,
     &  0.01189670994589177009505572,
     &  0.01110455975820691732662991,
     &  0.010411265261972096497478567,
     &  0.009799416126158803298389475,
     &  0.009255462182712732917728637,
     &  0.008768700134139385462952823,
     &  0.008330563433362871256469318,
     &  0.007934114564314020547248100,
     &  0.007573675487951840794972024,
     &  0.007244554301320383179543912,
     &  0.006942840107209529865664152,
     &  0.006665247032707682442354394,
     &  0.006408994188004207068439631,
     &  0.006171712263039457647532867,
     &  0.005951370112758847735624416,
     &  0.005746216513010115682023589,
     &  0.005554733551962801371038690/)

      if (n .le. 15.0) then
        nn = n + n
        if (nn .eq. INT(nn)) then
          stirlerr = sferr_halves (INT(nn))
        else
          stirlerr = LOG_GAMMA(n+1) - (n+0.5) * log(n) + n - LN_SQRT_2PI
        endif
        go to 10
      endif

      nn = n*n
      if (n .gt. 500) then
        stirlerr = (S0 - S1/nn) / n
      else if (n .gt. 80) then
        stirlerr = (S0 - (S1 - S2/nn) / nn) / n
      else if (n .gt. 35) then
        stirlerr = (S0 - (S1 - (S2 - S3/nn) / nn) / nn) / n
      else
        stirlerr = (S0 - (S1 - (S2 - (S3 - S4/nn) / nn) / nn) / nn) / n
      endif
10    return
      end

      real*8 function bd0(x, np)
      real*8 x, np
c
Cf2py intent(in) x
Cf2py intent(in) np
Cf2py intent(out) bd0
c
      real*8 ej, s, s1, v, nan64
      integer j
      logical is_infinite
      external is_infinite

      nan64 = transfer(Z'7ff8000000000000', nan64)

      if (is_infinite(x) .OR. is_infinite(np) .OR. (np .eq. 0)) then
        bd0 = nan64
        go to 100
      endif

      if (abs(x-np) .lt. 0.1 * (x+np)) then
        v = (x - np) / (x + np)
        s = (x - np) * v
        if (abs(s) .lt. TINY(s)) then
          bd0 = s
          go to 100
        endif
        ej = 2 * x * v
        v = v * v
        do 10 j = 1, 1000
          ej = ej * v
          s1 = s + ej / (SHIFTL(j, 1) + 1)
          if (s1 .eq. s) then
            bd0 = s1
            go to 100
          endif
          s = s1
10      continue
      endif
      bd0 = x * log(x / np) + np - x
100   return
      end
