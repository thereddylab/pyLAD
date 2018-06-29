      subroutine make_KMP_table(T, W, N)
      integer T(0:(N-1)), N
      character W(0:(N-1))
c
Cf2py intent(in,out) T
Cf2py intent(in) W
Cf2py integer intent(hide),depend(T) :: N=shape(T, 0)
c
      integer pos, cnd
      T(0) = -1
      pos = 1
      cnd = 0
      do 10 while (pos .lt. N)
        if (W(pos) .eq. W(cnd)) then
          T(pos) = T(cnd)
          pos = pos + 1
          cnd = cnd + 1
        else
          T(pos) = cnd
          cnd = T(cnd)
          do 20 while ((cnd .ge. 0) .AND. (W(pos) .ne. W(cnd)))
            cnd = T(cnd)
20        continue
          pos = pos + 1
          cnd = cnd + 1
        endif
10    continue
      return
      end

      integer function overlap(W, S, T, N, M)
      character W(0:(N-1)), S(0:(M-1))
      integer T(0:(N-1)), N, M
c
Cf2py intent(in) W
Cf2py intent(in) S
Cf2py intent(in) T
Cf2py integer intent(hide),depend(W) :: N=shape(W,0)
Cf2py integer intent(hide),depend(S) :: M=shape(S,0)
c
      integer j, k, l
      j = 0
      k = 0
      if (N .lt. M) then
        l = N
      else
        l = M
      endif
      do 10 while (j .lt. l)
        if (W(k) .eq. S(j)) then
          j = j + 1
          k = k + 1
        else
          k = T(k)
          if (k .lt. 0) then
            j = j + 1
            k = k + 1
          endif
        endif
10    continue
      overlap = k
      return
      end