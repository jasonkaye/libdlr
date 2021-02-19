c
c
c       dependencies: prini, dfft, and (for the debugging code) id_rand
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer n,k,l,ind(len),nblock
        real*8 r1,rnd,diffmax
        complex*16 v(len),v2(len),disp(len),disp2(len),
     1             wsave(len),w(len)
c
        r1 = 1
c
c
        call prini(6,13)
c
c
        print *,'Enter n:'
        read *,n
        call prinf('n = *',n,1)
c
c
        l = 15
        call prinf('l = *',l,1)
c
c
c       Choose at random locations to calculate.
c
        do k = 1,l
          call id_srand(1,rnd)
          ind(k) = 1+n*rnd
          if(ind(k) .gt. n) ind(k) = n
        enddo ! k
        call prinf('ind = *',ind,l)
c
c
c       Fill the real and imaginary parts of every entry
c       of the vector v being transformed with i.i.d. random variates,
c       drawn uniformly from [0,1].
c
        call id_srand(2*n,v)
c
c
c       Copy v into v2.
c
        do k = 1,n
          v2(k) = v(k)
        enddo ! k
c
c
c       Transform v2 via zfft2.
c
        call idz_ldiv(l,n,nblock)
        call zfftf2(nblock,n,v2,w)
c
        do k = 1,l
          disp2(k) = v2(ind(k))
        enddo ! k
c
        call prin2('disp2 = *',disp2,2*l)
c
c
c       Transform v via idz_sfft.
c
        call idz_sffti(l,ind,n,wsave)
        call idz_sfft(l,ind,n,wsave,v)
c
        do k = 1,l
          disp(k) = v(ind(k))
        enddo ! k
c
        call prin2('disp = *',disp,2*l)
c
c
c       Find the difference between disp and disp2.
c
        diffmax = 0
c
        do k = 1,l
          if(abs(disp(k)-disp2(k)) .gt. diffmax)
     1     diffmax = abs(disp(k)-disp2(k))
        enddo ! k
c
        call prin2('diffmax = *',diffmax,1)
c
c
        stop
        end
c
c
c
c
        subroutine zfftf2(l,n,v,w)
c
c       calculate via a full FFT what idz_sfft would compute
c       if it calculated every entry, rather than just a subset
c       of l entries.
c
c       input:
c       l -- number of entries in the output to compute
c       n -- length of v
c       v -- vector to be transformed
c
c       output:
c       v -- transformed vector
c
c       work:
c       w -- must be at least 2*n+15 complex*16 elements long
c
        implicit none
        integer n,l,k
        real*8 r1,fact
        complex*16 v(n),w(2*n+15)
c
        r1 = 1
c
c
c       Transpose v.
c
        call vectrans(l,n,v,w)
c
c
c       Transform v via zfftf and normalize it.
c
        call zffti(n,w)
        call zfftf(n,v,w)
c
        fact = 1/sqrt(r1*n)
        do k = 1,n
          v(k) = v(k)*fact
        enddo ! k
c
c
c       Transpose v.
c
        call vectrans(l,n,v,w)
c
c
        return
        end
c
c
c
c
        subroutine vectrans(l,n,v,w)
c
c       transposes v.
c
c       input:
c       l -- first stride
c       n -- length of v
c
c       output:
c       v -- transposed vector
c
c       work:
c       w -- must be at least n complex*16 elements long
c
        implicit none
        integer n,m,l,k,j
        complex*16 v(n),w(n)
c
c
        m = n/l
c
c
c       Transpose v to obtain w.
c
        do k = 1,m
          do j = 1,l
            w(m*(j-1)+k) = v(l*(k-1)+j)
          enddo ! j
        enddo ! k
c
c
c       Copy w into v.
c
        do k = 1,n
          v(k) = w(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
