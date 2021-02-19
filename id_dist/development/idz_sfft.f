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
c
c
c       routine idz_sffti initializes routine idz_sfft.
c
c       routine idz_sfft rapidly computes a subset of the entries
c       of the DFT of a vector, composed with permutation matrices
c       both on input and on output.
c
c       routine idz_ldiv finds the greatest integer less than or equal
c       to a specified integer, that is divisible by another (larger)
c       specified integer.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idz_ldiv(l,n,m)
c
c       finds the greatest integer less than or equal to l
c       that divides n.
c
c       input:
c       l -- integer at least as great as m
c       n -- integer divisible by m
c
c       output:
c       m -- greatest integer less than or equal to l that divides n
c
        implicit none
        integer n,l,m
c
c
        m = l
c
 1000   continue
        if(m*(n/m) .eq. n) goto 2000
c
          m = m-1
          goto 1000
c
 2000   continue
c
c
        return
        end
c
c
c
c
        subroutine idz_sffti(l,ind,n,wsave)
c
c       initializes wsave for use with routine idz_sfft.
c
c       input:
c       l -- number of entries in the output of idz_sfft to compute
c       ind -- indices of the entries in the output of idz_sfft
c              to compute
c       n -- length of the vector to be transformed
c
c       output:
c       wsave -- array needed by routine idz_sfft for processing
c
        implicit none
        integer l,ind(l),n,nblock,ii,m,idivm,imodm,i,j,k
        real*8 r1,twopi,fact
        complex*16 wsave(2*l+15+3*n),ci,twopii
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
        twopii = twopi*ci
c
c
c       Determine the block lengths for the FFTs.
c
        call idz_ldiv(l,n,nblock)
        m = n/nblock
c
c
c       Initialize wsave for use with routine zfftf.
c
        call zffti(nblock,wsave)
c
c
c       Calculate the coefficients in the linear combinations
c       needed for the direct portion of the calculation.
c
        fact = 1/sqrt(r1*n)
c
        ii = 2*l+15
c
        do j = 1,l
c
          i = ind(j)
c
          idivm = (i-1)/m
          imodm = (i-1)-m*idivm
c
          do k = 1,m
            wsave(ii+m*(j-1)+k) = exp(-twopii*imodm*(k-1)/(r1*m))
     1       * exp(-twopii*(k-1)*idivm/(r1*n)) * fact
          enddo ! k
c
        enddo ! j
c
c
        return
        end
c
c
c
c
        subroutine idz_sfft(l,ind,n,wsave,v)
c
c       computes a subset of the entries of the DFT of v,
c       composed with permutation matrices both on input and on output,
c       via a two-stage procedure (routine zfftf2 is supposed
c       to calculate the full vector from which idz_sfft returns
c       a subset of the entries, when zfftf2 has the same parameter
c       nblock as in the present routine).
c
c       input:
c       l -- number of entries in the output to compute
c       ind -- indices of the entries of the output to compute
c       n -- length of v
c       v -- vector to be transformed
c       wsave -- processing array initialized by routine idz_sffti
c
c       output:
c       v -- entries indexed by ind are given their appropriate
c            transformed values
c
c       _N.B._: The user has to boost the memory allocations
c               for wsave (and change iii accordingly) if s/he wishes
c               to use strange sizes of n; it's best to stick to powers
c               of 2.
c
c       references:
c       Sorensen and Burrus, "Efficient computation of the DFT with
c            only a subset of input or output points,"
c            IEEE Transactions on Signal Processing, 41 (3): 1184-1200,
c            1993.
c       Woolfe, Liberty, Rokhlin, Tygert, "A fast randomized algorithm
c            for the approximation of matrices," Applied and
c            Computational Harmonic Analysis, 25 (3): 335-366, 2008;
c            Section 3.3.
c
        implicit none
        integer n,m,l,k,j,ind(l),i,idivm,nblock,ii,iii
        real*8 r1,twopi
        complex*16 v(n),wsave(2*l+15+3*n),ci,sum
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
c       Determine the block lengths for the FFTs.
c
        call idz_ldiv(l,n,nblock)
c
c
        m = n/nblock
c
c
c       FFT each block of length nblock of v.
c
        do k = 1,m
          call zfftf(nblock,v(nblock*(k-1)+1),wsave)
        enddo ! k
c
c
c       Transpose v to obtain wsave(2*l+15+2*n+1 : 2*l+15+3*n).
c
        iii = 2*l+15+2*n
c
        do k = 1,m
          do j = 1,nblock
            wsave(iii+m*(j-1)+k) = v(nblock*(k-1)+j)
          enddo ! j
        enddo ! k
c
c
c       Directly calculate the desired entries of v.
c
        ii = 2*l+15
        iii = 2*l+15+2*n
c
        do j = 1,l
c
          i = ind(j)
c
          idivm = (i-1)/m
c
          sum = 0
c
          do k = 1,m
            sum = sum + wsave(ii+m*(j-1)+k) * wsave(iii+m*idivm+k)
          enddo ! k
c
          v(i) = sum
c
        enddo ! j
c
c
        return
        end
