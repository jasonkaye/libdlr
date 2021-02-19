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
        real*8 r1,diffmax,sum
        complex*16 v(len),v2(len),disp(len),disp2(len),
     1             wsave(len),w(len),v3(len)
c
        r1 = 1
c
c
        call prini(6,13)
c
c
        print *,'Enter n (n must be a positive integer power of 2):'
        read *,n
        call prinf('n = *',n,1)
c
c
        l = n/4
        call prinf('l = *',l,1)
c
c
c       Choose at random locations to calculate.
c
        call id_randperm(n/2,ind)
        ind(1) = 1
        ind(l) = n/2
        call prinf('ind = *',ind,l)
c
c
c       Fill the vector v to be transformed with N(0,1) variates.
c
        call normrand(n,v)
c
c
c       Copy v into v2.
c
        do k = 1,n/2
          v2(k) = v(k)
        enddo ! k
c
c
c       Transform v2 via dfftf2.
c
        call idd_ldiv(l,n,nblock)
        call dfftf2(nblock,n,v2,w)
c
        call prin2('v2 = *',v2,n)
c
c
        if(l .eq. 1) then
          disp2(1) = v2(ind(1))
        endif ! l .eq. 1
c
c
        if(l .gt. 1) then
c
c
c         Copy v into v3.
c
          do k = 1,n/2
            v3(k) = v(k)
          enddo ! k
c
c
c         Transform v3 via splitfft2.
c
          call splitfft2(nblock,n,v3,wsave,w)
c
          call prin2('w = *',w,n)
c
          do k = 1,l
            disp2(k) = w(ind(k))
          enddo ! k
c
c
c         Compare v2 and w (assuming that none of the complex-valued
c         entries in the same vector have the same absolute value).
c
          call matchconj(n/2,v2,w,sum)
          call prin2('v2 = *',v2,n)
          call prin2('w = *',w,n)
          call prin2('sum = *',sum,1)
c
c
        endif ! l .gt. 1
c
c
c       Initialize idd_sfft.
c
        call idd_sffti(l,ind,n,wsave)
c
c
c       Apply idd_sfft.
c
        call idd_sfft(l,ind,n,wsave,v)
        call prin2('v = *',v,n)
c
        do k = 1,l
          disp(k) = v(ind(k))
        enddo ! k
c
        call prin2('disp = *',disp,2*l)
        call prin2('disp2 = *',disp2,2*l)
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
        subroutine normrand(n,x)
c
c       constructs i.i.d. N(0,1) (pseudo)random variates via a simple
c       (but highly inefficient) scheme -- a stripped-down version
c       of the Box-Mueller-Marsaglia method.
c
c       input:
c       n -- number of i.i.d. draws to take
c
c       output:
c       x -- vector of i.i.d. N(0,1) random variates
c
        implicit none
        integer n,k
        real*8 x(n),a,b,twopi,r1
c
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
        do k = 1,n
c
          call id_srand(1,a)
          call id_srand(1,b)
c
          x(k) = sqrt(-2*log(a))*cos(twopi*b)
c
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine matchconj(n,u,v,sum)
c
c       sorts both u and v, and then checks whether any of their
c       entries differs by more than a conjugation. This routine
c       returns the sum of the absolute values of the differences
c       between the corresponding entries in the sorted arrays,
c       accumulating the smaller of the difference between each pair
c       of entries, or between one entry and the conjugate
c       of the other in the pair.
c
c       input:
c       n -- length of u and v
c       u -- array to be compared to v
c       v -- array to be compared to u
c
c       output:
c       sum -- the l^1 norm of the difference between the sorted
c              arrays u and v, adjusting for possible conjugations
c              to make sum as small as possible
c
c       _N.B._: This routine destroys u and v.
c
        implicit none
        integer n,k
        real*8 diff1,diff2,diff,sum
        complex*16 u(n),v(n)
c
c
c       Sort both u and v via the bubble process.
c
        call zbubble(n,u)
        call zbubble(n,v)
c
c
c       Calculate the l^1 norm of the difference between u and v,
c       after appropriate conjugations.
c
        sum = 0
c
        do k = 1,n
c
          diff1 = abs(u(k)-v(k))
          diff2 = abs(u(k)-conjg(v(k)))
          diff = min(diff1,diff2)
          sum = sum+diff
c
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine zbubble(n,v)
c
c       sort v via the bubble process, in order of increasing
c       absolute value.
c
c       input:
c       n -- length of v
c
c       output:
c       v -- sorted vector
c
        implicit none
        integer n,j,k
        complex*16 v(n),temp
c
c
        do j = 1,n-1
          do k = 1,n-j
c
            if(abs(v(k+1)) .lt. abs(v(k))) then
              temp = v(k)
              v(k) = v(k+1)
              v(k+1) = temp
            endif
c
          enddo ! k
        enddo ! j
c
c
        return
        end
c
c
c
c
        subroutine dfftf2(l,n,v,w)
c
c       calculate via a full FFT what splitfft2 is supposed to compute,
c       up to conjugation and permutation.
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
c       w -- must be at least 2*n+15 real*8 elements long
c
        implicit none
        integer n,l,k
        real*8 r1,fact,v(n),w(2*n+15)
c
        r1 = 1
c
c
c       Transpose v.
c
        call vectrans(l,n,v,w)
c
c
c       Transform v via dfftf and normalize it.
c
        call dffti(n,w)
        call dfftf(n,v,w)
c
        fact = 1/sqrt(r1*n)
        do k = 1,n
          v(k) = v(k)*fact
        enddo ! k
c
c
c       Shift most of v to the left by one element.
c
        call shifter(n,v)
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
c       v -- vector to be transposed
c
c       output:
c       v -- transposed vector
c
c       work:
c       w -- must be at least n real*8 elements long
c
        implicit none
        integer n,m,l,k,j
        real*8 v(n),w(n)
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
        subroutine splitfft2(l,n,v,wsave,w)
c
c       computes the DFT of v via a two-stage procedure,
c       composed with permutation matrices both on input and on output.
c       The output is actually half of the entries of the complex
c       output of a DFT of the real input. Which half is unclear,
c       but half of the entries of such a DFT are the conjugates
c       of the other half, and our chosen half contains only one
c       number from each conjugate pair. Thus, the output is some
c       permutation of the output of dfftf, followed by conjugations
c       of certain of the elements. The whole transform is therefore
c       orthogonal when viewed as a real-valued transform.
c
c       input:
c       l -- first stride; l must be at least 2
c       n -- length of v
c       v -- vector to be transformed
c
c       output:
c       w -- transformed vector
c
c       work:
c       wsave -- must be at least 2*(2*n+15) complex*16 elements long
c
c       _N.B._: l must be at least 2.
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
        integer n,m,l,k,j
        real*8 r1,fact,twopi,v(n)
        complex*16 w(n),wsave(2*(2*n+15)),ci
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
        m = n/l
c
c
c       FFT each block of length l of v.
c
        call dffti(l,wsave)
c
        do k = 1,m
          call dfftf(l,v(l*(k-1)+1),wsave)
        enddo ! k
c
c
c       Transpose v to obtain w.
c
        do k = 1,m
          do j = 1,l/2-1
            w(m*(j-1)+k) = v(l*(k-1)+2*j) + ci*v(l*(k-1)+2*j+1)
          enddo ! j
        enddo ! k
c
c       Handle the purely real frequency components separately.
c
        do k = 1,m
          w(m*(l/2-1)+k) = v(l*(k-1)+l)
        enddo ! k
c
c
c       Multiply by twiddle factors.
c
        do j = 1,l/2
          do k = 1,m
            w(m*(j-1)+k) = w(m*(j-1)+k)
     1                   * exp(-twopi*ci*(k-1)*(j)/(r1*n))
          enddo ! k
        enddo ! j
c
c
c       FFT each block of length m of w.
c
        call zffti(m,wsave(2*n+15+1))
c
        do j = 1,l/2
          call zfftf(m,w(m*(j-1)+1),wsave(2*n+15+1))
        enddo ! j
c
c       Handle the purely real frequency components separately.
c
        do k = 1,m/2
          w(m*(l/2-1)+m/2+k) = v(l*(2*k-2)+1) + ci*v(l*(2*k-1)+1)
        enddo ! k
c
        call dffti(m,wsave(2*n+15+1))
        call dfftf(m,w(m*(l/2-1)+m/2+1),wsave(2*n+15+1))
c
        call shifter(m,w(m*(l/2-1)+m/2+1))
c
c
c       Normalize w.
c
        fact = 1/sqrt(r1*n)
        do k = 1,n
          w(k) = w(k)*fact
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine shifter(m,v)
c
c       shifts all but the first and last element of v to the left
c       by one element, leaving the last element in place,
c       and placing the first element in the penultimate slot.
c
c       input:
c       m -- length of v
c       v -- vector to be shifted
c
c       output:
c       v -- shifted vector
c
        implicit none
        integer m,k
        real*8 v(m),temp
c
c
        temp = v(1)
c
        do k = 1,m-2
          v(k) = v(k+1)
        enddo ! k
c
        v(m-1) = temp
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
c       routine idd_sffti initializes routine idd_sfft.
c
c       routine idd_sfft rapidly computes a subset of the entries
c       of the DFT of a vector, composed with permutation matrices
c       both on input and on output.
c
c       routine idd_ldiv finds the greatest integer less than or equal
c       to a specified integer, that is divisible by another (larger)
c       specified integer.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idd_ldiv(l,n,m)
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
        subroutine idd_sffti(l,ind,n,wsave)
c
c       initializes wsave for using routine idd_sfft.
c
c       input:
c       l -- number of pairs of entries in the output of idd_sfft
c            to compute
c       ind -- indices of the pairs of entries in the output
c              of idd_sfft to compute; the indices must be chosen
c              in the range from 1 to n/2
c       n -- length of the vector to be transformed
c
c       output:
c       wsave -- array needed by routine idd_sfft for processing
c                (the present routine does not use the last n elements
c                 of wsave, but routine idd_sfft does)
c
        implicit none
        integer l,ind(l),n
        complex*16 wsave(2*l+15+4*n)
c
c
        if(l .eq. 1) call idd_sffti1(ind,n,wsave)
        if(l .gt. 1) call idd_sffti2(l,ind,n,wsave)
c
c
        return
        end
c
c
c
c
        subroutine idd_sffti1(ind,n,wsave)
c
c       routine idd_sffti serves as a wrapper around
c       the present routine; please see routine idd_sffti
c       for documentation.
c
        implicit none
        integer ind,n,k
        real*8 r1,twopi,wsave(2*(2+15+4*n)),fact
c
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
        fact = 1/sqrt(r1*n)
c
c
        do k = 1,n
          wsave(k) = cos(twopi*(k-1)*ind/(r1*n))*fact
        enddo ! k
c
        do k = 1,n
          wsave(n+k) = -sin(twopi*(k-1)*ind/(r1*n))*fact
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_sffti2(l,ind,n,wsave)
c
c       routine idd_sffti serves as a wrapper around
c       the present routine; please see routine idd_sffti
c       for documentation.
c
        implicit none
        integer l,ind(l),n,nblock,ii,m,idivm,imodm,i,j,k
        real*8 r1,twopi,fact
        complex*16 wsave(2*l+15+4*n),ci,twopii
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
        twopii = twopi*ci
c
c
c       Determine the block lengths for the FFTs.
c
        call idd_ldiv(l,n,nblock)
        m = n/nblock
c
c
c       Initialize wsave for using routine dfftf.
c
        call dffti(nblock,wsave)
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
c
          i = ind(j)
c
c
          if(i .le. n/2-m/2) then
c
            idivm = (i-1)/m
            imodm = (i-1)-m*idivm
c
            do k = 1,m
              wsave(ii+m*(j-1)+k) = exp(-twopii*(k-1)*imodm/(r1*m))
     1         * exp(-twopii*(k-1)*(idivm+1)/(r1*n)) * fact
            enddo ! k
c
          endif ! i .le. n/2-m/2
c
c
          if(i .gt. n/2-m/2) then
c
            idivm = i/(m/2)
            imodm = i-(m/2)*idivm
c
            do k = 1,m
              wsave(ii+m*(j-1)+k) = exp(-twopii*(k-1)*imodm/(r1*m))
     1                            * fact
            enddo ! k
c
          endif ! i .gt. n/2-m/2
c
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
        subroutine idd_sfft(l,ind,n,wsave,v)
c
c       computes a subset of the entries of the DFT of v,
c       composed with permutation matrices both on input and on output,
c       via a two-stage procedure (debugging code routine dfftf2 above
c       is supposed to calculate the full vector from which idd_sfft
c       returns a subset of the entries, when dfftf2 has
c       the same parameter nblock as in the present routine).
c
c       input:
c       l -- number of pairs of entries in the output to compute
c       ind -- indices of the pairs of entries in the output
c              to compute; the indices must be chosen
c              in the range from 1 to n/2
c       n -- length of v; n must be a positive integer power of 2
c       v -- vector to be transformed
c       wsave -- processing array initialized by routine idd_sffti
c
c       output:
c       v -- pairs of entries indexed by ind are given
c            their appropriately transformed values
c
c       _N.B._: n must be a positive integer power of 2.
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
        integer l,ind(l),n
        real*8 v(n)
        complex*16 wsave(2*l+15+4*n)
c
c
        if(l .eq. 1) call idd_sfft1(ind,n,v,wsave)
        if(l .gt. 1) call idd_sfft2(l,ind,n,v,wsave)
c
c
        return
        end
c
c
c
c
        subroutine idd_sfft1(ind,n,v,wsave)
c
c       routine idd_sfft serves as a wrapper around
c       the present routine; please see routine idd_sfft
c       for documentation.
c
        implicit none
        integer ind,n,k
        real*8 v(n),r1,twopi,sumr,sumi,fact,wsave(2*(2+15+4*n))
c
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
        if(ind .lt. n/2) then
c
c
          sumr = 0
c
          do k = 1,n
            sumr = sumr+wsave(k)*v(k)
          enddo ! k
c
c
          sumi = 0
c
          do k = 1,n
            sumi = sumi+wsave(n+k)*v(k)
          enddo ! k
c
c
        endif ! ind .lt. n/2
c
c
        if(ind .eq. n/2) then
c
c
          fact = 1/sqrt(r1*n)
c
c
          sumr = 0
c
          do k = 1,n
            sumr = sumr+v(k)
          enddo ! k
c
          sumr = sumr*fact
c
c
          sumi = 0
c
          do k = 1,n/2
            sumi = sumi+v(2*k-1)
            sumi = sumi-v(2*k)
          enddo ! k
c
          sumi = sumi*fact
c
c
        endif ! ind .eq. n/2
c
c
        v(2*ind-1) = sumr
        v(2*ind) = sumi
c
c
        return
        end
c
c
c
c
        subroutine idd_sfft2(l,ind,n,v,wsave)
c
c       routine idd_sfft serves as a wrapper around
c       the present routine; please see routine idd_sfft
c       for documentation.
c
        implicit none
        integer n,m,l,k,j,ind(l),i,idivm,nblock,ii,iii,imodm
        real*8 r1,twopi,v(n),rsum,fact
        complex*16 wsave(2*l+15+4*n),ci,sum
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
c       Determine the block lengths for the FFTs.
c
        call idd_ldiv(l,n,nblock)
c
c
        m = n/nblock
c
c
c       FFT each block of length nblock of v.
c
        do k = 1,m
          call dfftf(nblock,v(nblock*(k-1)+1),wsave)
        enddo ! k
c
c
c       Transpose v to obtain wsave(2*l+15+2*n+1 : 2*l+15+3*n).
c
        iii = 2*l+15+2*n
c
        do k = 1,m
          do j = 1,nblock/2-1
            wsave(iii+m*(j-1)+k) = v(nblock*(k-1)+2*j)
     1                           + ci*v(nblock*(k-1)+2*j+1)
          enddo ! j
        enddo ! k
c
c       Handle the purely real frequency components separately.
c
        do k = 1,m
          wsave(iii+m*(nblock/2-1)+k) = v(nblock*(k-1)+nblock)
          wsave(iii+m*(nblock/2)+k) = v(nblock*(k-1)+1)
        enddo ! k
c
c
c       Directly calculate the desired entries of v.
c
        ii = 2*l+15
c
        do j = 1,l
c
c
          i = ind(j)
c
c
          if(i .le. n/2-m/2) then
c
            idivm = (i-1)/m
            imodm = (i-1)-m*idivm
c
            sum = 0
c
            do k = 1,m
              sum = sum + wsave(iii+m*idivm+k) * wsave(ii+m*(j-1)+k)
            enddo ! k
c
            v(2*i-1) = sum
            v(2*i) = -ci*sum
c
          endif ! i .le. n/2-m/2
c
c
          if(i .gt. n/2-m/2) then
c
            if(i .lt. n/2) then
c
              idivm = i/(m/2)
              imodm = i-(m/2)*idivm
c
              sum = 0
c
              do k = 1,m
                sum = sum + wsave(iii+m*(nblock/2)+k)
     1              * wsave(ii+m*(j-1)+k)
              enddo ! k
c
              v(2*i-1) = sum
              v(2*i) = -ci*sum
c
            endif
c
            if(i .eq. n/2) then
c
              fact = 1/sqrt(r1*n)
c
c
              rsum = 0
c
              do k = 1,m
                rsum = rsum + wsave(iii+m*(nblock/2)+k)
              enddo ! k
c
              v(n-1) = rsum*fact
c
c
              rsum = 0
c
              do k = 1,m/2
                rsum = rsum + wsave(iii+m*(nblock/2)+2*k-1)
                rsum = rsum - wsave(iii+m*(nblock/2)+2*k)
              enddo ! k
c
              v(n) = rsum*fact
c
            endif
c
          endif ! i .gt. n/2-m/2
c
c
        enddo ! j
c
c
        return
        end
