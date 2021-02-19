c
c
c       dependencies: prini, dfft, id_rand, id_rtrans, idd_house,
c                     idd_qrpiv, idd_id, idd_sfft, lapack.a, blas.a,
c                     and (for the debugging code) idd_svd
c
c
        implicit none
c
        integer len
        parameter(len = 10 000 000)
c
        integer k,l,m,lm,ier,i,n,j,itype
        real*8 r(len),u(len),ru(len),v(len),s(len),cond,w(len)
c
c
        call prini(6,13)
c
c
        print *,'Enter k:'
        read *,k
        call prinf('k = *',k,1)
c
        l = k+4
        call prinf('l = *',l,1)
c
        print *,'Enter m:'
        read *,m
        call prinf('m = *',m,1)
c
c
c       Fill u with a matrix whose columns are orthonormal.
c
        itype = 1
        call fillortho(itype,m,k,u)
c
c
c       SVD u and display its greatest and least singular values
c       (both of which should be 1).
c
        do i = 1,m*k
          ru(i) = u(i)
        enddo ! i
c
        call iddr_svd(m,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        call prin2('s(1) = *',s(1),1)
        call prin2('s(k) = *',s(k),1)
c
c
c       Fill r with i.i.d. N(0,1) random variates.
c
        lm = l*m
        call normrand(lm,r)
c
c
c       Multiply r and u to obtain ru.
c
        call matmult(l,m,r,k,u,ru)
c
c
c       SVD ru and print its condition number.
c
        call iddr_svd(l,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        cond = s(1)/s(k)
        call prin2('cond = *',cond,1)
c
c
c       Initialize the random transforms.
c
        call idd_frmi(m,n,w)
c
c
c       Apply the random transforms to every column of u.
c
        do i = 1,k
          call idd_frm(m,n,w,u(1+m*(i-1)),v(1+n*(i-1)))
        enddo ! i
c
c
c       Copy the uppermost block of v into ru.
c
        do i = 1,k
          do j = 1,l
            ru(j+l*(i-1)) = v(j+n*(i-1))
          enddo ! j
        enddo ! i
c
c
c       SVD ru and print its condition number.
c
        call iddr_svd(l,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        cond = s(1)/s(k)
        call prin2('cond = *',cond,1)
c
c
c       Apply the subsampled random transforms to every column of u.
c
        call idd_sfrmi(l,m,n,w)
c
        do i = 1,k
          call idd_sfrm(l,m,n,w,u(1+m*(i-1)),ru(1+l*(i-1)))
        enddo ! i
c
c
c       SVD ru and print its condition number.
c
        call iddr_svd(l,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        cond = s(1)/s(k)
        call prin2('cond = *',cond,1)
c
c
        stop
        end
c
c
c
c
        subroutine fillortho(itype,m,k,u)
c
c       fills u with a matrix whose columns are orthonormal.
c
c       input:
c       itype -- specifies the matrix with which to fill u;
c                set to 1 for an identity matrix in the uppermost block
c                of u, and zeros elsewhere;
c                set to 2 for a subset of the columns of a random
c                unitary matrix;
c                set to 3 for a subset of the columns of a DCT-IV
c       m -- first dimension of u
c       k -- second dimension of u
c
c       output:
c       u -- matrix of the specified type whose columns are orthonormal
c
        implicit none
        integer m,k,i,j,itype,mk
        real*8 u(m,k),r1,pi
c
        r1 = 1
        pi = 4*atan(r1)
c
c
        if(itype .eq. 1) then
c
c         Put an identity matrix in the upper block of u,
c         and zeros elsewhere.
c
          do i = 1,k
            do j = 1,m
c
              u(j,i) = 0
c
            enddo ! j
          enddo ! i
c
          do i = 1,k
            u(i,i) = 1
          enddo ! i
c
        endif ! itype .eq. 1
c
c
        if(itype .eq. 2) then
c
c         Fill u with i.i.d. random variates drawn uniformly
c         from [0,1], and then orthonormalize the columns of u.
c
          mk = m*k
          call id_srand(mk,u)
          call orthonorm(m,k,u)
c
        endif ! itype .eq. 2
c
c
        if(itype .eq. 3) then
c
c         Fill u with part of the DCT-IV matrix.
c
          do i = 1,k
            do j = 1,m
              u(j,i) = cos(pi*(j-r1/2)*(i-r1/2)/m)*sqrt(r1*2/m)
            enddo ! j
          enddo ! i
c
        endif ! itype .eq. 3
c
c
        return
        end
c
c
c
c
        subroutine matmult(l,m,a,n,b,c)
c
c       multiplies a and b to obtain c.
c
c       input:
c       l -- first dimension of a and c
c       m -- second dimension of a, and first dimension of b
c       a -- leftmost matrix in the product c = a b
c       n -- second dimension of b and c
c       b -- rightmost matrix in the product c = a b
c
c       output:
c       c -- the product c = a b
c
        implicit none
        integer l,m,n,i,j,k
        real*8 a(l,m),b(m,n),c(l,n)
c
c
        do i = 1,l
          do k = 1,n
c
            c(i,k) = 0
            do j = 1,m
              c(i,k) = c(i,k) + a(i,j)*b(j,k)
            enddo ! j
c
          enddo ! k
        enddo ! i
c
c
        return
        end
c
c
c
c
        subroutine orthonorm(m,n,a)
c
c       orthonormalizes the columns of a
c       via the Gram-Schmidt process,
c       assuming that a has full rank
c       and does not require pivoting.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a
c       a -- matrix to orthonormalize
c
c       output:
c       a -- orthonormalized matrix
c
        implicit none
        integer m,n,j,k,l,loop
        real*8 a(m,n),rms,prod
c
c
        if(m .lt. n) then
          call prinf('bombing from orthonorm, since m < n....*',m,0)
          call prinf('m = *',m,1)
          call prinf('n = *',n,1)
          stop
        endif
c
c
        do k = 1,n
c
c
c         Calculate the root-mean-square of the entries
c         of the entries of column k.
c
          rms = 0
          do j = 1,m
            rms = rms + (a(j,k))**2
          enddo ! j
          rms = sqrt(rms)
c
c         Normalize column k.
c
          do j = 1,m
            a(j,k) = a(j,k)/rms
          enddo ! j
c
c
          if(k .lt. n) then
            do loop = 1,2
              do l = k+1,n
c
c               Compute the inner product of column k and column l.
c
                prod = 0
                do j = 1,m
                  prod = prod + a(j,k)*a(j,l)
                enddo ! j
c
c               Subtract off the component for column k in column l.
c
                do j = 1,m
                  a(j,l) = a(j,l) - prod*a(j,k)
                enddo ! j
c
              enddo ! l
            enddo ! loop
          endif ! k .lt. n
c
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
        subroutine normrand(n,x)
c
c       constructs i.i.d. N(0,1) (pseudo)random variates via a simple
c       (but highly inefficient) scheme -- a stripped-down version
c       of the Box-Mueller-Marsaglia method.
c
c       input:
c       n -- number of i.i.d. N(0,1) random variates to generate
c
c       output:
c       x -- vector whose entries are i.i.d. N(0,1) random variates
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine idd_frm transforms a vector via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c
c       routine idd_sfrm transforms a vector into a vector
c       of specified length via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c
c       routine idd_frmi initializes routine idd_frm.
c
c       routine idd_sfrmi initializes routine idd_sfrm.
c
c       routine idd_pairsamps calculates the indices of the pairs
c       of integers to which the individual integers
c       in a specified set belong.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idd_frm(m,n,w,x,y)
c
c       transforms x into y via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c       In contrast to routine idd_sfrm, the present routine works best
c       when the length of the transformed vector is the integer n
c       output by routine idd_frmi, or when the length
c       is not specified, but instead determined a posteriori
c       using the output of the present routine. The transformed vector
c       output by the present routine is randomly permuted.
c
c       input:
c       m -- length of x
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m, as obtained
c            from the routine idd_frmi; n is the length of y
c       w -- initialization array constructed by routine idd_frmi
c       x -- vector to be transformed
c
c       output:
c       y -- transform of x
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,iw,n,k
        real*8 w(17*m+70),x(m),y(n)
c
c
c       Apply Rokhlin's random transformation to x, obtaining
c       w(16*m+71 : 17*m+70).
c
        iw = w(3+m+n)
        call idd_random_transf(x,w(16*m+70+1),w(iw))
c
c
c       Subselect from  w(16*m+71 : 17*m+70)  to obtain y.
c
        call idd_subselect(n,w(3),m,w(16*m+70+1),y)
c
c
c       Copy y into  w(16*m+71 : 16*m+n+70).
c
        do k = 1,n
          w(16*m+70+k) = y(k)
        enddo ! k
c
c
c       Fourier transform  w(16*m+71 : 16*m+n+70).
c
        call dfftf(n,w(16*m+70+1),w(4+m+n))
c
c
c       Permute  w(16*m+71 : 16*m+n+70)  to obtain y.
c
        call idd_permute(n,w(3+m),w(16*m+70+1),y)
c
c
        return
        end
c
c
c
c
        subroutine idd_sfrm(l,m,n,w,x,y)
c
c       transforms x into y via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c       In contrast to routine idd_frm, the present routine works best
c       when the length l of the transformed vector is known a priori.
c
c       input:
c       l -- length of y; l must be less than or equal to n
c       m -- length of x
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m, as obtained
c            from the routine idd_sfrmi
c       w -- initialization array constructed by routine idd_sfrmi
c       x -- vector to be transformed
c
c       output:
c       y -- transform of x
c
c       _N.B._: l must be less than or equal to n.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,iw,n,l,l2
        real*8 w(27*m+90),x(m),y(l)
c
c
c       Retrieve the number of pairs of outputs to be calculated
c       via sfft.
c
        l2 = w(3)
c
c
c       Apply Rokhlin's random transformation to x, obtaining
c       w(25*m+91 : 26*m+90).
c
        iw = w(4+m+l+l2)
        call idd_random_transf(x,w(25*m+90+1),w(iw))
c
c
c       Subselect from  w(25*m+91 : 26*m+90)  to obtain
c       w(26*m+91 : 26*m+n+90).
c
        call idd_subselect(n,w(4),m,w(25*m+90+1),w(26*m+90+1))
c
c
c       Fourier transform  w(26*m+91 : 26*m+n+90).
c
        call idd_sfft(l2,w(4+m+l),n,w(5+m+l+l2),w(26*m+90+1))
c
c
c       Copy the desired entries from  w(26*m+91 : 26*m+n+90)
c       to y.
c
        call idd_subselect(l,w(4+m),n,w(26*m+90+1),y)
c
c
        return
        end
c
c
c
c
        subroutine idd_pairsamps(n,l,ind,l2,ind2,marker)
c
c       calculates the indices of the l2 pairs of integers
c       to which the l individual integers from ind belong.
c       The integers in ind may range from 1 to n.
c
c       input:
c       n -- upper bound on the integers in ind
c            (the number 1 must be a lower bound);
c            n must be even
c       l -- length of ind
c       ind -- integers selected from 1 to n
c
c       output:
c       l2 -- length of ind2
c       ind2 -- indices in the range from 1 to n/2 of the pairs
c               of integers to which the entries of ind belong
c
c       work:
c       marker -- must be at least n/2 integer elements long
c
c       _N.B._: n must be even.
c
        implicit none
        integer l,n,ind(l),ind2(l),marker(n/2),l2,k
c
c
c       Unmark all pairs.
c
        do k = 1,n/2
          marker(k) = 0
        enddo ! k
c
c
c       Mark the required pairs.
c
        do k = 1,l
          marker((ind(k)+1)/2) = marker((ind(k)+1)/2)+1
        enddo ! k
c
c
c       Record the required pairs in indpair.
c
        l2 = 0
c
        do k = 1,n/2
c
          if(marker(k) .ne. 0) then
            l2 = l2+1
            ind2(l2) = k
          endif
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
        subroutine idd_permute(n,ind,x,y)
c
c       copy the entries of x into y, rearranged according
c       to the permutation specified by ind.
c
c       input:
c       n -- length of ind, x, and y
c       ind -- permutation of n objects
c       x -- vector to be permuted
c
c       output:
c       y -- permutation of x
c
        implicit none
        integer n,ind(n),k
        real*8 x(n),y(n)
c
c
        do k = 1,n
          y(k) = x(ind(k))
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_subselect(n,ind,m,x,y)
c
c       copies into y the entries of x indicated by ind.
c
c       input:
c       n -- number of entries of x to copy into y
c       ind -- indices of the entries in x to copy into y
c       m -- length of x
c       x -- vector whose entries are to be copied
c
c       output:
c       y -- collection of entries of x specified by ind
c
        implicit none
        integer n,ind(n),m,k
        real*8 x(m),y(n)
c
c
        do k = 1,n
          y(k) = x(ind(k))
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_frmi(m,n,w)
c
c       initializes data for the routine idd_frm.
c
c       input:
c       m -- length of the vector to be transformed
c
c       output:
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m
c       w -- initialization array to be used by routine idd_frm
c
c
c       glossary for the fully initialized w:
c
c       w(1) = m
c       w(2) = n
c       w(3:2+m) stores a permutation of m objects
c       w(3+m:2+m+n) stores a permutation of n objects
c       w(3+m+n) = address in w of the initialization array
c                  for idd_random_transf
c       w(4+m+n:int(w(3+m+n))-1) stores the initialization array
c                                for dfft
c       w(int(w(3+m+n)):16*m+70) stores the initialization array
c                                for idd_random_transf
c
c
c       _N.B._: n is an output of the present routine;
c               this routine changes n.
c
c
        implicit none
        integer m,n,l,nsteps,keep,lw,ia
        real*8 w(17*m+70)
c
c
c       Find the greatest integer less than or equal to m
c       which is a power of two.
c
        call idd_poweroftwo(m,l,n)
c
c
c       Store m and n in w.
c
        w(1) = m
        w(2) = n
c
c
c       Store random permutations of m and n objects in w.
c
        call id_randperm(m,w(3))
        call id_randperm(n,w(3+m))
c
c
c       Store the address within w of the idd_random_transf_init
c       initialization data.
c
        ia = 4+m+n+2*n+15
        w(3+m+n) = ia
c
c
c       Store the initialization data for dfft in w.
c
        call dffti(n,w(4+m+n))
c
c
c       Store the initialization data for idd_random_transf_init in w.
c
        nsteps = 3
        call idd_random_transf_init(nsteps,m,w(ia),keep)
c
c
c       Calculate the total number of elements used in w.
c
        lw = 3+m+n+2*n+15 + 3*nsteps*m+2*m+m/4+50
c
        if(16*m+70 .lt. lw) then
          call prinf('lw = *',lw,1)
          call prinf('16m+70 = *',16*m+70,1)
          stop
        endif
c
c
        return
        end
c
c
c
c
        subroutine idd_sfrmi(l,m,n,w)
c
c       initializes data for the routine idd_sfrm.
c
c       input:
c       l -- length of the transformed (output) vector
c       m -- length of the vector to be transformed
c
c       output:
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m
c       w -- initialization array to be used by routine idd_sfrm
c
c
c       glossary for the fully initialized w:
c
c       w(1) = m
c       w(2) = n
c       w(3) = l2
c       w(4:3+m) stores a permutation of m objects
c       w(4+m:3+m+l) stores the indices of the l outputs which idd_sfft
c                    calculates
c       w(4+m+l:3+m+l+l2) stores the indices of the l2 pairs of outputs
c                         which idd_sfft calculates
c       w(4+m+l+l2) = address in w of the initialization array
c                     for idd_random_transf
c       w(5+m+l+l2:int(w(4+m+l+l2))-1) stores the initialization array
c                                      for idd_sfft
c       w(int(w(4+m+l+l2)):25*m+90) stores the initialization array
c                                   for idd_random_transf
c
c
c       _N.B._: n is an output of the present routine;
c               this routine changes n.
c
c
        implicit none
        integer l,m,n,idummy,nsteps,keep,lw,l2,ia
        real*8 w(27*m+90)
c
c
c       Find the greatest integer less than or equal to m
c       which is a power of two.
c
        call idd_poweroftwo(m,idummy,n)
c
c
c       Store m and n in w.
c
        w(1) = m
        w(2) = n
c
c
c       Store random permutations of m and n objects in w.
c
        call id_randperm(m,w(4))
        call id_randperm(n,w(4+m))
c
c
c       Find the pairs of integers covering the integers in
c       w(4+m : 3+m+(l+1)/2).
c
        call idd_pairsamps(n,l,w(4+m),l2,w(4+m+2*l),w(4+m+3*l))
        w(3) = l2
        call idd_copyints(l2,w(4+m+2*l),w(4+m+l))
c
c
c       Store the address within w of the idd_random_transf_init
c       initialization data.
c
        ia = 5+m+l+l2+4*l2+30+8*n
        w(4+m+l+l2) = ia
c
c
c       Store the initialization data for idd_sfft in w.
c
        call idd_sffti(l2,w(4+m+l),n,w(5+m+l+l2))
c
c
c       Store the initialization data for idd_random_transf_init in w.
c
        nsteps = 3
        call idd_random_transf_init(nsteps,m,w(ia),keep)
c
c
c       Calculate the total number of elements used in w.
c
        lw = 4+m+l+l2+4*l2+30+8*n + 3*nsteps*m+2*m+m/4+50
c
        if(25*m+90 .lt. lw) then
          call prinf('lw = *',lw,1)
          call prinf('25m+90 = *',25*m+90,1)
          stop
        endif
c
c
        return
        end
c
c
c
c
        subroutine idd_copyints(n,ia,ib)
c
c       copies ia into ib.
c
c       input:
c       n -- length of ia and ib
c       ia -- array to be copied
c
c       output:
c       ib -- copy of ia
c
        implicit none
        integer n,ia(n),ib(n),k
c
c
        do k = 1,n
          ib(k) = ia(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_poweroftwo(m,l,n)
c
c       computes l = floor(log_2(m)) and n = 2**l.
c
c       input:
c       m -- integer whose log_2 is to be taken
c
c       output:
c       l -- floor(log_2(m))
c       n -- 2**l
c
        implicit none
        integer l,m,n
c
c
        l = 0
        n = 1
c
 1000   continue
          l = l+1
          n = n*2
        if(n .le. m) goto 1000
c
        l = l-1
        n = n/2
c
c
        return
        end
