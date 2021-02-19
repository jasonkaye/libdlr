c
c
c       dependencies: prini, dfft, id_rand, id_rtrans, idz_house,
c                     idz_qrpiv, idz_id, idz_sfft, lapack.a, blas.a,
c                     and (for the debugging code) idz_svd
c
c
        implicit none
c
        integer len
        parameter(len = 10 000 000)
c
        integer k,l,m,lm,ier,i,n,j,itype
        real*8 s(len),cond
        complex*16 r(len),u(len),ru(len),v(len),w(len)
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
        itype = 3
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
        call idzr_svd(m,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        call prin2('s(1) = *',s(1),1)
        call prin2('s(k) = *',s(k),1)
c
c
c       Fill the real and imaginary parts of the entries of r
c       with i.i.d. N(0,1) random variates.
c
        lm = 2*l*m
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
        call idzr_svd(l,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        cond = s(1)/s(k)
        call prin2('cond = *',cond,1)
c
c
c       Initialize the random transforms.
c
        call idz_frmi(m,n,w)
c
c
c       Apply the random transforms to every column of u.
c
        do i = 1,k
          call idz_frm(m,n,w,u(1+m*(i-1)),v(1+n*(i-1)))
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
        call idzr_svd(l,k,ru,k,w,v,s,ier,r)
        call prinf('ier = *',ier,1)
c
        cond = s(1)/s(k)
        call prin2('cond = *',cond,1)
c
c
c       Apply the subsampled random transforms to every column of u.
c
        call idz_sfrmi(l,m,n,w)
c
        do i = 1,k
          call idz_sfrm(l,m,n,w,u(1+m*(i-1)),ru(1+l*(i-1)))
        enddo ! i
c
c
c       SVD ru and print its condition number.
c
        call idzr_svd(l,k,ru,k,w,v,s,ier,r)
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
c                set to 3 for a subset of the columns of a DFT
c       m -- first dimension of u
c       k -- second dimension of u
c
c       output:
c       u -- matrix of the specified type whose columns are orthonormal
c
        implicit none
        integer m,k,i,j,itype,mk
        real*8 r1
        complex*16 u(m,k),twopii
c
        r1 = 1
        twopii = 2*4*atan(r1)*(0,1)
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
c         Fill u with i.i.d. random variates whose real and imaginary
c         parts are drawn uniformly from [0,1], and then
c         orthonormalize the columns of u.
c
          mk = 2*m*k
          call id_srand(mk,u)
          call orthonorm(m,k,u)
c
        endif ! itype .eq. 2
c
c
        if(itype .eq. 3) then
c
c         Fill u with part of the DFT matrix.
c
          do i = 1,k
            do j = 1,m
              u(j,i) = exp(-twopii*(j-r1)*(i-r1)/m)*sqrt(r1/m)
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
        complex*16 a(l,m),b(m,n),c(l,n)
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
        real*8 rms
        complex*16 a(m,n),prod
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
            rms = rms + (a(j,k))*conjg(a(j,k))
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
                  prod = prod + conjg(a(j,k))*a(j,l)
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
