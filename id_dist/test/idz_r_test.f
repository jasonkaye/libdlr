c
c
c       dependencies: prini, idz_house, idz_qrpiv, idz_id, id_rand,
c                     idzr_rid, idzp_rid, idz_id2svd, idzr_rsvd,
c                     idzp_rsvd, idz_snorm, lapack.a, blas.a
c
c
        implicit none
c
        integer len
        parameter(len = 10 000 000)
c
        integer m,n,krank,krank2,krank3,its,ier,iu3,iv3,is3,lw,k
        real*8 s(len),s2(len),s3(len),diff2,diff3,eps,r1
        complex*16 u(len),v(len),w(len),u2(len),v2(len),
     1             u3(len),v3(len)
        external matvec,matveca
c
        r1 = 1
c
c
        call prini(6,13)
c
c
        print *,'Enter m:'
        read *,m
        call prinf('m = *',m,1)
c
        print *,'Enter n:'
        read *,n
        call prinf('n = *',n,1)
c
        krank = 12
        call prinf('krank = *',krank,1)
c
c
c       Fill the singular vectors and singular values.
c
        call fill(krank,m,n,u,v,s)
        call prin2('s = *',s,krank)
c
c
c
c       SVD the matrix  u diag(s) v^*.
c
        krank2 = krank-1
c
        call idzr_rsvd(m,n,matveca,krank,u,v,s,
     1                 matvec,krank,u,v,s,krank2,u2,v2,s2,ier,w)
c
        if(ier .ne. 0) then
          call prinf('ier = *',ier,1)
          stop 20
        endif
c
        call prin2('s2 = *',s2,krank2)
c
c
c       Estimate the spectral norm of
c       u diag(s) v^* - u2 diag(s2) (v2)^*.
c
        its = 100
c
        call idz_diffsnorm(m,n,matveca,krank,u,v,s,
     1                     matveca,krank2,u2,v2,s2,
     2                     matvec,krank,u,v,s,
     3                     matvec,krank2,u2,v2,s2,
     4                     its,diff2,w)
c
        call prin2('diff2 = *',diff2,1)
c
        if(diff2 .gt. sqrt(r1*m*n)*s(krank2+1)) then
          stop 200
        endif
c
c
c
c       SVD the matrix  u diag(s) v^*.
c
        eps = .5d-9/sqrt(sqrt(r1*m*n))
        lw = len
c
        call idzp_rsvd(lw,eps,m,n,matveca,krank,u,v,s,
     1                 matvec,krank,u,v,s,krank3,iu3,iv3,is3,w,ier)
c
        if(ier .ne. 0) then
          call prinf('ier = *',ier,1)
          stop 30
        endif
c
c
c       Copy u3, v3, and s3 from w.
c
        do k = 1,krank3*m
          u3(k) = w(iu3+k-1)
        enddo ! k
c
        do k = 1,krank3*n
          v3(k) = w(iv3+k-1)
        enddo ! k
c
        do k = 1,krank3
          s3(k) = w(is3+k-1)
        enddo ! k
        call prin2('s3 = *',s3,krank3)
c
c
c       Estimate the spectral norm of
c       u diag(s) v^* - u3 diag(s3) (v3)^*.
c
        its = 100
c
        call idz_diffsnorm(m,n,matveca,krank,u,v,s,
     1                     matveca,krank3,u3,v3,s3,
     2                     matvec,krank,u,v,s,
     3                     matvec,krank3,u3,v3,s3,
     4                     its,diff3,w)
c
        call prin2('diff3 = *',diff3,1)
c
        if(diff3 .gt. sqrt(r1*m*n)*eps) then
          stop 300
        endif
c
c
c
        stop
        end
c
c
c
c
        subroutine matveca(m,x,n,y,krank,u,v,s)
c
c       applies the adjoint of  u diag(s) v^*  to x, obtaining y.
c
c       input:
c       m -- length of x
c       x -- vector to which the adjoint of  u diag(s) v^*
c            is to be applied
c       n -- length of y
c       krank -- second dimension of u and v
c       u -- leftmost matrix in  u diag(s) v^*
c       v -- rightmost matrix in  u diag(s) v^*
c       s -- vector in  u diag(s) v^*
c
c       output:
c       y -- product of the adjoint of  u diag(s) v^*  and x
c
        implicit none
        integer m,n,krank,j,k
        real*8 s(krank)
        complex*16 x(m),y(n),u(m,krank),v(n,krank),
     1             sux(10 000 000),sum
c
c
c       Form sux = diag(s) u^* x.
c
        do k = 1,krank
c
          sum = 0
c
          do j = 1,m
            sum = sum+s(k)*conjg(u(j,k))*x(j)
          enddo ! j
c
          sux(k) = sum
c
        enddo ! k
c
c
c       Form y = v sux = v diag(s) u^* x.
c
        do k = 1,n
          y(k) = 0
        enddo ! k
c
        do j = 1,krank
          do k = 1,n
            y(k) = y(k)+v(k,j)*sux(j)
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
        subroutine matvec(n,y,m,x,krank,u,v,s)
c
c       applies u diag(s) v^* to y, obtaining x.
c
c       input:
c       n -- length of y
c       y -- vector to which  u diag(s) v^*  is to be applied
c       m -- length of x
c       krank -- second dimension of u and v
c       u -- leftmost matrix in  u diag(s) v^*
c       v -- rightmost matrix in  u diag(s) v^*
c       s -- vector in  u diag(s) v^*
c
c       output:
c       x -- product of  u diag(s) v^*  and y
c
        implicit none
        integer m,n,krank,j,k
        real*8 s(krank)
        complex*16 x(m),y(n),u(m,krank),v(n,krank),
     1             svy(10 000 000),sum
c
c
c       Form svy = diag(s) v^* y.
c
        do j = 1,krank
c
          sum = 0
c
          do k = 1,n
            sum = sum+s(j)*conjg(v(k,j))*y(k)
          enddo ! k
c
          svy(j) = sum
c
        enddo ! j
c
c
c       Form x = u svy = u diag(s) v^* y.
c
        do j = 1,m
          x(j) = 0
        enddo ! j
c
        do k = 1,krank
          do j = 1,m
            x(j) = x(j)+u(j,k)*svy(k)
          enddo ! j
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine fill(krank,m,n,u,v,s)
c
c       constructs singular vectors and singular values
c       for an m x n matrix.
c
c       input:
c       krank -- length of s, and second dimension of u and v
c       m -- first dimension of u
c       n -- first dimension of v
c
c       output:
c       u -- matrix of normalized left singular vectors
c       v -- matrix of normalized right singular vectors
c       s -- vector of singular values
c
        implicit none
        integer krank,m,n,m2krank,n2krank,j,k,l
        real*8 s(krank),r1
        complex*16 u(m,krank),v(n,krank)
c
        r1 = 1
c
c
c       Fill the real and imaginary parts of every entry of u and v
c       with i.i.d. random variables drawn uniformly from [-1,1].
c
        m2krank = m*2*krank
        call id_srand(m2krank,u)
c
        do l = 1,krank
          do j = 1,m
            u(j,l) = 2*u(j,l)-1
          enddo ! j
        enddo ! l
c
        n2krank = n*2*krank
        call id_srand(n2krank,v) 
c
        do l = 1,krank
          do k = 1,n
            v(k,l) = 2*v(k,l)-1
          enddo ! k
        enddo ! l
c
c
c       Orthonormalize the columns of u and v.
c
        call orthonorm(m,krank,u)
        call orthonorm(n,krank,v)
c
c
c       Fill s with exponentially decaying values.
c
        do k = 1,krank-1
          s(k) = exp(log(1d-10)*(k-1)/(krank-2))
        enddo ! k
c
        s(krank) = 1d-10
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
