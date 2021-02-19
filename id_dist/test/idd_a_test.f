c
c
c       dependencies: prini, idd_house, idd_qrpiv, idd_id, idd_svd,
c                     id_rand, idd_sfft, id_rtrans, idd_frm,
c                     iddr_aid, iddp_aid, idd_id2svd, iddr_asvd,
c                     iddp_asvd, idd_snorm, dfft, lapack.a, blas.a
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,krank2,krank3,ier,iu3,iv3,is3,lw,k,n2,
     1          krank4,krank5,iu5,iv5,is5,its
        real*8 s(len),s2(len),s3(len),eps,r1,b(len),a4(len),a5(len),
     1         errmax2,errrms2,errmax3,errrms3,
     2         errmax4,errrms4,errmax5,errrms5,
     3         u(len),v(len),w(len),u2(len),v2(len),a2(len),
     4         u3(len),v3(len),a(len),winitr(len),winitp(len),a3(len),
     5         u4(len),v4(len),s4(len),u5(len),v5(len),s5(len),
     6         diff2,diff3
        external matvec,matvect
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
        call fill(krank,m,n,a,u,v,s)
        call prin2('s = *',s,krank)
c
c
c
c       Initialize winitr for iddr_asvd.
c
        krank2 = krank-1
        call iddr_aidi(m,n,krank2,winitr)
c
c
c       Initialize winitp for iddp_asvd.
c
        call idd_frmi(m,n2,winitp)
c
c
c
c       SVD the matrix  a = u diag(s) v^T.
c
        call iddr_asvd(m,n,a,krank2,winitr,u2,v2,s2,ier)
c
        if(ier .ne. 0) then
          call prinf('ier = *',ier,1)
          stop 20
        endif
c
        call prin2('s2 = *',s2,krank2)
c
c
c       Efficiently estimate the spectral norm of
c       u diag(s) v^T - u2 diag(s2) (v2)^T.
c
        its = 10
c
        call idd_diffsnorm(m,n,matvect,krank,u,v,s,
     1                     matvect,krank2,u2,v2,s2,
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
c       Check the difference between a and u2 diag(s2) v2^T.
c
        call reconsvd(m,krank2,u2,s2,n,v2,a2)
        call materr(m,n,a,a2,errmax2,errrms2)
c
        call prin2('errmax2 = *',errmax2,1)
        call prin2('errrms2 = *',errrms2,1)
c
        if(errrms2 .gt. sqrt(r1*m*n)*s(krank2+1)) then
          stop 2000
        endif
c
c
c
c       SVD the matrix  a = u diag(s) v^T.
c
        eps = .5d-9
        lw = len
c
        call iddp_asvd(lw,eps,m,n,a,winitp,krank3,iu3,iv3,is3,w,ier)
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
c       Efficiently estimate the spectral norm of
c       u diag(s) v^T - u3 diag(s3) (v3)^T.
c
        its = 10
c
        call idd_diffsnorm(m,n,matvect,krank,u,v,s,
     1                     matvect,krank3,u3,v3,s3,
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
c       Check the difference between a and u3 diag(s3) v3^T.
c
        call reconsvd(m,krank3,u3,s3,n,v3,a3)
        call materr(m,n,a,a3,errmax3,errrms3)
c
        call prin2('errmax3 = *',errmax3,1)
        call prin2('errrms3 = *',errrms3,1)
c
        if(errrms3 .gt. sqrt(r1*m*n)*eps) then
          stop 3000
        endif
c
c
c
c       Copy a into b.
c
        do k = 1,m*n
          b(k) = a(k)
        enddo ! k
c
c
c       SVD b.
c
        krank4 = krank-1
c
        call iddr_svd(m,n,b,krank4,u4,v4,s4,ier,w)
c
        if(ier .ne. 0) then
          call prinf('ier = *',ier,1)
          stop 40
        endif
c
        call prin2('s4 = *',s4,krank4)
c
c
c       Check the difference between a and u4 diag(s4) v4^T.
c
        call reconsvd(m,krank4,u4,s4,n,v4,a4)
        call materr(m,n,a,a4,errmax4,errrms4)
c
        call prin2('errmax4 = *',errmax4,1)
        call prin2('errrms4 = *',errrms4,1)
c
        if(errrms4 .gt. sqrt(r1*m*n)*s(krank4+1)) then
          stop 400
        endif
c
c
c
c       Copy a into b.
c
        do k = 1,m*n
          b(k) = a(k)
        enddo ! k
c
c
c       SVD b.
c
        krank5 = krank-1
c
        call iddp_svd(lw,eps,m,n,b,krank5,iu5,iv5,is5,w,ier)
c
        if(ier .ne. 0) then
          call prinf('ier = *',ier,1)
          stop 50
        endif
c
c
c       Copy u5, v5, and s5 from w.
c
        do k = 1,krank5*m
          u5(k) = w(iu5+k-1)
        enddo ! k
c
        do k = 1,krank5*n
          v5(k) = w(iv5+k-1)
        enddo ! k
c
        do k = 1,krank5
          s5(k) = w(is5+k-1)
        enddo ! k
        call prin2('s5 = *',s5,krank5)
c
c
c       Check the difference between a and u5 diag(s5) v5^T.
c
        call reconsvd(m,krank5,u5,s5,n,v5,a5)
        call materr(m,n,a,a5,errmax5,errrms5)
c
        call prin2('errmax5 = *',errmax5,1)
        call prin2('errrms5 = *',errrms5,1)
c
        if(errrms5 .gt. sqrt(r1*m*n)*eps) then
          stop 500
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
        subroutine matvect(m,x,n,y,krank,u,v,s)
c
c       applies the transpose of  u diag(s) v^T  to x, obtaining y.
c
c       input:
c       m -- length of x
c       x -- vector to which the transpose of  u diag(s) v^T
c            is to be applied
c       n -- length of y
c       krank -- second dimension of u and v
c       u -- leftmost matrix in  u diag(s) v^T
c       v -- rightmost matrix in  u diag(s) v^T
c       s -- vector in  u diag(s) v^T
c
c       output:
c       y -- product of the transpose of  u diag(s) v^T  and x
c
        implicit none
        integer m,n,krank,j,k
        real*8 s(krank),x(m),y(n),u(m,krank),v(n,krank),
     1         sux(10 000 000),sum
c
c
c       Form sux = diag(s) u^T x.
c
        do k = 1,krank
c
          sum = 0
c
          do j = 1,m
            sum = sum+s(k)*u(j,k)*x(j)
          enddo ! j
c
          sux(k) = sum
c
        enddo ! k
c
c
c       Form y = v sux = v diag(s) u^T x.
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
c       applies u diag(s) v^T to y, obtaining x.
c
c       input:
c       n -- length of y
c       y -- vector to which  u diag(s) v^T  is to be applied
c       m -- length of x
c       krank -- second dimension of u and v
c       u -- leftmost matrix in  u diag(s) v^T
c       v -- rightmost matrix in  u diag(s) v^T
c       s -- vector in  u diag(s) v^T
c
c       output:
c       x -- product of  u diag(s) v^T  and y
c
        implicit none
        integer m,n,krank,j,k
        real*8 s(krank),x(m),y(n),u(m,krank),v(n,krank),
     1         svy(10 000 000),sum
c
c
c       Form svy = diag(s) v^T y.
c
        do j = 1,krank
c
          sum = 0
c
          do k = 1,n
            sum = sum+s(j)*v(k,j)*y(k)
          enddo ! k
c
          svy(j) = sum
c
        enddo ! j
c
c
c       Form x = u svy = u diag(s) v^T y.
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
        subroutine materr(m,n,a,b,errmax,errrms)
c
c       calculates the relative maximum and root-mean-square errors
c       corresponding to how much a and b differ.
c
c       input:
c       m -- first dimension of a and b
c       n -- second dimension of a and b
c       a -- matrix whose difference from b will be measured
c       b -- matrix whose difference from a will be measured
c
c       output:
c       errmax -- ratio of the maximum elementwise absolute difference
c                 between a and b to the maximum magnitude
c                 of all the elements of a
c       errrms -- ratio of the root-mean-square of the elements
c                 of the difference of a and b to the root-mean-square
c                 of all the elements of a
c
        implicit none
        integer m,n,j,k
        real*8 errmax,errrms,diff,amax,arss,a(m,n),b(m,n)
c
c
c       Calculate the maximum magnitude amax of the elements of a
c       and the root-sum-square arss of the elements of a.
c
        amax = 0
        arss = 0
c
        do k = 1,n
          do j = 1,m
c
            if(abs(a(j,k)) .gt. amax) amax = abs(a(j,k))
            arss = arss+a(j,k)**2
c
          enddo ! j
        enddo ! k
c
        arss = sqrt(arss)
c
c
c       Calculate the maximum elementwise absolute difference
c       between a and b, as well as the root-sum-square errrms
c       of the elements of the difference of a and b.
c
        errmax = 0
        errrms = 0
c
        do k = 1,n
          do j = 1,m
c
            diff = abs(a(j,k)-b(j,k))
c
            if(diff .gt. errmax) errmax = diff
            errrms = errrms+diff**2
c
          enddo ! j
        enddo ! k
c
        errrms = sqrt(errrms)
c
c
c       Calculate relative errors.
c
        errmax = errmax/amax
        errrms = errrms/arss
c
c
        return
        end
c
c
c
c
        subroutine reconsvd(m,krank,u,s,n,v,a)
c
c       forms a = u diag(s) v^T.
c
c       input:
c       m -- first dimension of u and a
c       krank -- size of s, and second dimension of u and v
c       u -- leftmost matrix in the product a = u diag(s) v^T
c       s -- entries on the diagonal in the middle matrix
c            in the product a = u diag(s) v^T
c       n -- second dimension of a, and first dimension of v
c       v -- rightmost matrix in the product a = u diag(s) v^T
c
c       output:
c       a -- matrix product u diag(s) v^T
c
        implicit none
        integer m,n,krank,j,k,l
        real*8 s(krank),u(m,krank),v(n,krank),a(m,n),sum
c
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+u(j,l)*s(l)*v(k,l)
            enddo ! l
c
            a(j,k) = sum
c
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
        subroutine fill(krank,m,n,a,u,v,s)
c
c       constructs singular vectors and singular values
c       for an m x n matrix a.
c
c       input:
c       krank -- length of s, and second dimension of u and v
c       m -- first dimension of u
c       n -- first dimension of v
c
c       output:
c       a -- u diag(s) v^T
c       u -- matrix of normalized left singular vectors
c       v -- matrix of normalized right singular vectors
c       s -- vector of singular values
c
        implicit none
        integer krank,m,n,mkrank,nkrank,j,k,l
        real*8 s(krank),r1,a(m,n),u(m,krank),v(n,krank),sum
c
        r1 = 1
c
c
c       Fill every entry of u and v
c       with i.i.d. random variables drawn uniformly from [-1,1].
c
        mkrank = m*krank
        call id_srand(mkrank,u)
c
        do l = 1,krank
          do j = 1,m
            u(j,l) = 2*u(j,l)-1
          enddo ! j
        enddo ! l
c
        nkrank = n*krank
        call id_srand(nkrank,v) 
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
c       Construct a = u diag(s) v^T.
c
        do l = 1,n
          do j = 1,m
c
            sum = 0
c
            do k = 1,krank
              sum = sum+u(j,k)*s(k)*v(l,k)
            enddo ! k
c
            a(j,l) = sum
c
          enddo ! j
        enddo ! l
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
        real*8 rms,a(m,n),prod
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
