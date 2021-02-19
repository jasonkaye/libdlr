c
c
c       dependencies: prini, idd_house, idd_qrpiv, idd_id, lapack.a,
c                     blas.a, and (for the debugging code) id_rand
c
c
        implicit none
c
        integer len
        parameter(len = 5 000 000)
c
        character*1 jobz
        integer m,n,krank,list(len),k,its,iterations,
     1          lda,lwork,info,ldu,ldvt,ier
        real*8 a(len),proj(len),eps,rnorms(len),
     1         b(len),a2(len),u(len),t(len),specnorm(len),
     2         work(len),v(len),s(len)
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
c
c       Fill the matrix a.
c
        call fill(m,n,a)
c
c
c       Copy a into a2.
c
        do k = 1,m*n
          a2(k) = a(k)
        enddo ! k
c
c
c       Use LAPACK to SVD a2.
c
        jobz = 'S'
        lda = m
        lwork = len
        ldu = m
        ldvt = n
c
        call dgesdd(jobz,m,n,a2,lda,s,u,ldu,v,ldvt,
     1              work,lwork,list,info)
c
        call prinf('info = *',info,1)
c
c
c       Determine how many entries in s are greater than eps
c       and display them.
c
        eps = .1d-12
        krank = 0
c
        do k = 1,n
          if(s(k) .gt. eps) krank = krank+1
        enddo ! k
c
        call prin2('s = *',s,krank)
c
c
c       ID a.
c
        do k = 1,m*n
          proj(k) = a(k)
        enddo ! k
c
        eps = .1d-12
        call iddp_id(eps,m,n,proj,krank,list,rnorms)
        call prinf('krank = *',krank,1)
c
c
c       Copy the selected columns of a into b
c       (in the order given by list).
c
        call idd_copycols(m,n,a,krank,list,b)
c
c
c       Convert the approximation to a in the form of an ID
c       to an approximation in the form of an SVD.
c
        call idd_id2svd(m,krank,b,n,list,proj,u,v,s,ier,work)
        call prinf('ier = *',ier,1)
        call prin2('s = *',s,krank)
c
c
c       Form a2 = u diag(s) v^T and compare it to a.
c
        call reconsvd(m,krank,u,s,n,v,a2)
c
        eps = .1d-3
        call id_srand(n,work)
        iterations = 100
c
        call powerchk(eps,m,n,a,a2,work,iterations,its,specnorm,t)
        call prin2('specnorm(its) = *',specnorm(its),1)
c
c
        stop
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
c       n -- second dimension of a and first dimension of v
c       v -- rightmost matrix in the product a = u diag(s) v^T
c
c       output:
c       a -- matrix product u diag(s) v^T
c
        implicit none
        integer m,n,krank,j,k,l
        real*8 u(m,krank),s(krank),v(n,krank),a(m,n),sum
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
        subroutine fill(m,n,a)
c
c       fills the matrix a.
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a
c
c       output:
c       a -- filled matrix
c
        implicit none
        integer m,n,j,k
        real*8 a(m,n),r1
c
        r1 = 1
c
c
        do k = 1,n
          do j = 1,m
            a(j,k) = 1/(r1+j+k)
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
        subroutine powerchk(eps,m,n,a,b,u,iterations,its,specnorm,t)
c
c       estimates the spectral norm || a - b ||,
c       using the power method.
c
c       input:
c       eps -- desired accuracy of the norm estimate
c       m -- first dimension of a and b
c       n -- second dimension of a and b
c       u -- vector to which a-b is to be applied;
c            u is destroyed by this routine
c       iterations -- length of specnorm; number of power method
c                     iterations to conduct
c
c       output:
c       its -- number of iterations actually conducted
c              (the lesser of iterations and the number required
c               to attain convergence to within eps)
c       specnorm -- estimates of the spectral norm || a - b ||,
c                   as provided by power method iterations
c
c       work:
c       t -- must be at least m real*8 elements long
c
        implicit none
        integer m,n,ktest,iteration,j,k,iterations,its
        real*8 a(m,n),b(m,n),eps,specnorm(iterations),
     1         u(n),t(m),rss1,rss2
c
c
        its = iterations
c
        do iteration = 1,iterations+1
c
c
c         Note that the first iteration normalizes u,
c         but that rss1 is meaningless during the first iteration,
c         unless the input u was normalized.
c
c
          ktest = 1
          do j = 1,m
c
            t(j) = 0
c
            do k = 1,n
              t(j) = t(j)+a(j,k)*u(k)-b(j,k)*u(k)
            enddo ! j
c
          enddo ! k
c
c
          rss1 = 0
c
          do j = 1,m
            rss1 = rss1 + (t(j))**2
          enddo ! j
c
          rss1 = sqrt(rss1)
c
c
          if(rss1 .gt. 0) then
            do j = 1,m
              t(j) = t(j) / rss1
            enddo ! j
          endif
c
c
          ktest = 1
          do k = 1,n
c
            u(k) = 0
c
            do j = 1,m
              u(k) = u(k)+a(j,k)*t(j)-b(j,k)*t(j)
            enddo ! j
c
          enddo ! k
c
c
          rss2 = 0
c
          do k = 1,n
            rss2 = rss2 + (u(k))**2
          enddo ! k
c
          rss2 = sqrt(rss2)
c
c
          if(rss2 .gt. 0) then
            do k = 1,n
              u(k) = u(k) / rss2
            enddo ! k
          endif
c
c
          if(iteration .gt. 1) then
            specnorm(iteration-1) = sqrt(rss1*rss2)
          endif
c
          if(iteration .gt. 2) then
            if(abs(specnorm(iteration-2)-specnorm(iteration-1))
     1       .le. eps*specnorm(iteration-1)) then
             its = iteration-1
             goto 1000
            endif
          endif
c
c
        enddo ! iteration
c
 1000   continue
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
