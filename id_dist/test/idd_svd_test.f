c
c
c       dependencies: prini, idd_house, idd_qrpiv, lapack.a, blas.a,
c                     and (for the debugging code) id_rand
c
c
        implicit none
c
        integer len
        parameter(len = 4 000 000)
c
        integer m,n,k,krank,ier,ifdisp,iterations,its,n2,iu,iv,is,
     1          lwork,loop
        real*8 specnorm(10 000),eps,s(len),a(len),work(len),
     1         v(len),u(len),a2(len),a0(len)
c
c
        ifdisp = 0
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
c       Fill the matrix to be SVD'd.
c
        krank = 11
        call fill(krank,m,n,a0,s)
        call prin2('s = *',s,krank+1)
c
        if(ifdisp .eq. 1) then
c
          call prin2('a0 = *',a0,0)
          do k = 1,n
            call prin2('*',a0(n*(k-1)+1),2*m)
          enddo ! k
c
        endif
c
c
        do loop = 1,2
c
c
          do k = 1,m*n
            a(k) = a0(k)
          enddo ! k
c
c
          if(loop .eq. 1) then
c
c           SVD a.
c
            call iddr_svd(m,n,a,krank,u,v,s,ier,work)
            call prin2('s = *',s,krank)
c
          endif ! loop .eq. 1
c
c
          if(loop .eq. 2) then
c
c           SVD a.
c
            eps = .1d-8
            lwork = len
c
            call iddp_svd(lwork,eps,m,n,a,krank,iu,iv,is,work,ier)
c
c           Copy u, v, and s from work.
c
            do k = 1,krank*m
              u(k) = work(iu+k-1)
            enddo ! k
c
            do k = 1,krank*n
              v(k) = work(iv+k-1)
            enddo ! k
c
            do k = 1,krank
              s(k) = work(is+k-1)
            enddo ! k
            call prin2('s = *',s,krank)
c
          endif ! loop .eq. 2
c
c
c         Form a2 = U S V^T.
c
          call svdrecon(m,krank,u,s,n,v,a2)
c
          if(ifdisp .eq. 1) then
c
            call prin2('a2 = *',a2,0)
            do k = 1,n
              call prin2('*',a2(n*(k-1)+1),2*m)
            enddo ! k
c
          endif
c
c
c         Calculate the spectral norm of a0-a2.
c
          n2 = 2*n
          call id_srand(n2,u)
c
          eps = .1d-1
          iterations = 100
          call powerchk(eps,m,n,a0,a2,u,iterations,
     1                  its,specnorm,work)
c
          call prin2('specnorm(its) = *',specnorm(its),1)
c
c
        enddo ! loop
c
c
        stop
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
        real*8 eps,specnorm(iterations),rss1,rss2,a(m,n),b(m,n),
     1         u(n),t(m)
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
            rss1 = rss1 + t(j)**2
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
            rss2 = rss2 + u(k)**2
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
        subroutine svdrecon(m,krank,u,s,n,v,a)
c
c       forms a = u diag(s) v^T.
c
c       input:
c       m -- first dimension of u and a
c       krank -- size of s, second dimension of u,
c                and second dimension of v
c       u -- leftmost matrix in the product a = u diag(s) v^T
c       s -- entries on the diagonal in the middle matrix
c            in the product a = u diag(s) v^T
c       n -- second dimension of v^T and a
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
        subroutine fill(krank,m,n,a,s)
c
c       fills an m x n matrix with suitably decaying singular values,
c       and left and right singular vectors taken from the DCT-IV.
c
c       input:
c       krank -- one less than the rank of the matrix to be constructed
c       m -- first dimension of a
c       n -- second dimension of a
c
c       output:
c       a -- filled matrix
c       s -- singular values of a
c
        implicit none
        integer krank,j,k,l,m,n
        real*8 r1,pi,a(m,n),sum,s(krank+1)
c
        r1 = 1
        pi = 4*atan(r1)
c
c
c       Specify the singular values.
c
        do k = 1,krank
          s(k) = exp(log(1d-10)*(k-1)/(krank-1))
        enddo ! k
c
        s(krank+1) = 1d-10
c
c
c       Construct a.
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+cos(pi*(j-r1/2)*(l-r1/2)/m)*sqrt(r1*2/m)
     1                 *cos(pi*(k-r1/2)*(l-r1/2)/n)*sqrt(r1*2/n)*s(l)
            enddo ! l
c
            l = krank+1
            sum = sum+cos(pi*(j-r1/2)*(l-r1/2)/m)*sqrt(r1*2/m)
     1               *cos(pi*(k-r1/2)*(l-r1/2)/n)*sqrt(r1*2/n)*s(l)
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
