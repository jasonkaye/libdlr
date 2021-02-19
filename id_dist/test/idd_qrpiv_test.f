c
c
c       dependencies: prini, idd_house
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,ifdisp
        real*8 a(len),r(len),c(len),diff(len),a0(len),v(len),
     1         q(len),q2(len),qtr(len)
c
c
        call prini(6,13)
c
c
        print *,
     1   'To display full matrices, enter 1; otherwise, enter 0:'
        read *,ifdisp
        call prinf('ifdisp = *',ifdisp,1)
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
        call check(ifdisp,m,n,a,a0,v,r,c,diff,q,q2,qtr)
c
c
        stop
        end
c
c
c
c
        subroutine check(ifdisp,m,n,a,a0,v,r,c,diff,q,q2,qtr)
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,ind(len),k,j,indprod(len),ifdisp,
     1          iftranspose,loop
        real*8 a(m,n),ss(len),eps,r1,q(m,m),q2(m,m),qtr(m,m),v(m),
     1         c(m,n),a0(m,n),r(m,n),pi,diff(m,n),diffmax,
     2         errmax,errrms,errmaxtr,errrmstr,rswap,
     3         diffrms,diffrmstr,difffrob,difffrobtr
c
        r1 = 1
        pi = 4*atan(r1)
c
c
c       Fill a with something.
c
        do k = 1,n
          do j = 1,m
            a(j,k) = sin(j*k/(r1*m+1))
          enddo ! j
        enddo ! k
c
        if(n .ge. 6) then
c
          do k = 4,6
            do j = 1,m
              a(j,k) = ( a(j,k-3)+a(j,1) )/5
            enddo ! j        
          enddo ! k
c
        endif
c
c
c       Duplicate a into a0.
c
        do k = 1,n
          do j = 1,m
            a0(j,k) = a(j,k)
          enddo ! j
        enddo ! k
c
c
c       Compute the pivoted QR factorization of a.
c
        eps = .1d-13
        call iddp_qrpiv(eps,m,n,a,krank,ind,ss)
        if(ifdisp .eq. 1)
     1   call rectdisp('after iddp_qrpiv, a = *',a,m,n)
        call prinf('krank = *',krank,1)
        call prinf('ind = *',ind,krank)
c
c
        call idd_qinqr(m,n,a,krank,q)
        if(ifdisp .eq. 1) call rectdisp('q = *',q,m,m)
c
c
c       Check that q is orthogonal.
c
        call checkuni(m,m,q,errmax,errrms)
c
        call prin2('errmax = *',errmax,1)
        call prin2('errrms = *',errrms,1)
c
        call transpose(m,q)
        call checkuni(m,m,q,errmaxtr,errrmstr)
        call transpose(m,q)
c
        call prin2('errmaxtr = *',errmaxtr,1)
        call prin2('errrmstr = *',errrmstr,1)
c
c
c       Apply the q in the QR factorization of a
c       to the identity matrix, using qmatmat,
c       and check that it matches the q computed using qinqr.
c
c       Fill q2 with the identity matrix.
c
        do k = 1,m
          do j = 1,m
            q2(j,k) = 0
          enddo ! k
        enddo ! j
c
        do k = 1,m
          q2(k,k) = 1
        enddo ! k
c
c       apply q to q2
c
        iftranspose = 0
        call idd_qmatmat(iftranspose,m,n,a,krank,m,q2,v)
        if(ifdisp .eq. 1) call rectdisp('q2 = *',q2,m,m)
c
        difffrob = 0
c
        do k = 1,m
          do j = 1,m
            difffrob = difffrob + (q(j,k)-q2(j,k))**2
          enddo ! j
        enddo ! k
c
        difffrob = sqrt(difffrob)
c
        call prin2('diffrob = *',difffrob,1)
c
c
c       Apply the transpose of q in the QR factorization of a
c       to the identity matrix, using qmatmat,
c       and check that it matches the transpose
c       of the q computed using qinqr.
c
c       Fill qtr with the identity matrix.
c
        do k = 1,m
          do j = 1,m
            qtr(j,k) = 0
          enddo ! k
        enddo ! j
c
        do k = 1,m
          qtr(k,k) = 1
        enddo ! k
c
c       Apply the transpose of q to qtr.
c
        iftranspose = 1
        call idd_qmatmat(iftranspose,m,n,a,krank,m,qtr,v)
        if(ifdisp .eq. 1) call rectdisp('qtr = *',qtr,m,m)
c
        difffrobtr = 0
c
        do k = 1,m
          do j = 1,m
            difffrobtr = difffrobtr + (q(k,j)-qtr(j,k))**2
          enddo ! j
        enddo ! k
c
        difffrobtr = sqrt(difffrobtr)
c
        call prin2('diffrobtr = *',difffrobtr,1)
c
c
c       Apply the q in the QR factorization of a
c       to the first unit basis vector, using qmatvec,
c       and check that it matches the first column of q
c       computed using qinqr.
c
        v(1) = 1
        do k = 2,m
          v(k) = 0
        enddo ! k
c
        iftranspose = 0
        call idd_qmatvec(iftranspose,m,n,a,krank,v)
        call prin2('v = *',v,m)
c
        diffrms = 0
c
        do k = 1,m
          diffrms = diffrms + (q(k,1)-v(k))**2
        enddo ! k
c
        diffrms = sqrt(diffrms)
c
        call prin2('diffrms = *',diffrms,1)
c
c
c       Apply the transpose of the q in the QR factorization of a
c       to the first unit basis vector, using qmatvec,
c       and check that it matches the first row of q
c       computed using qinqr.
c
        v(1) = 1
        do k = 2,m
          v(k) = 0
        enddo ! k
c
        iftranspose = 1
        call idd_qmatvec(iftranspose,m,n,a,krank,v)
        call prin2('v = *',v,m)
c
        diffrmstr = 0
c
        do k = 1,m
          diffrmstr = diffrmstr + (q(1,k)-v(k))**2
        enddo ! k
c
        diffrmstr = sqrt(diffrmstr)
c
        call prin2('diffrmstr = *',diffrmstr,1)
c
c
        do loop = 1,2
c
c
          if(loop .eq. 2) then
c
c
c           Duplicate a0 into a.
c
            do k = 1,n
              do j = 1,m
                a(j,k) = a0(j,k)
              enddo ! j
            enddo ! k
c
c
c           Compute the pivoted QR factorization of a.
c
            call iddr_qrpiv(m,n,a,krank,ind,ss)
            if(ifdisp .eq. 1)
     1       call rectdisp('after iddr_qrpiv, a = *',a,m,n)
c
c
            call idd_qinqr(m,n,a,krank,q)
            if(ifdisp .eq. 1) call rectdisp('q = *',q,m,m)
c
c
          endif
c
c
c         Check that the qr factorization of a0
c         given by q and the triangular factor stored in a
c         correctly reconstructs a0.
c
c         Copy a into r and zero out the appropriate Householder
c         vectors that are stored in one triangle of a.
c
          do k = 1,n
            do j = 1,m
              r(j,k) = a(j,k)
            enddo ! j
          enddo ! k
c
          do k = 1,n
            if(k .lt. m) then
              do j = k+1,m
                r(j,k) = 0
              enddo ! j
            endif
          enddo ! k
c
          if(ifdisp .eq. 1) call rectdisp('r = *',r,m,n)
c
c         Multiply q and r to get c.
c
          call rectmult(m,m,q,n,r,c)
          if(ifdisp .eq. 1) call rectdisp('c = *',c,m,n)
          if(ifdisp .eq. 1) call rectdisp('a0 = *',a0,m,n)
c
c         Rearrange c according to ind.
c
          do k = krank,1,-1
            do j = 1,m
c
              rswap = c(j,k)
              c(j,k) = c(j,ind(k))
              c(j,ind(k)) = rswap
c
            enddo ! j
          enddo ! k
c
          if(ifdisp .eq. 1) call rectdisp('after rearrangement, c = *',
     1                                    c,m,n)
c
c         Display the difference of c from a0.
c
          diffmax = 0
c
          do k = 1,n
            do j = 1,m
              diff(j,k) = c(j,k)-a0(j,k)
              if(abs(diff(j,k)) .gt. diffmax) diffmax = abs(diff(j,k))
            enddo ! j
          enddo ! k
c
          if(ifdisp .eq. 1) call rectdisp('diff = *',diff,m,n)
c
          call prin2('diffmax = *',diffmax,1)
c
c
        enddo ! loop
c
c
c       Form the product of the permutations in ind.
c
        call idd_permmult(krank,ind,n,indprod)
        call prinf('indprod = *',indprod,n)
c
c
        return
        end
c
c
c
c
        subroutine rectdisp(str,a,m,n)
c
c       displays a real rectangular matrix a via prini,
c       with the first index of a ascending as you read the rows
c       from left to right,
c       and the second index of a ascending as you read the columns
c       from top to bottom.
c
c       input:
c       str -- message for prin2
c       a -- matrix to display
c       m -- first dimension of a
c       n -- second dimension of a
c
c       _N.B._: You must call prini for initialization
c               before calling this routine.
c
        implicit none
        integer m,n,k
        real*8 a(m,n)
        character*1 str(1)
c
c
        call prin2(str,a,0)
        do k = 1,n
          call prin2('*',a(1,k),m)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine rectmult(l,m,a,n,b,c)
c
c       multiplies a and b to get c.
c
c       input:
c       l -- first dimensions of a and c
c       m -- second dimension of a; first dimension of b
c       a -- matrix to multiply
c       n -- second dimension of b
c       b -- matrix to multiply
c
c       output:
c       c -- product of a and b
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
              c(i,k) = c(i,k)+a(i,j)*b(j,k)
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
        subroutine matcopy(n,a,b)
c
c       copies the n x n matrix a into b.
c
c       input:
c       n -- size of a and b
c       a -- matrix to duplicate
c
c       output:
c       b -- duplicate of a
c
        implicit none
        integer n,j,k
        real*8 a(n,n),b(n,n)
c
c
        do j = 1,n
          do k = 1,n
            b(k,j) = a(k,j)
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
        subroutine transpose(n,a)
c
c       transposes the n x n matrix a in place.
c
c       input:
c       n -- size of a
c       a -- matrix to be transposed
c
c       output:
c       a -- transposition of the matrix that was input
c            via the same array a
c
        implicit none
        integer n,j,k
        real*8 a(n,n),rswap
c
c
        do j = 1,n
          do k = 1,j
c
c           Swap a(k,j) and a(j,k).
c
            rswap = a(k,j)
            a(k,j) = a(j,k)
            a(j,k) = rswap
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
        subroutine checkuni(n,m,a,errmax,errrms)
c
c       calculates the relative maximum and root-mean-square errors
c       corresponding to how much the product of a and its adjoint
c       differs from the identity.
c
c       input:
c       n -- first dimension of a
c       m -- second dimension of a
c       a -- matrix whose orthogonality will be checked
c
c       output:
c       errmax -- ratio of the maximum error
c                 in the product of a and its adjoint from the identity
c                 to the maximum among all elements in the identity
c                 (a maximum which is, of course, equal to 1)
c       errrms -- ratio of the root-mean-square error
c                 in the product of a and its adjoint from the identity
c                 to the root-mean-square value of all elements
c                 in the identity
c
        implicit none
        integer m,n,j,k,l
        real*8 a(n,m),errmax,errrms,diff,r1,prod
c
        r1 = 1
c
c
        errmax = 0
        errrms = 0
c
        do j = 1,m
          do l = 1,m
c
            prod = 0
            do k = 1,n
              prod = prod+a(k,j)*a(k,l)
            enddo ! k
c
            if(j .eq. l) diff = abs(prod-1)
            if(j .ne. l) diff = abs(prod)
c
            if(diff .gt. errmax) errmax = diff
            errrms = errrms+diff**2
c
          enddo ! l
        enddo ! j
c
        errmax = errmax/1
        errrms = errrms/(r1*m)
        errrms = sqrt(errrms)
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
