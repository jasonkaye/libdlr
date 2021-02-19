c
c
c       dependencies: prini, idz_house
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,ifdisp
        real*8 work(len)
        complex*16 a(len),r(len),c(len),a0(len),v(len),
     1             q(len),q2(len),qadj(len),diff(len)
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
        call ccheck(ifdisp,m,n,a,a0,v,r,c,diff,q,q2,qadj,work)
c
c
        stop
        end
c
c
c
c
        subroutine ccheck(ifdisp,m,n,a,a0,v,r,c,diff,q,q2,qadj,work)
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,ind(len),k,j,indprod(len),ifdisp,ifadjoint,
     1          loop
        real*8 ss(len),eps,r1,pi,diffmax,
     1         errmax,errrms,errmaxadj,errrmsadj,
     2         diffrms,diffrmsadj,difffrob,difffrobadj,work(m)
        complex*16 a(m,n),q(m,m),q2(m,m),qadj(m,m),v(m),c(m,n),
     1             a0(m,n),r(m,n),cswap,diff(m,n),ci
c
        r1 = 1
        ci = (0,1)
        pi = 4*atan(r1)
c
c
c       Fill a with something.
c
        do k = 1,n
          do j = 1,m
            a(j,k) = sin(j*k/(r1*m+1)) - ci*cos(j**2*k/(r1*m+2))
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
        call idzp_qrpiv(eps,m,n,a,krank,ind,ss)
        if(ifdisp .eq. 1)
     1   call rectdisp('after idzp_qrpiv, a = *',a,2*m,n)
        call prinf('krank = *',krank,1)
        call prinf('ind = *',ind,krank)
c
c
        call idz_qinqr(m,n,a,krank,q)
        if(ifdisp .eq. 1) call rectdisp('q = *',q,2*m,m)
c
c
c       Check that q is unitary.
c
        call ccheckuni(m,m,q,errmax,errrms)
c
        call prin2('errmax = *',errmax,1)
        call prin2('errrms = *',errrms,1)
c
        call cadjoint(m,q)
        call ccheckuni(m,m,q,errmaxadj,errrmsadj)
        call cadjoint(m,q)
c
        call prin2('errmaxadj = *',errmaxadj,1)
        call prin2('errrmsadj = *',errrmsadj,1)
c
c
c       Apply the q in the QR factorization of a
c       to the identity matrix, using cqmatmat,
c       and check that it matches the q computed using cqinqr.
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
c       Apply q to q2.
c
        ifadjoint = 0
        call idz_qmatmat(ifadjoint,m,n,a,krank,m,q2,work)
        if(ifdisp .eq. 1) call rectdisp('q2 = *',q2,2*m,m)
c
        difffrob = 0
c
        do k = 1,m
          do j = 1,m
            difffrob = difffrob
     1               + (q(j,k)-q2(j,k)) * conjg(q(j,k)-q2(j,k))
          enddo ! j
        enddo ! k
c
        difffrob = sqrt(difffrob)
c
        call prin2('diffrob = *',difffrob,1)
c
c
c       Apply the adjoint of q in the QR factorization of a
c       to the identity matrix, using cqmatmat,
c       and check that it matches the adjoint
c       of the q computed using cqinqr.
c
c       Fill qadj with the identity matrix.
c
        do k = 1,m
          do j = 1,m
            qadj(j,k) = 0
          enddo ! k
        enddo ! j
c
        do k = 1,m
          qadj(k,k) = 1
        enddo ! k
c
c       Apply the adjoint of q to qadj.
c
        ifadjoint = 1
        call idz_qmatmat(ifadjoint,m,n,a,krank,m,qadj,work)
        if(ifdisp .eq. 1) call rectdisp('qadj = *',qadj,2*m,m)
c
        difffrobadj = 0
c
        do k = 1,m
          do j = 1,m
            difffrobadj = difffrobadj
     1                  + (q(k,j)-conjg(qadj(j,k)))
     2                  * (conjg(q(k,j))-qadj(j,k))
          enddo ! j
        enddo ! k
c
        difffrobadj = sqrt(difffrobadj)
c
        call prin2('diffrobadj = *',difffrobadj,1)
c
c
c       Apply the q in the QR factorization of a
c       to the first unit basis vector, using cqmatvec,
c       and check that it matches the first column of q
c       computed using cqinqr.
c
        v(1) = 1
        do k = 2,m
          v(k) = 0
        enddo ! k
c
        ifadjoint = 0
        call idz_qmatvec(ifadjoint,m,n,a,krank,v)
        call prin2('v = *',v,2*m)
c
        diffrms = 0
c
        do k = 1,m
          diffrms = diffrms + (q(k,1)-v(k)) * conjg(q(k,1)-v(k))
        enddo ! k
c
        diffrms = sqrt(diffrms)
c
        call prin2('diffrms = *',diffrms,1)
c
c
c       Apply the adjoint of the q in the QR factorization of a
c       to the first unit basis vector, using cqmatvec,
c       and check that it matches the conjugate of the first row of q
c       computed using cqinqr.
c
        v(1) = 1
        do k = 2,m
          v(k) = 0
        enddo ! k
c
        ifadjoint = 1
        call idz_qmatvec(ifadjoint,m,n,a,krank,v)
        call prin2('v = *',v,2*m)
c
        diffrmsadj = 0
c
        do k = 1,m
          diffrmsadj = diffrmsadj
     1               + (conjg(q(1,k))-v(k)) * (q(1,k)-conjg(v(k)))
        enddo ! k
c
        diffrmsadj = sqrt(diffrmsadj)
c
        call prin2('diffrmsadj = *',diffrmsadj,1)
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
            call idzr_qrpiv(m,n,a,krank,ind,ss)
            if(ifdisp .eq. 1)
     1       call rectdisp('after idzr_qrpiv, a = *',a,m,n)
c
c
            call idz_qinqr(m,n,a,krank,q)
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
          if(ifdisp .eq. 1) call rectdisp('r = *',r,2*m,n)
c
c         Multiply q and r to get c.
c
          call crectmult(m,m,q,n,r,c)
          if(ifdisp .eq. 1) call rectdisp('c = *',c,2*m,n)
          if(ifdisp .eq. 1) call rectdisp('a0 = *',a0,2*m,n)
c
c         Rearrange c according to ind.
c
          do k = krank,1,-1
            do j = 1,m
c
              cswap = c(j,k)
              c(j,k) = c(j,ind(k))
              c(j,ind(k)) = cswap
c
            enddo ! j
          enddo ! k
c
          if(ifdisp .eq. 1) call rectdisp('after rearrangement, c = *',
     1                                    c,2*m,n)
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
          if(ifdisp .eq. 1) call rectdisp('diff = *',diff,2*m,n)
c
          call prin2('diffmax = *',diffmax,1)
c
c
        enddo ! loop
c
c
c       Form the product of the permutations in ind.
c
        call idz_permmult(krank,ind,n,indprod)
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
        subroutine crectmult(l,m,a,n,b,c)
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
        complex*16 a(l,m),b(m,n),c(l,n)
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
        subroutine cmatcopy(n,a,b)
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
        complex*16 a(n,n),b(n,n)
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
        subroutine cadjoint(n,a)
c
c       conjugate transposes the n x n matrix a in place.
c
c       input:
c       n -- size of a
c       a -- matrix to be conjugate transposed
c
c       output:
c       a -- adjoint of the matrix that was input
c            via the same array a
c
        implicit none
        integer n,j,k
        complex*16 a(n,n),cswap
c
c
        do j = 1,n
          do k = 1,j
c
c           Swap a(k,j) and a(j,k), and conjugate them.
c
            cswap = a(k,j)
            a(k,j) = conjg(a(j,k))
            a(j,k) = conjg(cswap)
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
        subroutine ccheckuni(n,m,a,errmax,errrms)
c
c       calculates the relative maximum and root-mean-square errors
c       corresponding to how much the product of a and its adjoint
c       differs from the identity.
c
c       input:
c       n -- first dimension of a
c       m -- second dimension of a
c       a -- matrix whose unitarity will be checked
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
        real*8 errmax,errrms,diff,r1
        complex*16 a(n,m),prod
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
              prod = prod+a(k,j)*conjg(a(k,l))
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
