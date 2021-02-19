c
c
c       dependencies: prini, idz_house, idz_qrpiv
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,ifdisp
        complex*16 a(len),a0(len),col(len)
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
        call ccheck(ifdisp,m,n,a,a0,col)
c
c
        stop
        end
c
c
c
c
        subroutine ccheck(ifdisp,m,n,a,a0,col)
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,j,k,krank,list(len),ifdisp,loop
        real*8 r1,pi,work(len),errmax,errrms,eps
        complex*16 a(m,n),a0(m,n),approx(len),col(len),ci
c
        r1 = 1
        ci = (0,1)
        pi = 4*atan(r1)
c
c
c       Fill a0 with something.
c
        do k = 1,n
          do j = 1,m
            a0(j,k) = sin(j*k/(r1*m+1)) - ci*cos(j**2*k/(r1*m+3))
          enddo ! j
        enddo ! k
c
        if(n .ge. 6) then
c
          do k = 4,6
            do j = 1,m
              a0(j,k) = ( a0(j,k-3)+a0(j,1) )/5
            enddo ! j
          enddo ! k
c
        endif
c
        if(ifdisp .eq. 1) call rectdisp('a0 = *',a0,2*m,n)
c
c
        do loop = 1,2
c
c
c         Duplicate a0 into a.
c
          do k = 1,n
            do j = 1,m
              a(j,k) = a0(j,k)
            enddo ! j
          enddo ! k
c
c
          if(loop .eq. 1) then
c
c
c           ID a.
c
            eps = .1d-13
c
            call idzp_id(eps,m,n,a,krank,list,work)
c
            call prinf('krank = *',krank,1)
            call prinf('list = *',list,n)
            if(ifdisp .eq. 1)
     1       call rectdisp('a (proj) = *',a,2*krank,n-krank)
c
c
          endif ! loop .eq. 1
c
c
          if(loop .eq. 2) then
c
c
c           ID a.
c
            call idzr_id(m,n,a,krank,list,work)
            call prinf('list = *',list,n)
            if(ifdisp .eq. 1)
     1       call rectdisp('a (proj) = *',a,2*krank,n-krank)
c
c
          endif ! loop .eq. 2
c
c
c         Copy the selected columns of a0 into col
c         (in the order given by list).
c
          call idz_copycols(m,n,a0,krank,list,col)
c
c
c         Reconstruct a0 from col and the proj in a.
c
          call idz_reconid(m,krank,col,n,list,a,approx)
          if(ifdisp .eq. 1) call rectdisp('approx = *',approx,2*m,n)
c
c
          if(krank .gt. 0) then
c
c           Calculate the relative maximum and root-mean-square errors
c           corresponding to how much a0 and approx differ.
c
            call cmaterr(m,n,a0,approx,errmax,errrms)
            call prin2('errmax = *',errmax,1)
            call prin2('errrms = *',errrms,1)
c
          endif
c
c
        enddo ! loop
c
c
        return
        end
c
c
c
c
        subroutine cmaterr(m,n,a,b,errmax,errrms)
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
        real*8 errmax,errrms,diff,amax,arss
        complex*16 a(m,n),b(m,n)
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
            arss = arss+a(j,k)*conjg(a(j,k))
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
