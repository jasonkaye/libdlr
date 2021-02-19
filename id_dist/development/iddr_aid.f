c
c
c       dependencies: prini, idd_house, idd_qrpiv, idd_id, id_rand,
c                     idd_sfft, id_rtrans, idd_frm, dfft
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,list(len)
        real*8 a(len),b(len),proj(len),col(len),errmax,errrms
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
        krank = 5
        call prinf('krank = *',krank,1)
c
c
c       Fill a.
c
        call fill(krank,m,n,a)
        call prin2('a = *',a,m*n)
c
c
c       Initialize b as the work array for iddr_aid.
c
        call iddr_aidi(m,n,krank,b)
c
c
c       ID a.
c
        call iddr_aid(m,n,a,krank,b,list,proj)
        call prinf('list = *',list,krank)
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idd_copycols(m,n,a,krank,list,col)
c
c
c       Reconstruct a, obtaining b.
c
        call idd_reconid(m,krank,col,n,list,proj,b)
c
c
c       Compute the difference between a and b.
c
        call materr(m,n,a,b,errmax,errrms)
        call prin2('errmax = *',errmax,1)
        call prin2('errrms = *',errrms,1)
c
c
        stop
        end
c
c
c
c
        subroutine fill(krank,m,n,a)
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
c
        implicit none
        integer krank,j,k,l,m,n
        real*8 r1,pi,a(m,n),sum
c
        r1 = 1
        pi = 4*atan(r1)
c
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+cos(pi*(j-r1/2)*(l-r1/2)/m)*sqrt(r1*2/m)
     1                 *cos(pi*(k-r1/2)*(l-r1/2)/n)*sqrt(r1*2/n)
     2                 *exp(log(1d-10)*(l-1)/(krank-1))
            enddo ! l
c
            l = krank+1
            sum = sum+cos(pi*(j-r1/2)*(l-r1/2)/m)*sqrt(r1*2/m)
     1               *cos(pi*(k-r1/2)*(l-r1/2)/n)*sqrt(r1*2/n)*1d-10
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
        real*8 a(m,n),b(m,n),errmax,errrms,diff,amax,arss
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine iddr_aid computes the ID, to a specified rank,
c       of an arbitrary matrix. This routine is randomized.
c
c       routine iddr_aidi initializes routine iddr_aid.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddr_aid(m,n,a,krank,w,list,proj)
c
c       computes the ID of the matrix a, i.e., lists in list
c       the indices of krank columns of a such that 
c
c       a(j,list(k))  =  a(j,list(k))
c
c       for all j = 1, ..., m; k = 1, ..., krank, and
c
c                       min(m,n,krank)
c       a(j,list(k))  =     Sigma      a(j,list(l)) * proj(l,k-krank)(*)
c                            l=1
c
c                     +  epsilon(j,k-krank)
c
c       for all j = 1, ..., m; k = krank+1, ..., n,
c
c       for some matrix epsilon, dimensioned epsilon(m,n-krank),
c       whose norm is (hopefully) minimized by the pivoting procedure.
c
c       input:
c       m -- number of rows in a
c       n -- number of columns in a
c       a -- matrix to be ID'd; the present routine does not alter a
c       krank -- rank of the ID to be constructed
c       w -- initialization array that routine iddr_aidi
c            has constructed
c
c       output:
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd
c
c       _N.B._: The algorithm used by this routine is randomized.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,n,krank,list(n),lw,ir,lr,lw2,iw
        real*8 a(m,n),proj(krank*(n-krank)),w((2*krank+17)*n+27*m+100)
c
c
c       Allocate memory in w.
c
        lw = 0
c
        iw = lw+1
        lw2 = 27*m+100+n
        lw = lw+lw2
c
        ir = lw+1
        lr = (krank+8)*2*n
        lw = lw+lr
c
c
        call iddr_aid0(m,n,a,krank,w(iw),list,proj,w(ir))
c
c
        return
        end
c
c
c
c
        subroutine iddr_aid0(m,n,a,krank,w,list,proj,r)
c
c       routine iddr_aid serves as a memory wrapper
c       for the present routine
c       (see iddr_aid for further documentation).
c
        implicit none
        integer k,l,m,n2,n,krank,list(n),mn,lproj
        real*8 a(m,n),r(krank+8,2*n),proj(krank,n-krank),
     1         w(27*m+100+n)
c
c       Please note that the second dimension of r is 2*n
c       (instead of n) so that if krank+8 >= m/2, then
c       we can copy the whole of a into r.
c
c
c       Retrieve the number of random test vectors
c       and the greatest integer less than m that is
c       a positive integer power of two.
c
        l = w(1)
        n2 = w(2)
c
c
        if(l .lt. n2 .and. l .le. m) then
c
c         Apply the random matrix.
c
          do k = 1,n
            call idd_sfrm(l,m,n2,w(11),a(1,k),r(1,k))
          enddo ! k
c
c         ID r.
c
          call iddr_id(l,n,r,krank,list,w(26*m+101))
c
c         Retrieve proj from r.
c
          lproj = krank*(n-krank)
          call iddr_copydarr(lproj,r,proj)
c
        endif
c
c
        if(l .ge. n2 .or. l .gt. m) then
c
c         ID a directly.
c
          mn = m*n
          call iddr_copydarr(mn,a,r)
          call iddr_id(m,n,r,krank,list,w(26*m+101))
c
c         Retrieve proj from r.
c
          lproj = krank*(n-krank)
          call iddr_copydarr(lproj,r,proj)
c
        endif
c
c
        return
        end
c
c
c
c
        subroutine iddr_copydarr(n,a,b)
c
c       copies a into b.
c
c       input:
c       n -- length of a and b
c       a -- array to copy into b
c
c       output:
c       b -- copy of a
c
        implicit none
        integer n,k
        real*8 a(n),b(n)
c
c
        do k = 1,n
          b(k) = a(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine iddr_aidi(m,n,krank,w)
c
c       initializes the array w for using routine iddr_aid.
c
c       input:
c       m -- number of rows in the matrix to be ID'd
c       n -- number of columns in the matrix to be ID'd
c       krank -- rank of the ID to be constructed
c
c       output:
c       w -- initialization array for using routine iddr_aid
c
        implicit none
        integer m,n,krank,l,n2
        real*8 w((2*krank+17)*n+27*m+100)
c
c
c       Set the number of random test vectors to 8 more than the rank.
c
        l = krank+8
        w(1) = l
c
c
c       Initialize the rest of the array w.
c
        n2 = 0
        if(l .le. m) call idd_sfrmi(l,m,n2,w(11))
        w(2) = n2
c
c
        return
        end
