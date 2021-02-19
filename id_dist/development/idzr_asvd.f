c
c
c       dependencies: prini, idz_house, idz_qrpiv, idz_id, id_rand,
c                     idz_sfft, id_rtrans, idz_frm, idzr_aid,
c                     idz_id2svd, dfft, lapack.a, blas.a
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,ier
        real*8 errmax,errrms,s(100 000)
        complex*16 a(len),b(len),u(len),v(len),w(len)
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
        call fill(krank,m,n,a,s)
        call prin2('a = *',a,m*n)
        call prin2('s = *',s,krank+1)
c
c
c       Calculate an SVD approximating a.
c
        call idzr_aidi(m,n,krank,w)
        call idzr_asvd(m,n,a,krank,w,u,v,s,ier)
c
c
c       Construct b = u diag(s) v^*.
c
        call reconsvd(m,krank,u,s,n,v,b)
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
        subroutine fill(krank,m,n,a,s)
c
c       fills an m x n matrix with suitably decaying singular values,
c       and left and right singular vectors taken from the DFT.
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
        real*8 r1,pi,s(krank+1)
        complex*16 a(m,n),sum,ci
c
        r1 = 1
        pi = 4*atan(r1)
        ci = (0,1)
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
              sum = sum+exp(2*pi*ci*(j-r1)*(l-r1)/m)*sqrt(r1/m)
     1                 *exp(2*pi*ci*(k-r1)*(l-r1)/n)*sqrt(r1/n)*s(l)
            enddo ! l
c
            l = krank+1
            sum = sum+exp(2*pi*ci*(j-r1)*(l-r1)/m)*sqrt(r1/m)
     1               *exp(2*pi*ci*(k-r1)*(l-r1)/n)*sqrt(r1/n)*s(l)
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
        subroutine reconsvd(m,krank,u,s,n,v,a)
c
c       forms a = u diag(s) v^*.
c
c       input:
c       m -- first dimension of u and a
c       krank -- size of s, and second dimension of u and v
c       u -- leftmost matrix in the product a = u diag(s) v^*
c       s -- entries on the diagonal in the middle matrix
c            in the product a = u diag(s) v^*
c       n -- second dimension of a and first dimension of v
c       v -- rightmost matrix in the product a = u diag(s) v^*
c
c       output:
c       a -- matrix product u diag(s) v^*
c
        implicit none
        integer m,n,krank,j,k,l
        real*8 s(krank)
        complex*16 u(m,krank),v(n,krank),a(m,n),sum
c
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+u(j,l)*s(l)*conjg(v(k,l))
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine idzr_aid computes the SVD, to a specified rank,
c       of an arbitrary matrix. This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idzr_asvd(m,n,a,krank,w,u,v,s,ier)
c
c       constructs a rank-krank SVD  u diag(s) v^*  approximating a,
c       where u is an m x krank matrix whose columns are orthonormal,
c       v is an n x krank matrix whose columns are orthonormal,
c       and diag(s) is a diagonal krank x krank matrix whose entries
c       are all nonnegative. This routine uses a randomized algorithm.
c
c       input:
c       m -- number of rows in a
c       n -- number of columns in a 
c       a -- matrix to be decomposed; the present routine does not
c            alter a
c       krank -- rank of the SVD being constructed
c       w -- initialization array that routine idzr_aidi
c            has constructed (for use in the present routine,
c            w must be at least
c            (2*krank+22)*m+(6*krank+21)*n+8*krank**2+10*krank+90
c            complex*16 elements long)
c
c       output:
c       u -- matrix of orthonormal left singular vectors of a
c       v -- matrix of orthonormal right singular vectors of a
c       s -- array of singular values of a
c       ier -- 0 when the routine terminates successfully;
c              nonzero otherwise
c
c       _N.B._: The algorithm used by this routine is randomized.
c
        implicit none
        integer m,n,krank,lw,ilist,llist,iproj,lproj,icol,lcol,
     1          iwork,lwork,iwinit,lwinit,ier
        real*8 s(krank)
        complex*16 a(m,n),u(m,krank),v(n,krank),
     1             w((2*krank+22)*m+(6*krank+21)*n+8*krank**2
     2              +10*krank+90)
c
c
c       Allocate memory in w.
c
        lw = 0
c
        iwinit = lw+1
        lwinit = (2*krank+17)*n+21*m+80
        lw = lw+lwinit
c
        ilist = lw+1
        llist = n
        lw = lw+llist
c
        iproj = lw+1
        lproj = krank*(n-krank)
        lw = lw+lproj
c
        icol = lw+1
        lcol = m*krank
        lw = lw+lcol
c
        iwork = lw+1
        lwork = (krank+1)*(m+3*n+10)+9*krank**2
        lw = lw+lwork
c
c
        call idzr_asvd0(m,n,a,krank,w(iwinit),u,v,s,ier,
     1                  w(ilist),w(iproj),w(icol),w(iwork))
c
c
        return
        end
c
c
c
c
        subroutine idzr_asvd0(m,n,a,krank,winit,u,v,s,ier,
     1                        list,proj,col,work)
c
c       routine idzr_asvd serves as a memory wrapper
c       for the present routine (please see routine idzr_asvd
c       for further documentation).
c
        implicit none
        integer m,n,krank,list(n),ier
        real*8 s(krank)
        complex*16 a(m,n),u(m,krank),v(n,krank),
     1             proj(krank,n-krank),col(m*krank),
     2             winit((2*krank+17)*n+21*m+80),
     3             work((krank+1)*(m+3*n+10)+9*krank**2)
c
c
c       ID a.
c
        call idzr_aid(m,n,a,krank,winit,list,proj)
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idz_copycols(m,n,a,krank,list,col)
c
c
c       Convert the ID to an SVD.
c
        call idz_id2svd(m,krank,col,n,list,proj,u,v,s,ier,work)
c
c
        return
        end
