c
c
c       dependencies: prini, idd_house, idd_qrpiv, idd_id, id_rand,
c                     iddr_rid, idd_id2svd, lapack.a, blas.a
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,ier
        real*8 a(len),b(len),dummy,work(len),
     1         errmax,errrms,u(len),v(len),s(100 000)
        external matvect,matvec
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
        call iddr_rsvd(m,n,matvect,a,dummy,dummy,dummy,
     1                 matvec,a,dummy,dummy,dummy,krank,
     2                 u,v,s,ier,work)
c
c
c       Construct b = u diag(s) v^T.
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
        subroutine matvect(m,x,n,y,a,p2,p3,p4)
c
c       applies the transpose of a to x, obtaining y.
c
c       input:
c       m -- first dimension of a, and length of x
c       x -- vector to which a^T is to be applied
c       n -- second dimension of a, and length of y
c       a -- matrix whose transpose is to be applied to x
c            in order to create y
c       p2 -- dummy input
c       p3 -- dummy input
c       p4 -- dummy input
c
c       output:
c       y -- product of a^T and x
c
        implicit none
        integer m,n,j,k
        real*8 a(m,n),p2,p3,p4,x(m),y(n),sum
c
c
        do k = 1,n
c
          sum = 0
c
          do j = 1,m
            sum = sum+a(j,k)*x(j)
          enddo ! j
c
          y(k) = sum
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
        subroutine matvec(n,x,m,y,a,p2,p3,p4)
c
c       applies a to x, obtaining y.
c
c       input:
c       m -- first dimension of a, and length of x
c       x -- vector to which a is to be applied
c       n -- second dimension of a, and length of y
c       a -- matrix to be applied to x in order to create y
c       p2 -- dummy input
c       p3 -- dummy input
c       p4 -- dummy input
c
c       output:
c       y -- product of a and x
c
        implicit none
        integer m,n,j,k
        real*8 a(m,n),p2,p3,p4,x(n),y(m),sum
c
c
        do j = 1,m
c
          sum = 0
c
          do k = 1,n
            sum = sum+a(j,k)*x(k)
          enddo ! k
c
          y(j) = sum
c
        enddo ! j
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine iddr_rsvd computes the SVD, to a specified rank,
c       of a matrix specified by routines for applying the matrix
c       and its transpose to arbitrary vectors.
c       This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddr_rsvd(m,n,matvect,p1t,p2t,p3t,p4t,
     1                       matvec,p1,p2,p3,p4,krank,u,v,s,ier,w)
c
c       constructs a rank-krank SVD  u diag(s) v^T  approximating a,
c       where matvect is a routine which applies a^T
c       to an arbitrary vector, and matvec is a routine
c       which applies a to an arbitrary vector;
c       u is an m x krank matrix whose columns are orthonormal,
c       v is an n x krank matrix whose columns are orthonormal,
c       and diag(s) is a diagonal krank x krank matrix whose entries
c       are all nonnegative. This routine uses a randomized algorithm.
c
c       input:
c       m -- number of rows in a
c       n -- number of columns in a 
c       matvect -- routine which applies the transpose
c                  of the matrix to be SVD'd
c                  to an arbitrary vector; this routine must have
c                  a calling sequence of the form
c
c                  matvect(m,x,n,y,p1t,p2t,p3t,p4t),
c
c                  where m is the length of x,
c                  x is the vector to which the transpose
c                  of the matrix is to be applied,
c                  n is the length of y,
c                  y is the product of the transposed matrix and x,
c                  and p1t, p2t, p3t, and p4t are user-specified
c                  parameters
c       p1t -- parameter to be passed to routine matvect
c       p2t -- parameter to be passed to routine matvect
c       p3t -- parameter to be passed to routine matvect
c       p4t -- parameter to be passed to routine matvect
c       matvec -- routine which applies the matrix to be SVD'd
c                 to an arbitrary vector; this routine must have
c                 a calling sequence of the form
c
c                 matvec(n,x,m,y,p1,p2,p3,p4),
c
c                 where n is the length of x,
c                 x is the vector to which the matrix is to be applied,
c                 m is the length of y,
c                 y is the product of the matrix and x,
c                 and p1, p2, p3, and p4 are user-specified parameters
c       p1 -- parameter to be passed to routine matvec
c       p2 -- parameter to be passed to routine matvec
c       p3 -- parameter to be passed to routine matvec
c       p4 -- parameter to be passed to routine matvec
c       krank -- rank of the SVD being constructed
c
c       output:
c       u -- matrix of orthonormal left singular vectors of a
c       v -- matrix of orthonormal right singular vectors of a
c       s -- array of singular values of a
c       ier -- 0 when the routine terminates successfully;
c              nonzero otherwise
c
c       work:
c       w -- must be at least (krank+1)*(2*m+4*n)+25*krank**2
c            real*8 elements long
c
c       _N.B._: The algorithm used by this routine is randomized.
c
        implicit none
        integer m,n,krank,lw,ilist,llist,iproj,lproj,icol,lcol,
     1          iwork,lwork,ier
        real*8 p1t,p2t,p3t,p4t,p1,p2,p3,p4,u(m,krank),v(n,krank),
     1         s(krank),w((krank+1)*(2*m+4*n)+25*krank**2)
        external matvect,matvec
c
c
c       Allocate memory in w.
c
        lw = 0
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
        lwork = (krank+1)*(m+3*n)+26*krank**2
        lw = lw+lwork
c
c
        call iddr_rsvd0(m,n,matvect,p1t,p2t,p3t,p4t,
     1                  matvec,p1,p2,p3,p4,krank,u,v,s,ier,
     2                  w(ilist),w(iproj),w(icol),w(iwork))
c
c
        return
        end
c
c
c
c
        subroutine iddr_rsvd0(m,n,matvect,p1t,p2t,p3t,p4t,
     1                        matvec,p1,p2,p3,p4,krank,u,v,s,ier,
     2                        list,proj,col,work)
c
c       routine iddr_rsvd serves as a memory wrapper
c       for the present routine (please see routine iddr_rsvd
c       for further documentation).
c
        implicit none
        integer m,n,krank,list(n),ier,k
        real*8 p1t,p2t,p3t,p4t,p1,p2,p3,p4,u(m,krank),v(n,krank),
     1         s(krank),proj(krank*(n-krank)),col(m*krank),
     2         work((krank+1)*(m+3*n)+26*krank**2)
        external matvect,matvec
c
c
c       ID a.
c
        call iddr_rid(m,n,matvect,p1t,p2t,p3t,p4t,krank,list,work)
c
c
c       Retrieve proj from work.
c
        do k = 1,krank*(n-krank)
          proj(k) = work(k)
        enddo ! k
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idd_getcols(m,n,matvec,p1,p2,p3,p4,krank,list,col,work)
c
c
c       Convert the ID to an SVD.
c
        call idd_id2svd(m,krank,col,n,list,proj,u,v,s,ier,work)
c
c
        return
        end
