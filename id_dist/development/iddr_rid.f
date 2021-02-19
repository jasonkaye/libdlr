c
c
c       dependencies: prini, idd_house, idd_qrpiv, idd_id, id_rand
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,list(len)
        real*8 a(len),p2,p3,p4,b(len),proj(len),
     1         col(len),errmax,errrms
        external matvect
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
c       ID a.
c
        call iddr_rid(m,n,matvect,a,p2,p3,p4,krank,list,proj)
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine iddr_rid computes the ID, to a specified rank,
c       of a matrix specified by a routine for applying its transpose
c       to arbitrary vectors. This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddr_rid(m,n,matvect,p1,p2,p3,p4,krank,list,proj)
c
c       computes the ID of a matrix "a" specified by
c       the routine matvect -- matvect must apply the transpose
c       of the matrix being ID'd to an arbitrary vector --
c       i.e., the present routine lists in list the indices
c       of krank columns of a such that 
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
c       m -- number of rows in the matrix to be ID'd
c       n -- number of columns in the matrix to be ID'd
c       matvect -- routine which applies the transpose
c                  of the matrix to be ID'd to an arbitrary vector;
c                  this routine must have a calling sequence
c                  of the form
c
c                  matvect(m,x,n,y,p1,p2,p3,p4),
c
c                  where m is the length of x,
c                  x is the vector to which the transpose
c                  of the matrix is to be applied,
c                  n is the length of y,
c                  y is the product of the transposed matrix and x,
c                  and p1, p2, p3, and p4 are user-specified parameters
c       p1 -- parameter to be passed to routine matvect
c       p2 -- parameter to be passed to routine matvect
c       p3 -- parameter to be passed to routine matvect
c       p4 -- parameter to be passed to routine matvect
c       krank -- rank of the ID to be constructed
c
c       output:
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd;
c               proj doubles as a work array in the present routine, so
c               proj must be at least m+(krank+3)*n real*8 elements
c               long
c
c       _N.B._: The algorithm used by this routine is randomized.
c               proj must be at least m+(krank+3)*n real*8 elements
c               long.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,n,krank,list(n),lw,ix,lx,iy,ly,ir,lr
        real*8 p1,p2,p3,p4,proj(m+(krank+3)*n)
        external matvect
c
c
c       Allocate memory in w.
c
        lw = 0
c
        ir = lw+1
        lr = (krank+2)*n
        lw = lw+lr
c
        ix = lw+1
        lx = m
        lw = lw+lx
c
        iy = lw+1
        ly = n
        lw = lw+ly
c
c
        call iddr_ridall0(m,n,matvect,p1,p2,p3,p4,krank,
     1                    list,proj(ir),proj(ix),proj(iy))
c
c
        return
        end
c
c
c
c
        subroutine iddr_ridall0(m,n,matvect,p1,p2,p3,p4,krank,
     1                          list,r,x,y)
c
c       routine iddr_ridall serves as a memory wrapper
c       for the present routine
c       (see iddr_ridall for further documentation).
c
        implicit none
        integer j,k,l,m,n,krank,list(n)
        real*8 x(m),y(n),p1,p2,p3,p4,r(krank+2,n)
        external matvect
c
c
c       Set the number of random test vectors to 2 more than the rank.
c
        l = krank+2
c
c       Apply the transpose of the original matrix to l random vectors.
c
        do j = 1,l
c
c         Generate a random vector.
c
          call id_srand(m,x)
c
c         Apply the transpose of the matrix to x, obtaining y.
c
          call matvect(m,x,n,y,p1,p2,p3,p4)
c
c         Copy y into row j of r.
c
          do k = 1,n
            r(j,k) = y(k)
          enddo ! k
c
        enddo ! j
c
c
c       ID r.
c
        call iddr_id(l,n,r,krank,list,y)
c
c
        return
        end
