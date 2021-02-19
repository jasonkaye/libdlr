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
        integer m,n,krank,list(len),ier,lproj
        real*8 a(len),p2,p3,p4,eps,col(len),
     1         errmax,errrms,proj(len),b(len)
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
c
c
c       ID a via a randomized algorithm.
c
        eps = .1d-12
        lproj = len
c
        call iddp_rid(lproj,eps,m,n,matvect,a,p2,p3,p4,
     1                krank,list,proj,ier)
c
        call prinf('ier = *',ier,1)
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
c       krank -- rank of the matrix to be constructed
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
c       routine iddp_rid computes the ID, to a specified precision,
c       of a matrix specified by a routine for applying its transpose
c       to arbitrary vectors. This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddp_rid(lproj,eps,m,n,matvect,p1,p2,p3,p4,
     1                      krank,list,proj,ier)
c
c       computes the ID of a, i.e., lists in list the indices
c       of krank columns of a such that
c
c       a(j,list(k))  =  a(j,list(k))
c
c       for all j = 1, ..., m; k = 1, ..., krank, and
c
c                        krank
c       a(j,list(k))  =  Sigma  a(j,list(l)) * proj(l,k-krank)       (*)
c                         l=1
c
c                     +  epsilon(j,k-krank)
c
c       for all j = 1, ..., m; k = krank+1, ..., n,
c
c       for some matrix epsilon dimensioned epsilon(m,n-krank)
c       such that the greatest singular value of epsilon
c       <= the greatest singular value of a * eps.
c
c       input:
c       lproj -- maximum usable length (in real*8 elements)
c                of the array proj
c       eps -- precision to which the ID is to be computed
c       m -- first dimension of a
c       n -- second dimension of a
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
c
c       output:
c       krank -- numerical rank
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd;
c               the present routine uses proj as a work array, too, so
c               proj must be at least m+1 + 2*n*(krank+1) real*8
c               elements long, where krank is the rank output
c               by the present routine
c       ier -- 0 when the routine terminates successfully;
c              -1000 when lproj is too small
c
c       _N.B._: The algorithm used by this routine is randomized.
c               proj must be at least m+1 + 2*n*(krank+1) real*8
c               elements long, where krank is the rank output
c               by the present routine.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,n,list(n),krank,lw,iwork,lwork,ira,kranki,lproj,
     1          lra,ier,k
        real*8 eps,p1,p2,p3,p4,proj(*)
        external matvect
c
c
        ier = 0
c
c
c       Allocate memory in proj.
c
        lw = 0
c
        iwork = lw+1
        lwork = m+2*n+1
        lw = lw+lwork
c
        ira = lw+1
c
c
c       Find the rank of a.
c
        lra = lproj-lwork
        call idd_findrank(lra,eps,m,n,matvect,p1,p2,p3,p4,
     1                    kranki,proj(ira),ier,proj(iwork))
        if(ier .ne. 0) return
c
c
        if(lproj .lt. lwork+2*kranki*n) then
          ier = -1000
          return
        endif
c
c
c       Transpose ra.
c
        call idd_rtransposer(n,kranki,proj(ira),proj(ira+kranki*n))
c
c
c       Move the tranposed matrix to the beginning of proj.
c
        do k = 1,kranki*n
          proj(k) = proj(ira+kranki*n+k-1)
        enddo ! k
c
c
c       ID the transposed matrix.
c
        call iddp_id(eps,kranki,n,proj,krank,list,proj(1+kranki*n))
c
c
        return
        end
c
c
c
c
        subroutine idd_findrank(lra,eps,m,n,matvect,p1,p2,p3,p4,
     1                          krank,ra,ier,w)
c
c       estimates the numerical rank krank of a matrix a to precision
c       eps, where the routine matvect applies the transpose of a
c       to an arbitrary vector. This routine applies the transpose of a
c       to krank random vectors, and returns the resulting vectors
c       as the columns of ra.
c
c       input:
c       lra -- maximum usable length (in real*8 elements) of array ra
c       eps -- precision defining the numerical rank
c       m -- first dimension of a
c       n -- second dimension of a
c       matvect -- routine which applies the transpose
c                  of the matrix whose rank is to be estimated
c                  to an arbitrary vector; this routine must have
c                  a calling sequence of the form
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
c
c       output:
c       krank -- estimate of the numerical rank of a
c       ra -- product of the transpose of a and a matrix whose entries
c             are pseudorandom realizations of i.i.d. random numbers,
c             uniformly distributed on [0,1];
c             ra must be at least 2*n*krank real*8 elements long
c       ier -- 0 when the routine terminates successfully;
c              -1000 when lra is too small
c
c       work:
c       w -- must be at least m+2*n+1 real*8 elements long
c
c       _N.B._: ra must be at least 2*n*krank real*8 elements long.
c               Also, the algorithm used by this routine is randomized.
c
        implicit none
        integer m,n,lw,krank,ix,lx,iy,ly,iscal,lscal,lra,ier
        real*8 eps,p1,p2,p3,p4,ra(n,*),w(m+2*n+1)
        external matvect
c
c
        lw = 0
c
        ix = lw+1
        lx = m
        lw = lw+lx
c
        iy = lw+1
        ly = n
        lw = lw+ly
c
        iscal = lw+1
        lscal = n+1
        lw = lw+lscal
c
c
        call idd_findrank0(lra,eps,m,n,matvect,p1,p2,p3,p4,
     1                     krank,ra,ier,w(ix),w(iy),w(iscal))
c
c
        return
        end
c
c
c
c
        subroutine idd_findrank0(lra,eps,m,n,matvect,p1,p2,p3,p4,
     1                           krank,ra,ier,x,y,scal)
c
c       routine idd_findrank serves as a memory wrapper
c       for the present routine. (Please see routine idd_findrank
c       for further documentation.)
c
        implicit none
        integer m,n,krank,ifrescal,k,lra,ier
        real*8 x(m),ra(n,2,*),p1,p2,p3,p4,scal(n+1),y(n),eps,residual,
     1         enorm
        external matvect
c
c
        ier = 0
c
c
        krank = 0
c
c
c       Loop until the relative residual is greater than eps,
c       or krank = m or krank = n.
c
 1000   continue
c
c
          if(lra .lt. n*2*(krank+1)) then
            ier = -1000
            return
          endif
c
c
c         Apply the transpose of a to a random vector.
c
          call id_srand(m,x)
          call matvect(m,x,n,ra(1,1,krank+1),p1,p2,p3,p4)
c
          do k = 1,n
            y(k) = ra(k,1,krank+1)
          enddo ! k
c
c
          if(krank .eq. 0) then
c
c           Compute the Euclidean norm of y.
c
            enorm = 0
c
            do k = 1,n
              enorm = enorm + y(k)**2
            enddo ! k
c
            enorm = sqrt(enorm)
c
          endif ! krank .eq. 0
c
c
          if(krank .gt. 0) then
c
c           Apply the previous Householder transformations to y.
c
            ifrescal = 0
c
            do k = 1,krank
              call idd_houseapp(n-k+1,ra(1,2,k),y(k),
     1                          ifrescal,scal(k),y(k))
            enddo ! k
c
          endif ! krank .gt. 0
c
c
c         Compute the Householder vector associated with y.
c
          call idd_house(n-krank,y(krank+1),
     1                   residual,ra(1,2,krank+1),scal(krank+1))
          residual = abs(residual)
c
c
          krank = krank+1
c
c
        if(residual .gt. eps*enorm
     1   .and. krank .lt. m .and. krank .lt. n)
     2   goto 1000
c
c
c       Delete the Householder vectors from the array ra.
c
        call idd_crunch(n,krank,ra)
c
c
        return
        end
c
c
c
c
        subroutine idd_crunch(n,l,a)
c
c       removes every other block of n entries from a vector.
c
c       input:
c       n -- length of each block to remove
c       l -- half of the total number of blocks
c       a -- original array
c
c       output:
c       a -- array with every other block of n entries removed
c
        implicit none
        integer j,k,n,l
        real*8 a(n,2*l)
c
c
        do j = 2,l
          do k = 1,n
c
            a(k,j) = a(k,2*j-1)
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
        subroutine idd_rtransposer(m,n,a,at)
c
c       transposes a to obtain at.
c
c       input:
c       m -- first dimension of a, and second dimension of at
c       n -- second dimension of a, and first dimension of at
c       a -- matrix to be transposed
c
c       output:
c       at -- transpose of a
c
        implicit none
        integer m,n,j,k
        real*8 a(m,n),at(n,m)
c
c
        do k = 1,n
          do j = 1,m
c
            at(k,j) = a(j,k)
c
          enddo ! j
        enddo ! k
c
c
        return
        end
