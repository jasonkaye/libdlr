c
c
c       dependencies: prini, idz_house, idz_qrpiv, idz_id, id_rand,
c                     idz_sfft, id_rtrans, idz_frm, dfft
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,list(len),n2
        real*8 errmax,errrms,eps
        complex*16 a(len),col(len),work(len),proj(len),b(len)
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
c       Initialize the array work for use in idzp_aid.
c
        call idz_frmi(m,n2,work)
c
c
c       ID a via a randomized algorithm.
c
        eps = .1d-10
        call idzp_aid(eps,m,n,a,work,krank,list,proj)
        call prinf('list = *',list,krank)
c
c
c       Collect together the columns of a indexed by list into col.
c
        call idz_copycols(m,n,a,krank,list,col)
c
c
c       Reconstruct a, obtaining b.
c
        call idz_reconid(m,krank,col,n,list,proj,b)
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
c       and left and right singular vectors taken from the DFT.
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
        real*8 r1,pi
        complex*16 a(m,n),sum,ci
c
        r1 = 1
        pi = 4*atan(r1)
        ci = (0,1)
c
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+exp(2*pi*ci*(j-r1)*(l-r1)/m)*sqrt(r1/m)
     1                 *exp(2*pi*ci*(k-r1)*(l-r1)/n)*sqrt(r1/n)
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine idzp_aid computes the ID, to a specified precision,
c       of an arbitrary matrix. This routine is randomized.
c
c       routine idz_estrank estimates the numerical rank,
c       to a specified precision, of an arbitrary matrix.
c       This routine is randomized.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idzp_aid(eps,m,n,a,work,krank,list,proj)
c
c       computes the ID of the matrix a, i.e., lists in list
c       the indices of krank columns of a such that
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
c       eps -- precision to which the ID is to be computed
c       m -- first dimension of a
c       n -- second dimension of a
c       a -- matrix to be decomposed; the present routine does not
c            alter a
c       work -- initialization array that has been constructed
c               by routine idz_frmi
c
c       output:
c       krank -- numerical rank of a to precision eps
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd;
c               proj doubles as a work array in the present routine, so
c               proj must be at least n*(2*n2+1)+n2+1 complex*16
c               elements long, where n2 is the greatest integer
c               less than or equal to m, such that n2 is
c               a positive integer power of two.
c
c       _N.B._: The algorithm used by this routine is randomized.
c               proj must be at least n*(2*n2+1)+n2+1 complex*16
c               elements long, where n2 is the greatest integer
c               less than or equal to m, such that n2 is
c               a positive integer power of two.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,n,list(n),krank,kranki,n2
        real*8 eps
        complex*16 a(m,n),proj(*),work(17*m+70)
c
c
c       Allocate memory in proj.
c
        n2 = work(2)
c
c
c       Find the rank of a.
c
        call idz_estrank(eps,m,n,a,work,kranki,proj)
c
c
        if(kranki .eq. 0) call idzp_aid0(eps,m,n,a,krank,list,proj,
     1                                   proj(m*n+1))
c
        if(kranki .ne. 0) call idzp_aid1(eps,n2,n,kranki,proj,
     1                                   krank,list,proj(n2*n+1))
c
c
        return
        end
c
c
c
c
        subroutine idzp_aid0(eps,m,n,a,krank,list,proj,rnorms)
c
c       uses routine idzp_id to ID a without modifying its entries
c       (in contrast to the usual behavior of idzp_id).
c
c       input:
c       eps -- precision of the decomposition to be constructed
c       m -- first dimension of a
c       n -- second dimension of a
c
c       output:
c       krank -- numerical rank of the ID
c       list -- indices of the columns in the ID
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns in a;
c               proj doubles as a work array in the present routine, so
c               must be at least m*n complex*16 elements long
c
c       work:
c       rnorms -- must be at least n real*8 elements long
c
c       _N.B._: proj must be at least m*n complex*16 elements long
c
        implicit none
        integer m,n,krank,list(n),j,k
        real*8 eps,rnorms(n)
        complex*16 a(m,n),proj(m,n)
c
c
c       Copy a into proj.
c
        do k = 1,n
          do j = 1,m
            proj(j,k) = a(j,k)
          enddo ! j
        enddo ! k
c
c
c       ID proj.
c
        call idzp_id(eps,m,n,proj,krank,list,rnorms)
c
c
        return
        end
c
c
c
c
        subroutine idzp_aid1(eps,n2,n,kranki,proj,krank,list,rnorms)
c
c       IDs the uppermost kranki x n block of the n2 x n matrix
c       input as proj.
c
c       input:
c       eps -- precision of the decomposition to be constructed
c       n2 -- first dimension of proj as input
c       n -- second dimension of proj as input
c       kranki -- number of rows to extract from proj
c       proj -- matrix containing the kranki x n block to be ID'd
c
c       output:
c       proj -- matrix of coefficients needed to interpolate
c               from the selected columns to the other columns
c               in the original matrix being ID'd
c       krank -- numerical rank of the ID
c       list -- indices of the columns in the ID
c
c       work:
c       rnorms -- must be at least n real*8 elements long
c
        implicit none
        integer n,n2,kranki,krank,list(n),j,k
        real*8 eps,rnorms(n)
        complex*16 proj(n2*n)
c
c
c       Move the uppermost kranki x n block of the n2 x n matrix proj
c       to the beginning of proj.
c
        do k = 1,n
          do j = 1,kranki
            proj(j+kranki*(k-1)) = proj(j+n2*(k-1))
          enddo ! j
        enddo ! k
c
c
c       ID proj.
c
        call idzp_id(eps,kranki,n,proj,krank,list,rnorms)
c
c
        return
        end
c
c
c
c
        subroutine idz_estrank(eps,m,n,a,w,krank,ra)
c
c       estimates the numerical rank krank of an m x n matrix a
c       to precision eps. This routine applies n2 random vectors
c       to a, obtaining ra, where n2 is the greatest integer
c       less than or equal to m such that n2 is a positive integer
c       power of two. krank is typically about 8 higher than
c       the actual numerical rank.
c
c       input:
c       eps -- precision defining the numerical rank
c       m -- first dimension of a
c       n -- second dimension of a
c       a -- matrix whose rank is to be estimated
c       w -- initialization array that has been constructed
c            by routine idz_frmi
c
c       output:
c       krank -- estimate of the numerical rank of a;
c                this routine returns krank = 0 when the actual
c                numerical rank is nearly full (that is,
c                greater than n - 8 or n2 - 8)
c       ra -- product of an n2 x m random matrix and the m x n matrix
c             a, where n2 is the greatest integer less than or equal
c             to m such that n2 is a positive integer power of two;
c             ra doubles as a work array in the present routine, and so
c             must be at least n*n2+(n+1)*(n2+1) complex*16 elements
c             long
c
c       _N.B._: ra must be at least n*n2+(n2+1)*(n+1) complex*16
c               elements long for use in the present routine
c               (here, n2 is the greatest integer less than or equal
c               to m, such that n2 is a positive integer power of two).
c               This routine returns krank = 0 when the actual
c               numerical rank is nearly full.
c
        implicit none
        integer m,n,krank,n2,irat,lrat,iscal,lscal,ira,lra,lra2
        real*8 eps
        complex*16 a(m,n),ra(*),w(17*m+70)
c
c
c       Extract from the array w initialized by routine idz_frmi
c       the greatest integer less than or equal to m that is
c       a positive integer power of two.
c
        n2 = w(2)
c
c
c       Allocate memory in ra.
c
        lra = 0
c
        ira = lra+1
        lra2 = n2*n
        lra = lra+lra2
c
        irat = lra+1
        lrat = n*(n2+1)
        lra = lra+lrat
c
        iscal = lra+1
        lscal = n2+1
        lra = lra+lscal
c
        call idz_estrank0(eps,m,n,a,w,n2,krank,ra(ira),ra(irat),
     1                    ra(iscal))
c
c
        return
        end
c
c
c
c
        subroutine idz_estrank0(eps,m,n,a,w,n2,krank,ra,rat,scal)
c
c       routine idz_estrank serves as a memory wrapper
c       for the present routine. (Please see routine idz_estrank
c       for further documentation.)
c
        implicit none
        integer m,n,n2,krank,ifrescal,k,nulls,j
        real*8 eps,scal(n2+1),ss,ssmax
        complex*16 a(m,n),ra(n2,n),residual,w(17*m+70),rat(n,n2+1)
c
c
c       Apply the random matrix to every column of a, obtaining ra.
c
        do k = 1,n
          call idz_frm(m,n2,w,a(1,k),ra(1,k))
        enddo ! k
c
c
c       Compute the sum of squares of the entries in each column of ra
c       and the maximum of all such sums.
c
        ssmax = 0
c
        do k = 1,n
c
          ss = 0
          do j = 1,m
            ss = ss+a(j,k)*conjg(a(j,k))
          enddo ! j
c
          if(ss .gt. ssmax) ssmax = ss
c
        enddo ! k
c
c
c       Transpose ra to obtain rat.
c
        call idz_transposer(n2,n,ra,rat)
c
c
        krank = 0
        nulls = 0
c
c
c       Loop until nulls = 7, krank+nulls = n2, or krank+nulls = n.
c
 1000   continue
c
c
          if(krank .gt. 0) then
c
c           Apply the previous Householder transformations
c           to rat(:,krank+1).
c
            ifrescal = 0
c
            do k = 1,krank
              call idz_houseapp(n-k+1,rat(1,k),rat(k,krank+1),
     1                          ifrescal,scal(k),rat(k,krank+1))
            enddo ! k
c
          endif ! krank .gt. 0
c
c
c         Compute the Householder vector associated
c         with rat(krank+1:*,krank+1).
c
          call idz_house(n-krank,rat(krank+1,krank+1),
     1                   residual,rat(1,krank+1),scal(krank+1))
c
c
          krank = krank+1
          if(abs(residual) .le. eps*sqrt(ssmax)) nulls = nulls+1
c
c
        if(nulls .lt. 7 .and. krank+nulls .lt. n2
     1   .and. krank+nulls .lt. n)
     2   goto 1000
c
c
        if(nulls .lt. 7) krank = 0
c
c
        return
        end
c
c
c
c
        subroutine idz_transposer(m,n,a,at)
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
        complex*16 a(m,n),at(n,m)
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
