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
c
c
c       routine iddp_qrpiv computes the pivoted QR decomposition
c       of a matrix via Householder transformations,
c       stopping at a specified precision of the decomposition.
c
c       routine iddr_qrpiv computes the pivoted QR decomposition
c       of a matrix via Householder transformations,
c       stopping at a specified rank of the decomposition.
c
c       routine idd_qmatvec applies to a single vector
c       the Q matrix (or its transpose) in the QR decomposition
c       of a matrix, as described by the output of iddp_qrpiv
c       or iddr_qrpiv. If you're concerned about efficiency
c       and want to apply Q (or its transpose) to multiple vectors,
c       use idd_qmatmat instead.
c
c       routine idd_qmatmat applies
c       to multiple vectors collected together
c       as a matrix the Q matrix (or its transpose)
c       in the QR decomposition of a matrix, as described
c       by the output of iddp_qrpiv or iddr_qrpiv. If you don't want
c       to provide a work array and want to apply Q (or its transpose)
c       to a single vector, use idd_qmatvec instead.
c
c       routine idd_qinqr reconstructs the Q matrix
c       in a QR decomposition from the data generated
c       by iddp_qrpiv or iddr_qrpiv.
c
c       routine idd_permmult multiplies together a bunch
c       of permutations.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c

        subroutine idd_permmult(m,ind,n,indprod)
c
c       multiplies together the series of permutations in ind.
c
c       input:
c       m -- length of ind
c       ind(k) -- number of the slot with which to swap
c                 the k^th slot
c       n -- length of indprod and indprodinv
c
c       output:
c       indprod -- product of the permutations in ind,
c                  with the permutation swapping 1 and ind(1)
c                  taken leftmost in the product,
c                  that swapping 2 and ind(2) taken next leftmost,
c                  ..., that swapping krank and ind(krank)
c                  taken rightmost; indprod(k) is the number
c                  of the slot with which to swap the k^th slot
c                  in the product permutation
c
        implicit none
        integer m,n,ind(m),indprod(n),k,iswap
c
c
        do k = 1,n
          indprod(k) = k
        enddo ! k 
c
        do k = m,1,-1
c
c         Swap indprod(k) and indprod(ind(k)).
c
          iswap = indprod(k)
          indprod(k) = indprod(ind(k))
          indprod(ind(k)) = iswap
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
        subroutine idd_qinqr(m,n,a,krank,q)
c
c       constructs the matrix q from iddp_qrpiv or iddr_qrpiv
c       (see the routine iddp_qrpiv or iddr_qrpiv
c       for more information).
c
c       input:
c       m -- first dimension of a; also, right now, q is m x m
c       n -- second dimension of a
c       a -- matrix output by iddp_qrpiv or iddr_qrpiv
c            (and denoted the same there)
c       krank -- numerical rank output by iddp_qrpiv or iddr_qrpiv
c                (and denoted the same there)
c
c       output:
c       q -- orthogonal matrix implicitly specified by the data in a
c            from iddp_qrpiv or iddr_qrpiv
c
c       Note:
c       Right now, this routine simply multiplies
c       one after another the krank Householder matrices
c       in the full QR decomposition of a,
c       in order to obtain the complete m x m Q factor in the QR.
c       This routine should instead use the following 
c       (more elaborate but more efficient) scheme
c       to construct a q dimensioned q(krank,m); this scheme
c       was introduced by Robert Schreiber and Charles Van Loan
c       in "A Storage-Efficient _WY_ Representation
c       for Products of Householder Transformations,"
c       _SIAM Journal on Scientific and Statistical Computing_,
c       Vol. 10, No. 1, pp. 53-57, January, 1989:
c
c       Theorem 1. Suppose that Q = _1_ + YTY^T is
c       an m x m orthogonal real matrix,
c       where Y is an m x k real matrix
c       and T is a k x k upper triangular real matrix.
c       Suppose also that P = _1_ - 2 v v^T is
c       a real Householder matrix and Q_+ = QP,
c       where v is an m x 1 real vector,
c       normalized so that v^T v = 1.
c       Then, Q_+ = _1_ + Y_+ T_+ Y_+^T,
c       where Y_+ = (Y v) is the m x (k+1) matrix
c       formed by adjoining v to the right of Y,
c                 ( T   z )
c       and T_+ = (       ) is
c                 ( 0  -2 )
c       the (k+1) x (k+1) upper triangular matrix
c       formed by adjoining z to the right of T
c       and the vector (0 ... 0 -2) with k zeroes below (T z),
c       where z = -2 T Y^T v.
c
c       Now, suppose that A is a (rank-deficient) matrix
c       whose complete QR decomposition has
c       the blockwise partioned form
c           ( Q_11 Q_12 ) ( R_11 R_12 )   ( Q_11 )
c       A = (           ) (           ) = (      ) (R_11 R_12).
c           ( Q_21 Q_22 ) (  0    0   )   ( Q_21 )
c       Then, the only blocks of the orthogonal factor
c       in the above QR decomposition of A that matter are
c                                                        ( Q_11 )
c       Q_11 and Q_21, _i.e._, only the block of columns (      )
c                                                        ( Q_21 )
c       interests us.
c       Suppose in addition that Q_11 is a k x k matrix,
c       Q_21 is an (m-k) x k matrix, and that
c       ( Q_11 Q_12 )
c       (           ) = _1_ + YTY^T, as in Theorem 1 above.
c       ( Q_21 Q_22 )
c       Then, Q_11 = _1_ + Y_1 T Y_1^T
c       and Q_21 = Y_2 T Y_1^T,
c       where Y_1 is the k x k matrix and Y_2 is the (m-k) x k matrix
c                   ( Y_1 )
c       so that Y = (     ).
c                   ( Y_2 )
c
c       So, you can calculate T and Y via the above recursions,
c       and then use these to compute the desired Q_11 and Q_21.
c
c
        implicit none
        integer m,n,krank,j,k,mm,ifrescal
        real*8 a(m,n),q(m,m),scal
c
c
c       Zero all of the entries of q.
c
        do k = 1,m
          do j = 1,m
            q(j,k) = 0
          enddo ! j
        enddo ! k
c
c
c       Place 1's along the diagonal of q.
c
        do k = 1,m
          q(k,k) = 1
        enddo ! k
c
c
c       Apply the krank Householder transformations stored in a.
c
        do k = krank,1,-1
          do j = k,m
            mm = m-k+1
            ifrescal = 1
            if(k .lt. m)
     1       call idd_houseapp(mm,a(k+1,k),q(k,j),ifrescal,scal,q(k,j))
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
        subroutine idd_qmatvec(iftranspose,m,n,a,krank,v)
c
c       applies to a single vector the Q matrix (or its transpose)
c       which the routine iddp_qrpiv or iddr_qrpiv has stored
c       in a triangle of the matrix it produces (stored, incidentally,
c       as data for applying a bunch of Householder reflections).
c       Use the routine qmatmat to apply the Q matrix
c       (or its transpose)
c       to a bunch of vectors collected together as a matrix,
c       if you're concerned about efficiency.
c
c       input:
c       iftranspose -- set to 0 for applying Q;
c                      set to 1 for applying the transpose of Q
c       m -- first dimension of a and length of v
c       n -- second dimension of a
c       a -- data describing the qr decomposition of a matrix,
c            as produced by iddp_qrpiv or iddr_qrpiv
c       krank -- numerical rank
c       v -- vector to which Q (or its transpose) is to be applied
c
c       output:
c       v -- vector to which Q (or its transpose) has been applied
c
        implicit none
        save
        integer m,n,krank,k,ifrescal,mm,iftranspose
        real*8 a(m,n),v(m),scal
c
c
        ifrescal = 1
c
c
        if(iftranspose .eq. 0) then
c
          do k = krank,1,-1
            mm = m-k+1
            if(k .lt. m)
     1       call idd_houseapp(mm,a(k+1,k),v(k),ifrescal,scal,v(k))
          enddo ! k
c
        endif
c
c
        if(iftranspose .eq. 1) then
c
          do k = 1,krank
            mm = m-k+1
            if(k .lt. m)
     1       call idd_houseapp(mm,a(k+1,k),v(k),ifrescal,scal,v(k))
          enddo ! k
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
        subroutine idd_qmatmat(iftranspose,m,n,a,krank,l,b,work)
c
c       applies to a bunch of vectors collected together as a matrix
c       the Q matrix (or its transpose) which the routine iddp_qrpiv or
c       iddr_qrpiv has stored in a triangle of the matrix it produces
c       (stored, incidentally, as data for applying a bunch
c       of Householder reflections).
c       Use the routine qmatvec to apply the Q matrix
c       (or its transpose)
c       to a single vector, if you'd rather not provide a work array.
c
c       input:
c       iftranspose -- set to 0 for applying Q;
c                      set to 1 for applying the transpose of Q
c       m -- first dimension of both a and b
c       n -- second dimension of a
c       a -- data describing the qr decomposition of a matrix,
c            as produced by iddp_qrpiv or iddr_qrpiv
c       krank -- numerical rank
c       l -- second dimension of b
c       b -- matrix to which Q (or its transpose) is to be applied
c
c       output:
c       b -- matrix to which Q (or its transpose) has been applied
c
c       work:
c       work -- must be at least krank real*8 elements long
c
        implicit none
        save
        integer l,m,n,krank,j,k,ifrescal,mm,iftranspose
        real*8 a(m,n),b(m,l),work(krank)
c
c
        if(iftranspose .eq. 0) then
c
c
c         Handle the first iteration, j = 1,
c         calculating all scals (ifrescal = 1).
c
          ifrescal = 1
c
          j = 1
c
          do k = krank,1,-1
            if(k .lt. m) then
              mm = m-k+1
              call idd_houseapp(mm,a(k+1,k),b(k,j),ifrescal,
     1                          work(k),b(k,j))
            endif
          enddo ! k
c
c
          if(l .gt. 1) then
c
c           Handle the other iterations, j > 1,
c           using the scals just computed (ifrescal = 0).
c
            ifrescal = 0
c
            do j = 2,l
c
              do k = krank,1,-1
                if(k .lt. m) then
                  mm = m-k+1
                  call idd_houseapp(mm,a(k+1,k),b(k,j),ifrescal,
     1                              work(k),b(k,j))
                endif
              enddo ! k
c
            enddo ! j
c
          endif ! j .gt. 1
c
c
        endif ! iftranspose .eq. 0
c
c
        if(iftranspose .eq. 1) then
c
c
c         Handle the first iteration, j = 1,
c         calculating all scals (ifrescal = 1).
c
          ifrescal = 1
c
          j = 1
c
          do k = 1,krank
            if(k .lt. m) then
              mm = m-k+1
              call idd_houseapp(mm,a(k+1,k),b(k,j),ifrescal,
     1                          work(k),b(k,j))
            endif
          enddo ! k
c
c
          if(l .gt. 1) then
c
c           Handle the other iterations, j > 1,
c           using the scals just computed (ifrescal = 0).
c
            ifrescal = 0
c
            do j = 2,l
c
              do k = 1,krank
                if(k .lt. m) then
                  mm = m-k+1
                  call idd_houseapp(mm,a(k+1,k),b(k,j),ifrescal,
     1                              work(k),b(k,j))
                endif
              enddo ! k
c
            enddo ! j
c
          endif ! j .gt. 1
c
c
        endif ! iftranspose .eq. 1
c
c
        return
        end
c
c
c
c
        subroutine iddp_qrpiv(eps,m,n,a,krank,ind,ss)
c
c       computes the pivoted QR decomposition
c       of the matrix input into a, using Householder transformations,
c       _i.e._, transforms the matrix a from its input value in
c       to the matrix out with entry
c
c                               m
c       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k),
c                              l=1
c
c       for all j = 1, ..., krank, and k = 1, ..., n,
c
c       where in = the a from before the routine runs,
c       out = the a from after the routine runs,
c       out(j,k) = 0 when j > k (so that out is triangular),
c       q(1:m,1), ..., q(1:m,krank) are orthonormal,
c       indprod is the product of the permutations given by ind,
c       (as computable via the routine permmult,
c       with the permutation swapping 1 and ind(1) taken leftmost
c       in the product, that swapping 2 and ind(2) taken next leftmost,
c       ..., that swapping krank and ind(krank) taken rightmost),
c       and with the matrix out satisfying
c
c                   krank
c       in(j,k)  =  Sigma  q(j,l) * out(l,indprod(k))  +  epsilon(j,k),
c                    l=1
c
c       for all j = 1, ..., m, and k = 1, ..., n,
c
c       for some matrix epsilon such that
c       the root-sum-square of the entries of epsilon
c       <= the root-sum-square of the entries of in * eps.
c       Well, technically, this routine outputs the Householder vectors
c       (or, rather, their second through last entries)
c       in the part of a that is supposed to get zeroed, that is,
c       in a(j,k) with m >= j > k >= 1.
c
c       input:
c       eps -- relative precision of the resulting QR decomposition
c       m -- first dimension of a and q
c       n -- second dimension of a
c       a -- matrix whose QR decomposition gets computed
c
c       output:
c       a -- triangular (R) factor in the QR decompositon
c            of the matrix input into the same storage locations, 
c            with the Householder vectors stored in the part of a
c            that would otherwise consist entirely of zeroes, that is,
c            in a(j,k) with m >= j > k >= 1
c       krank -- numerical rank
c       ind(k) -- index of the k^th pivot vector;
c                 the following code segment will correctly rearrange
c                 the product b of q and the upper triangle of out
c                 so that b matches the input matrix in
c                 to relative precision eps:
c
c                 copy the non-rearranged product of q and out into b
c                 set k to krank
c                 [start of loop]
c                   swap b(1:m,k) and b(1:m,ind(k))
c                   decrement k by 1
c                 if k > 0, then go to [start of loop]
c
c       work:
c       ss -- must be at least n real*8 words long
c
c       _N.B._: This routine outputs the Householder vectors
c       (or, rather, their second through last entries)
c       in the part of a that is supposed to get zeroed, that is,
c       in a(j,k) with m >= j > k >= 1.
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        integer n,m,ind(n),krank,k,j,kpiv,mm,nupdate,ifrescal
        real*8 a(m,n),ss(n),eps,feps,ssmax,scal,ssmaxin,rswap
c
c
        feps = .1d-16
c
c
c       Compute the sum of squares of the entries in each column of a,
c       the maximum of all such sums, and find the first pivot
c       (column with the greatest such sum).
c
        ssmax = 0
        kpiv = 1
c
        do k = 1,n
c
          ss(k) = 0
          do j = 1,m
            ss(k) = ss(k)+a(j,k)**2
          enddo ! j
c
          if(ss(k) .gt. ssmax) then
            ssmax = ss(k)
            kpiv = k
          endif
c
        enddo ! k
c
        ssmaxin = ssmax
c
        nupdate = 0
c
c
c       While ssmax > eps**2*ssmaxin, krank < m, and krank < n,
c       do the following block of code,
c       which ends at the statement labeled 2000.
c
        krank = 0
 1000   continue
c
        if(ssmax .le. eps**2*ssmaxin
     1   .or. krank .ge. m .or. krank .ge. n) goto 2000
        krank = krank+1
c
c
          mm = m-krank+1
c
c
c         Perform the pivoting.
c
          ind(krank) = kpiv
c
c         Swap a(1:m,krank) and a(1:m,kpiv).
c
          do j = 1,m
            rswap = a(j,krank)
            a(j,krank) = a(j,kpiv)
            a(j,kpiv) = rswap
          enddo ! j
c
c         Swap ss(krank) and ss(kpiv).
c
          rswap = ss(krank)
          ss(krank) = ss(kpiv)
          ss(kpiv) = rswap
c
c
          if(krank .lt. m) then
c
c
c           Compute the data for the Householder transformation
c           which will zero a(krank+1,krank), ..., a(m,krank)
c           when applied to a, replacing a(krank,krank)
c           with the first entry of the result of the application
c           of the Householder matrix to a(krank:m,krank),
c           and storing entries 2 to mm of the Householder vector
c           in a(krank+1,krank), ..., a(m,krank)
c           (which otherwise would get zeroed upon application
c           of the Householder transformation).
c
            call idd_house(mm,a(krank,krank),a(krank,krank),
     1                     a(krank+1,krank),scal)
            ifrescal = 0
c
c
c           Apply the Householder transformation
c           to the lower right submatrix of a
c           with upper leftmost entry at position (krank,krank+1).
c
            if(krank .lt. n) then
              do k = krank+1,n
                call idd_houseapp(mm,a(krank+1,krank),a(krank,k),
     1                            ifrescal,scal,a(krank,k))
              enddo ! k
            endif
c
c
c           Update the sums-of-squares array ss.
c
            do k = krank,n
              ss(k) = ss(k)-a(krank,k)**2
            enddo ! k
c
c
c           Find the pivot (column with the greatest sum of squares
c           of its entries).
c
            ssmax = 0
            kpiv = krank+1
c
            if(krank .lt. n) then
c
              do k = krank+1,n
c
                if(ss(k) .gt. ssmax) then
                  ssmax = ss(k)
                  kpiv = k
                endif
c
              enddo ! k
c
            endif ! krank .lt. n
c
c
c           Recompute the sums-of-squares and the pivot
c           when ssmax first falls below
c           sqrt((1000*feps)^2) * ssmaxin
c           and when ssmax first falls below
c           ((1000*feps)^2) * ssmaxin.
c
            if(
     1       (ssmax .lt. sqrt((1000*feps)**2) * ssmaxin
     2        .and. nupdate .eq. 0) .or.
     3       (ssmax .lt. ((1000*feps)**2) * ssmaxin
     4        .and. nupdate .eq. 1)
     5      ) then
c
              nupdate = nupdate+1
c
              ssmax = 0
              kpiv = krank+1
c
              if(krank .lt. n) then
c
                do k = krank+1,n
c
                  ss(k) = 0
                  do j = krank+1,m
                    ss(k) = ss(k)+a(j,k)**2
                  enddo ! j
c
                  if(ss(k) .gt. ssmax) then
                    ssmax = ss(k)
                    kpiv = k
                  endif
c
                enddo ! k
c
              endif ! krank .lt. n
c
            endif
c
c
          endif ! krank .lt. m
c
c
        goto 1000
 2000   continue
c
c
        return
        end
c
c
c
c
        subroutine iddr_qrpiv(m,n,a,krank,ind,ss)
c
c       computes the pivoted QR decomposition
c       of the matrix input into a, using Householder transformations,
c       _i.e._, transforms the matrix a from its input value in
c       to the matrix out with entry
c
c                               m
c       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k),
c                              l=1
c
c       for all j = 1, ..., krank, and k = 1, ..., n,
c
c       where in = the a from before the routine runs,
c       out = the a from after the routine runs,
c       out(j,k) = 0 when j > k (so that out is triangular),
c       q(1:m,1), ..., q(1:m,krank) are orthonormal,
c       indprod is the product of the permutations given by ind,
c       (as computable via the routine permmult,
c       with the permutation swapping 1 and ind(1) taken leftmost
c       in the product, that swapping 2 and ind(2) taken next leftmost,
c       ..., that swapping krank and ind(krank) taken rightmost),
c       and with the matrix out satisfying
c
c                  min(krank,m,n)
c       in(j,k)  =     Sigma      q(j,l) * out(l,indprod(k))
c                       l=1
c
c                +  epsilon(j,k),
c
c       for all j = 1, ..., m, and k = 1, ..., n,
c
c       for some matrix epsilon whose norm is (hopefully) minimized
c       by the pivoting procedure.
c       Well, technically, this routine outputs the Householder vectors
c       (or, rather, their second through last entries)
c       in the part of a that is supposed to get zeroed, that is,
c       in a(j,k) with m >= j > k >= 1.
c
c       input:
c       m -- first dimension of a and q
c       n -- second dimension of a
c       a -- matrix whose QR decomposition gets computed
c       krank -- desired rank of the output matrix
c                (please note that if krank > m or krank > n,
c                then the rank of the output matrix will be
c                less than krank)
c
c       output:
c       a -- triangular (R) factor in the QR decompositon
c            of the matrix input into the same storage locations, 
c            with the Householder vectors stored in the part of a
c            that would otherwise consist entirely of zeroes, that is,
c            in a(j,k) with m >= j > k >= 1
c       ind(k) -- index of the k^th pivot vector;
c                 the following code segment will correctly rearrange
c                 the product b of q and the upper triangle of out
c                 so that b best matches the input matrix in:
c
c                 copy the non-rearranged product of q and out into b
c                 set k to krank
c                 [start of loop]
c                   swap b(1:m,k) and b(1:m,ind(k))
c                   decrement k by 1
c                 if k > 0, then go to [start of loop]
c
c       work:
c       ss -- must be at least n real*8 words long
c
c       _N.B._: This routine outputs the Householder vectors
c       (or, rather, their second through last entries)
c       in the part of a that is supposed to get zeroed, that is,
c       in a(j,k) with m >= j > k >= 1.
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        integer n,m,ind(n),krank,k,j,kpiv,mm,nupdate,ifrescal,
     1          loops,loop
        real*8 a(m,n),ss(n),ssmax,scal,ssmaxin,rswap,feps
c
c
        feps = .1d-16
c
c
c       Compute the sum of squares of the entries in each column of a,
c       the maximum of all such sums, and find the first pivot
c       (column with the greatest such sum).
c
        ssmax = 0
        kpiv = 1
c
        do k = 1,n
c
          ss(k) = 0
          do j = 1,m
            ss(k) = ss(k)+a(j,k)**2
          enddo ! j
c
          if(ss(k) .gt. ssmax) then
            ssmax = ss(k)
            kpiv = k
          endif
c
        enddo ! k
c
        ssmaxin = ssmax
c
        nupdate = 0
c
c
c       Set loops = min(krank,m,n).
c
        loops = krank
        if(m .lt. loops) loops = m
        if(n .lt. loops) loops = n
c
        do loop = 1,loops
c
c
          mm = m-loop+1
c
c
c         Perform the pivoting.
c
          ind(loop) = kpiv
c
c         Swap a(1:m,loop) and a(1:m,kpiv).
c
          do j = 1,m
            rswap = a(j,loop)
            a(j,loop) = a(j,kpiv)
            a(j,kpiv) = rswap
          enddo ! j
c
c         Swap ss(loop) and ss(kpiv).
c
          rswap = ss(loop)
          ss(loop) = ss(kpiv)
          ss(kpiv) = rswap
c
c
          if(loop .lt. m) then
c
c
c           Compute the data for the Householder transformation
c           which will zero a(loop+1,loop), ..., a(m,loop)
c           when applied to a, replacing a(loop,loop)
c           with the first entry of the result of the application
c           of the Householder matrix to a(loop:m,loop),
c           and storing entries 2 to mm of the Householder vector
c           in a(loop+1,loop), ..., a(m,loop)
c           (which otherwise would get zeroed upon application
c           of the Householder transformation).
c
            call idd_house(mm,a(loop,loop),a(loop,loop),
     1                     a(loop+1,loop),scal)
            ifrescal = 0
c
c
c           Apply the Householder transformation
c           to the lower right submatrix of a
c           with upper leftmost entry at position (loop,loop+1).
c
            if(loop .lt. n) then
              do k = loop+1,n
                call idd_houseapp(mm,a(loop+1,loop),a(loop,k),
     1                            ifrescal,scal,a(loop,k))
              enddo ! k
            endif
c
c
c           Update the sums-of-squares array ss.
c
            do k = loop,n
              ss(k) = ss(k)-a(loop,k)**2
            enddo ! k
c
c
c           Find the pivot (column with the greatest sum of squares
c           of its entries).
c
            ssmax = 0
            kpiv = loop+1
c
            if(loop .lt. n) then
c
              do k = loop+1,n
c
                if(ss(k) .gt. ssmax) then
                  ssmax = ss(k)
                  kpiv = k
                endif
c
              enddo ! k
c
            endif ! loop .lt. n
c
c
c           Recompute the sums-of-squares and the pivot
c           when ssmax first falls below
c           sqrt((1000*feps)^2) * ssmaxin
c           and when ssmax first falls below
c           ((1000*feps)^2) * ssmaxin.
c
            if(
     1       (ssmax .lt. sqrt((1000*feps)**2) * ssmaxin
     2        .and. nupdate .eq. 0) .or.
     3       (ssmax .lt. ((1000*feps)**2) * ssmaxin
     4        .and. nupdate .eq. 1)
     5      ) then
c
              nupdate = nupdate+1
c
              ssmax = 0
              kpiv = loop+1
c
              if(loop .lt. n) then
c
                do k = loop+1,n
c
                  ss(k) = 0
                  do j = loop+1,m
                    ss(k) = ss(k)+a(j,k)**2
                  enddo ! j
c
                  if(ss(k) .gt. ssmax) then
                    ssmax = ss(k)
                    kpiv = k
                  endif
c
                enddo ! k
c
              endif ! loop .lt. n
c
            endif
c
c
          endif ! loop .lt. m
c
c
        enddo ! loop
c
c
        return
        end
