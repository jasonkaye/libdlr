      !
      !
      ! This file contains subroutines from the id_dist library used to
      ! implement pivoted QR. It also contains a subroutine,
      ! ind_rearrange, which is not in the original library.
      !
      !
       
        
        
        subroutine iddp_qrpiv(eps,m,n,a,krank,ind,ss)
!
!       computes the pivoted QR decomposition
!       of the matrix input into a, using Householder transformations,
!       _i.e._, transforms the matrix a from its input value in
!       to the matrix out with entry
!
!                               m
!       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k),
!                              l=1
!
!       for all j = 1, ..., krank, and k = 1, ..., n,
!
!       where in = the a from before the routine runs,
!       out = the a from after the routine runs,
!       out(j,k) = 0 when j > k (so that out is triangular),
!       q(1:m,1), ..., q(1:m,krank) are orthonormal,
!       indprod is the product of the permutations given by ind,
!       (as computable via the routine permmult,
!       with the permutation swapping 1 and ind(1) taken leftmost
!       in the product, that swapping 2 and ind(2) taken next leftmost,
!       ..., that swapping krank and ind(krank) taken rightmost),
!       and with the matrix out satisfying
!
!                   krank
!       in(j,k)  =  Sigma  q(j,l) * out(l,indprod(k))  +  epsilon(j,k),
!                    l=1
!
!       for all j = 1, ..., m, and k = 1, ..., n,
!
!       for some matrix epsilon such that
!       the root-sum-square of the entries of epsilon
!       <= the root-sum-square of the entries of in * eps.
!       Well, technically, this routine outputs the Householder vectors
!       (or, rather, their second through last entries)
!       in the part of a that is supposed to get zeroed, that is,
!       in a(j,k) with m >= j > k >= 1.
!
!       input:
!       eps -- relative precision of the resulting QR decomposition
!       m -- first dimension of a and q
!       n -- second dimension of a
!       a -- matrix whose QR decomposition gets computed
!
!       output:
!       a -- triangular (R) factor in the QR decompositon
!            of the matrix input into the same storage locations, 
!            with the Householder vectors stored in the part of a
!            that would otherwise consist entirely of zeroes, that is,
!            in a(j,k) with m >= j > k >= 1
!       krank -- numerical rank
!       ind(k) -- index of the k^th pivot vector;
!                 the following code segment will correctly rearrange
!                 the product b of q and the upper triangle of out
!                 so that b matches the input matrix in
!                 to relative precision eps:
!
!                 copy the non-rearranged product of q and out into b
!                 set k to krank
!                 [start of loop]
!                   swap b(1:m,k) and b(1:m,ind(k))
!                   decrement k by 1
!                 if k > 0, then go to [start of loop]
!
!       work:
!       ss -- must be at least n real*8 words long
!
!       _N.B._: This routine outputs the Householder vectors
!       (or, rather, their second through last entries)
!       in the part of a that is supposed to get zeroed, that is,
!       in a(j,k) with m >= j > k >= 1.
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        integer n,m,ind(n),krank,k,j,kpiv,mm,nupdate,ifrescal
        real*8 a(m,n),ss(n),eps,feps,ssmax,scal,ssmaxin,rswap
!
!
        feps = .1d-16
!
!
!       Compute the sum of squares of the entries in each column of a,
!       the maximum of all such sums, and find the first pivot
!       (column with the greatest such sum).
!
        ssmax = 0
        kpiv = 1
!
        do k = 1,n
!
          ss(k) = 0
          do j = 1,m
            ss(k) = ss(k)+a(j,k)**2
          enddo ! j
!
          if(ss(k) .gt. ssmax) then
            ssmax = ss(k)
            kpiv = k
          endif
!
        enddo ! k
!
        ssmaxin = ssmax
!
        nupdate = 0
!
!
!       While ssmax > eps**2*ssmaxin, krank < m, and krank < n,
!       do the following block of code,
!       which ends at the statement labeled 2000.
!
        krank = 0
 1000   continue
!
        if(ssmax .le. eps**2*ssmaxin
     1   .or. krank .ge. m .or. krank .ge. n) goto 2000
        krank = krank+1
!
!
          mm = m-krank+1
!
!
!         Perform the pivoting.
!
          ind(krank) = kpiv
!
!         Swap a(1:m,krank) and a(1:m,kpiv).
!
          do j = 1,m
            rswap = a(j,krank)
            a(j,krank) = a(j,kpiv)
            a(j,kpiv) = rswap
          enddo ! j
!
!         Swap ss(krank) and ss(kpiv).
!
          rswap = ss(krank)
          ss(krank) = ss(kpiv)
          ss(kpiv) = rswap
!
!
          if(krank .lt. m) then
!
!
!           Compute the data for the Householder transformation
!           which will zero a(krank+1,krank), ..., a(m,krank)
!           when applied to a, replacing a(krank,krank)
!           with the first entry of the result of the application
!           of the Householder matrix to a(krank:m,krank),
!           and storing entries 2 to mm of the Householder vector
!           in a(krank+1,krank), ..., a(m,krank)
!           (which otherwise would get zeroed upon application
!           of the Householder transformation).
!
            call idd_house(mm,a(krank,krank),a(krank,krank),
     1                     a(krank+1,krank),scal)
            ifrescal = 0
!
!
!           Apply the Householder transformation
!           to the lower right submatrix of a
!           with upper leftmost entry at position (krank,krank+1).
!
            if(krank .lt. n) then
              do k = krank+1,n
                call idd_houseapp(mm,a(krank+1,krank),a(krank,k),
     1                            ifrescal,scal,a(krank,k))
              enddo ! k
            endif
!
!
!           Update the sums-of-squares array ss.
!
            do k = krank,n
              ss(k) = ss(k)-a(krank,k)**2
            enddo ! k
!
!
!           Find the pivot (column with the greatest sum of squares
!           of its entries).
!
            ssmax = 0
            kpiv = krank+1
!
            if(krank .lt. n) then
!
              do k = krank+1,n
!
                if(ss(k) .gt. ssmax) then
                  ssmax = ss(k)
                  kpiv = k
                endif
!
              enddo ! k
!
            endif ! krank .lt. n
!
!
!           Recompute the sums-of-squares and the pivot
!           when ssmax first falls below
!           sqrt((1000*feps)^2) * ssmaxin
!           and when ssmax first falls below
!           ((1000*feps)^2) * ssmaxin.
!
            if(
     1       (ssmax .lt. sqrt((1000*feps)**2) * ssmaxin
     2        .and. nupdate .eq. 0) .or.
     3       (ssmax .lt. ((1000*feps)**2) * ssmaxin
     4        .and. nupdate .eq. 1)
     5      ) then
!
              nupdate = nupdate+1
!
              ssmax = 0
              kpiv = krank+1
!
              if(krank .lt. n) then
!
                do k = krank+1,n
!
                  ss(k) = 0
                  do j = krank+1,m
                    ss(k) = ss(k)+a(j,k)**2
                  enddo ! j
!
                  if(ss(k) .gt. ssmax) then
                    ssmax = ss(k)
                    kpiv = k
                  endif
!
                enddo ! k
!
              endif ! krank .lt. n
!
            endif
!
!
          endif ! krank .lt. m
!
!
        goto 1000
 2000   continue
!
!
        return
        end
!
!
!
!
        subroutine iddr_qrpiv(m,n,a,krank,ind,ss)
!
!       computes the pivoted QR decomposition
!       of the matrix input into a, using Householder transformations,
!       _i.e._, transforms the matrix a from its input value in
!       to the matrix out with entry
!
!                               m
!       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k),
!                              l=1
!
!       for all j = 1, ..., krank, and k = 1, ..., n,
!
!       where in = the a from before the routine runs,
!       out = the a from after the routine runs,
!       out(j,k) = 0 when j > k (so that out is triangular),
!       q(1:m,1), ..., q(1:m,krank) are orthonormal,
!       indprod is the product of the permutations given by ind,
!       (as computable via the routine permmult,
!       with the permutation swapping 1 and ind(1) taken leftmost
!       in the product, that swapping 2 and ind(2) taken next leftmost,
!       ..., that swapping krank and ind(krank) taken rightmost),
!       and with the matrix out satisfying
!
!                  min(krank,m,n)
!       in(j,k)  =     Sigma      q(j,l) * out(l,indprod(k))
!                       l=1
!
!                +  epsilon(j,k),
!
!       for all j = 1, ..., m, and k = 1, ..., n,
!
!       for some matrix epsilon whose norm is (hopefully) minimized
!       by the pivoting procedure.
!       Well, technically, this routine outputs the Householder vectors
!       (or, rather, their second through last entries)
!       in the part of a that is supposed to get zeroed, that is,
!       in a(j,k) with m >= j > k >= 1.
!
!       input:
!       m -- first dimension of a and q
!       n -- second dimension of a
!       a -- matrix whose QR decomposition gets computed
!       krank -- desired rank of the output matrix
!                (please note that if krank > m or krank > n,
!                then the rank of the output matrix will be
!                less than krank)
!
!       output:
!       a -- triangular (R) factor in the QR decompositon
!            of the matrix input into the same storage locations, 
!            with the Householder vectors stored in the part of a
!            that would otherwise consist entirely of zeroes, that is,
!            in a(j,k) with m >= j > k >= 1
!       ind(k) -- index of the k^th pivot vector;
!                 the following code segment will correctly rearrange
!                 the product b of q and the upper triangle of out
!                 so that b best matches the input matrix in:
!
!                 copy the non-rearranged product of q and out into b
!                 set k to krank
!                 [start of loop]
!                   swap b(1:m,k) and b(1:m,ind(k))
!                   decrement k by 1
!                 if k > 0, then go to [start of loop]
!
!       work:
!       ss -- must be at least n real*8 words long
!
!       _N.B._: This routine outputs the Householder vectors
!       (or, rather, their second through last entries)
!       in the part of a that is supposed to get zeroed, that is,
!       in a(j,k) with m >= j > k >= 1.
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        integer n,m,ind(n),krank,k,j,kpiv,mm,nupdate,ifrescal,
     1          loops,loop
        real*8 a(m,n),ss(n),ssmax,scal,ssmaxin,rswap,feps
!
!
        feps = .1d-16
!
!
!       Compute the sum of squares of the entries in each column of a,
!       the maximum of all such sums, and find the first pivot
!       (column with the greatest such sum).
!
        ssmax = 0
        kpiv = 1
!
        do k = 1,n
!
          ss(k) = 0
          do j = 1,m
            ss(k) = ss(k)+a(j,k)**2
          enddo ! j
!
          if(ss(k) .gt. ssmax) then
            ssmax = ss(k)
            kpiv = k
          endif
!
        enddo ! k
!
        ssmaxin = ssmax
!
        nupdate = 0
!
!
!       Set loops = min(krank,m,n).
!
        loops = krank
        if(m .lt. loops) loops = m
        if(n .lt. loops) loops = n
!
        do loop = 1,loops
!
!
          mm = m-loop+1
!
!
!         Perform the pivoting.
!
          ind(loop) = kpiv
!
!         Swap a(1:m,loop) and a(1:m,kpiv).
!
          do j = 1,m
            rswap = a(j,loop)
            a(j,loop) = a(j,kpiv)
            a(j,kpiv) = rswap
          enddo ! j
!
!         Swap ss(loop) and ss(kpiv).
!
          rswap = ss(loop)
          ss(loop) = ss(kpiv)
          ss(kpiv) = rswap
!
!
          if(loop .lt. m) then
!
!
!           Compute the data for the Householder transformation
!           which will zero a(loop+1,loop), ..., a(m,loop)
!           when applied to a, replacing a(loop,loop)
!           with the first entry of the result of the application
!           of the Householder matrix to a(loop:m,loop),
!           and storing entries 2 to mm of the Householder vector
!           in a(loop+1,loop), ..., a(m,loop)
!           (which otherwise would get zeroed upon application
!           of the Householder transformation).
!
            call idd_house(mm,a(loop,loop),a(loop,loop),
     1                     a(loop+1,loop),scal)
            ifrescal = 0
!
!
!           Apply the Householder transformation
!           to the lower right submatrix of a
!           with upper leftmost entry at position (loop,loop+1).
!
            if(loop .lt. n) then
              do k = loop+1,n
                call idd_houseapp(mm,a(loop+1,loop),a(loop,k),
     1                            ifrescal,scal,a(loop,k))
              enddo ! k
            endif
!
!
!           Update the sums-of-squares array ss.
!
            do k = loop,n
              ss(k) = ss(k)-a(loop,k)**2
            enddo ! k
!
!
!           Find the pivot (column with the greatest sum of squares
!           of its entries).
!
            ssmax = 0
            kpiv = loop+1
!
            if(loop .lt. n) then
!
              do k = loop+1,n
!
                if(ss(k) .gt. ssmax) then
                  ssmax = ss(k)
                  kpiv = k
                endif
!
              enddo ! k
!
            endif ! loop .lt. n
!
!
!           Recompute the sums-of-squares and the pivot
!           when ssmax first falls below
!           sqrt((1000*feps)^2) * ssmaxin
!           and when ssmax first falls below
!           ((1000*feps)^2) * ssmaxin.
!
            if(
     1       (ssmax .lt. sqrt((1000*feps)**2) * ssmaxin
     2        .and. nupdate .eq. 0) .or.
     3       (ssmax .lt. ((1000*feps)**2) * ssmaxin
     4        .and. nupdate .eq. 1)
     5      ) then
!
              nupdate = nupdate+1
!
              ssmax = 0
              kpiv = loop+1
!
              if(loop .lt. n) then
!
                do k = loop+1,n
!
                  ss(k) = 0
                  do j = loop+1,m
                    ss(k) = ss(k)+a(j,k)**2
                  enddo ! j
!
                  if(ss(k) .gt. ssmax) then
                    ssmax = ss(k)
                    kpiv = k
                  endif
!
                enddo ! k
!
              endif ! loop .lt. n
!
            endif
!
!
          endif ! loop .lt. m
!
!
        enddo ! loop
!
!
        return
        end



        subroutine idd_houseapp(n,vn,u,ifrescal,scal,v)
!
!       applies the Householder matrix
!       identity_matrix - scal * vn * transpose(vn)
!       to the vector u, yielding the vector v;
!
!       scal = 2/(1 + vn(2)^2 + ... + vn(n)^2)
!       when vn(2), ..., vn(n) don't all vanish;
!
!       scal = 0
!       when vn(2), ..., vn(n) do all vanish
!       (including when n = 1).
!
!       input:
!       n -- size of vn, u, and v, though the indexing on vn goes
!            from 2 to n
!       vn -- components 2 to n of the Householder vector vn;
!             vn(1) is assumed to be 1 
!       u -- vector to be transformed
!       ifrescal -- set to 1 to recompute scal from vn(2), ..., vn(n);
!                   set to 0 to use scal as input
!       scal -- see the entry for ifrescal in the decription
!               of the input
!
!       output:
!       scal -- see the entry for ifrescal in the decription
!               of the input
!       v -- result of applying the Householder matrix to u;
!            it's O.K. to have v be the same as u
!            in order to apply the matrix to the vector in place
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        save
        integer n,k,ifrescal
        real*8 vn(2:*),scal,u(n),v(n),fact,sum
!
!
!       Get out of this routine if n = 1.
!
        if(n .eq. 1) then
          v(1) = u(1)
          return
        endif
!
!
        if(ifrescal .eq. 1) then
!
!
!         Calculate (vn(2))^2 + ... + (vn(n))^2.
!
          sum = 0
          do k = 2,n
            sum = sum+vn(k)**2
          enddo ! k
!
!
!         Calculate scal.
!
          if(sum .eq. 0) scal = 0
          if(sum .ne. 0) scal = 2/(1+sum)
!
!
        endif
!
!
!       Calculate fact = scal * transpose(vn) * u.
!
        fact = u(1)
!
        do k = 2,n
          fact = fact+vn(k)*u(k)
        enddo ! k
!
        fact = fact*scal
!
!
!       Subtract fact*vn from u, yielding v.
!      
        v(1) = u(1) - fact
!
        do k = 2,n
          v(k) = u(k) - fact*vn(k)
        enddo ! k
!
!
        return
        end
!
!
!
!
        subroutine idd_house(n,x,rss,vn,scal)
!
!       constructs the vector vn with vn(1) = 1
!       and the scalar scal such that
!       H := identity_matrix - scal * vn * transpose(vn) is orthogonal
!       and Hx = +/- e_1 * the root-sum-square of the entries of x
!       (H is the Householder matrix corresponding to x).
!
!       input:
!       n -- size of x and vn, though the indexing on vn goes
!            from 2 to n
!       x -- vector to reflect into its first component
!
!       output:
!       rss -- first entry of the vector resulting from the application
!              of the Householder matrix to x;
!              its absolute value is the root-sum-square
!              of the entries of x
!       vn -- entries 2 to n of the Householder vector vn;
!             vn(1) is assumed to be 1
!       scal -- scalar multiplying vn * transpose(vn);
!
!               scal = 2/(1 + vn(2)^2 + ... + vn(n)^2)
!               when vn(2), ..., vn(n) don't all vanish;
!
!               scal = 0
!               when vn(2), ..., vn(n) do all vanish
!               (including when n = 1)
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        save
        integer n,k
        real*8 x(n),rss,sum,v1,scal,vn(2:*),x1
!
!
        x1 = x(1)
!
!
!       Get out of this routine if n = 1.
!
        if(n .eq. 1) then
          rss = x1
          scal = 0
          return
        endif
!
!
!       Calculate (x(2))^2 + ... (x(n))^2
!       and the root-sum-square value of the entries in x.
!
!
        sum = 0
        do k = 2,n
          sum = sum+x(k)**2
        enddo ! k
!
!
!       Get out of this routine if sum = 0;
!       flag this case as such by setting v(2), ..., v(n) all to 0.
!
        if(sum .eq. 0) then
!
          rss = x1
          do k = 2,n
            vn(k) = 0
          enddo ! k
          scal = 0
!
          return
!
        endif
!
!
        rss = x1**2 + sum
        rss = sqrt(rss)
!
!
!       Determine the first component v1
!       of the unnormalized Householder vector
!       v = x - rss * (1 0 0 ... 0 0)^T.
!
!       If x1 <= 0, then form x1-rss directly,
!       since that expression cannot involve any cancellation.
!
        if(x1 .le. 0) v1 = x1-rss
!
!       If x1 > 0, then use the fact that
!       x1-rss = -sum / (x1+rss),
!       in order to avoid potential cancellation.
!
        if(x1 .gt. 0) v1 = -sum / (x1+rss)
!
!
!       Compute the vector vn and the scalar scal such that vn(1) = 1
!       in the Householder transformation
!       identity_matrix - scal * vn * transpose(vn).
!
        do k = 2,n
          vn(k) = x(k)/v1
        enddo ! k
!
!       scal = 2
!            / ( vn(1)^2 + vn(2)^2 + ... + vn(n)^2 )
!
!            = 2
!            / ( 1 + vn(2)^2 + ... + vn(n)^2 )
!
!            = 2*v(1)^2
!            / ( v(1)^2 + (v(1)*vn(2))^2 + ... + (v(1)*vn(n))^2 )
!
!            = 2*v(1)^2
!            / ( v(1)^2 + (v(2)^2 + ... + v(n)^2) )
!
        scal = 2*v1**2 / (v1**2+sum)
!
!
        return
        end



        subroutine idzr_qrpiv(m,n,a,krank,ind,ss)
!
!       computes the pivoted QR decomposition
!       of the matrix input into a, using Householder transformations,
!       _i.e._, transforms the matrix a from its input value in
!       to the matrix out with entry
!
!                               m
!       out(j,indprod(k))  =  Sigma  q(l,j) * in(l,k),
!                              l=1
!
!       for all j = 1, ..., krank, and k = 1, ..., n,
!
!       where in = the a from before the routine runs,
!       out = the a from after the routine runs,
!       out(j,k) = 0 when j > k (so that out is triangular),
!       q(1:m,1), ..., q(1:m,krank) are orthonormal,
!       indprod is the product of the permutations given by ind,
!       (as computable via the routine permmult,
!       with the permutation swapping 1 and ind(1) taken leftmost
!       in the product, that swapping 2 and ind(2) taken next leftmost,
!       ..., that swapping krank and ind(krank) taken rightmost),
!       and with the matrix out satisfying
!
!                  min(m,n,krank)
!       in(j,k)  =     Sigma      q(j,l) * out(l,indprod(k))
!                       l=1
!
!                +  epsilon(j,k),
!
!       for all j = 1, ..., m, and k = 1, ..., n,
!
!       for some matrix epsilon whose norm is (hopefully) minimized
!       by the pivoting procedure.
!       Well, technically, this routine outputs the Householder vectors
!       (or, rather, their second through last entries)
!       in the part of a that is supposed to get zeroed, that is,
!       in a(j,k) with m >= j > k >= 1.
!
!       input:
!       m -- first dimension of a and q
!       n -- second dimension of a
!       a -- matrix whose QR decomposition gets computed
!       krank -- desired rank of the output matrix
!                (please note that if krank > m or krank > n,
!                then the rank of the output matrix will be
!                less than krank)
!
!       output:
!       a -- triangular (R) factor in the QR decompositon
!            of the matrix input into the same storage locations, 
!            with the Householder vectors stored in the part of a
!            that would otherwise consist entirely of zeroes, that is,
!            in a(j,k) with m >= j > k >= 1
!       ind(k) -- index of the k^th pivot vector;
!                 the following code segment will correctly rearrange
!                 the product b of q and the upper triangle of out
!                 so that b matches the input matrix in
!                 to relative precision eps:
!
!                 copy the non-rearranged product of q and out into b
!                 set k to krank
!                 [start of loop]
!                   swap b(1:m,k) and b(1:m,ind(k))
!                   decrement k by 1
!                 if k > 0, then go to [start of loop]
!
!       work:
!       ss -- must be at least n real*8 words long
!
!       _N.B._: This routine outputs the Householder vectors
!       (or, rather, their second through last entries)
!       in the part of a that is supposed to get zeroed, that is,
!       in a(j,k) with m >= j > k >= 1.
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        integer n,m,ind(n),krank,k,j,kpiv,mm,nupdate,ifrescal,
     1          loops,loop
        real*8 ss(n),ssmax,scal,ssmaxin,rswap,feps
        complex*16 a(m,n),cswap
!
!
        feps = .1d-16
!
!
!       Compute the sum of squares of the entries in each column of a,
!       the maximum of all such sums, and find the first pivot
!       (column with the greatest such sum).
!
        ssmax = 0
        kpiv = 1
!
        do k = 1,n
!
          ss(k) = 0
          do j = 1,m
            ss(k) = ss(k)+a(j,k)*conjg(a(j,k))
          enddo ! j
!
          if(ss(k) .gt. ssmax) then
            ssmax = ss(k)
            kpiv = k
          endif
!
        enddo ! k
!
        ssmaxin = ssmax
!
        nupdate = 0
!
!
!       Set loops = min(krank,m,n).
!
        loops = krank
        if(m .lt. loops) loops = m
        if(n .lt. loops) loops = n
!
        do loop = 1,loops
!
!
          mm = m-loop+1
!
!
!         Perform the pivoting.
!
          ind(loop) = kpiv
!
!         Swap a(1:m,loop) and a(1:m,kpiv).
!
          do j = 1,m
            cswap = a(j,loop)
            a(j,loop) = a(j,kpiv)
            a(j,kpiv) = cswap
          enddo ! j
!
!         Swap ss(loop) and ss(kpiv).
!
          rswap = ss(loop)
          ss(loop) = ss(kpiv)
          ss(kpiv) = rswap
!
!
          if(loop .lt. m) then
!
!
!           Compute the data for the Householder transformation
!           which will zero a(loop+1,loop), ..., a(m,loop)
!           when applied to a, replacing a(loop,loop)
!           with the first entry of the result of the application
!           of the Householder matrix to a(loop:m,loop),
!           and storing entries 2 to mm of the Householder vector
!           in a(loop+1,loop), ..., a(m,loop)
!           (which otherwise would get zeroed upon application
!           of the Householder transformation).
!
            call idz_house(mm,a(loop,loop),a(loop,loop),
     1                     a(loop+1,loop),scal)
            ifrescal = 0
!
!
!           Apply the Householder transformation
!           to the lower right submatrix of a
!           with upper leftmost entry at position (loop,loop+1).
!
            if(loop .lt. n) then
              do k = loop+1,n
                call idz_houseapp(mm,a(loop+1,loop),a(loop,k),
     1                            ifrescal,scal,a(loop,k))
              enddo ! k
            endif
!
!
!           Update the sums-of-squares array ss.
!
            do k = loop,n
              ss(k) = ss(k)-a(loop,k)*conjg(a(loop,k))
            enddo ! k
!
!
!           Find the pivot (column with the greatest sum of squares
!           of its entries).
!
            ssmax = 0
            kpiv = loop+1
!
            if(loop .lt. n) then
!
              do k = loop+1,n
!
                if(ss(k) .gt. ssmax) then
                  ssmax = ss(k)
                  kpiv = k
                endif
!
              enddo ! k
!
            endif ! loop .lt. n
!
!
!           Recompute the sums-of-squares and the pivot
!           when ssmax first falls below
!           sqrt((1000*feps)^2) * ssmaxin
!           and when ssmax first falls below
!           ((1000*feps)^2) * ssmaxin.
!
            if(
     1       (ssmax .lt. sqrt((1000*feps)**2) * ssmaxin
     2        .and. nupdate .eq. 0) .or.
     3       (ssmax .lt. ((1000*feps)**2) * ssmaxin
     4        .and. nupdate .eq. 1)
     5      ) then
!
              nupdate = nupdate+1
!
              ssmax = 0
              kpiv = loop+1
!
              if(loop .lt. n) then
!
                do k = loop+1,n
!
                  ss(k) = 0
                  do j = loop+1,m
                    ss(k) = ss(k)+a(j,k)*conjg(a(j,k))
                  enddo ! j
!
                  if(ss(k) .gt. ssmax) then
                    ssmax = ss(k)
                    kpiv = k
                  endif
!
                enddo ! k
!
              endif ! loop .lt. n
!
            endif
!
!
          endif ! loop .lt. m
!
!
        enddo ! loop
!
!
        return
        end



        subroutine idz_houseapp(n,vn,u,ifrescal,scal,v)
!
!       applies the Householder matrix
!       identity_matrix - scal * vn * adjoint(vn)
!       to the vector u, yielding the vector v;
!
!       scal = 2/(1 + |vn(2)|^2 + ... + |vn(n)|^2)
!       when vn(2), ..., vn(n) don't all vanish;
!
!       scal = 0
!       when vn(2), ..., vn(n) do all vanish
!       (including when n = 1).
!
!       input:
!       n -- size of vn, u, and v, though the indexing on vn goes
!            from 2 to n
!       vn -- components 2 to n of the Householder vector vn;
!             vn(1) is assumed to be 1 
!       u -- vector to be transformed
!       ifrescal -- set to 1 to recompute scal from vn(2), ..., vn(n);
!                   set to 0 to use scal as input
!       scal -- see the entry for ifrescal in the decription
!               of the input
!
!       output:
!       scal -- see the entry for ifrescal in the decription
!               of the input
!       v -- result of applying the Householder matrix to u;
!            it's O.K. to have v be the same as u
!            in order to apply the matrix to the vector in place
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        save
        integer n,k,ifrescal
        real*8 scal,sum
        complex*16 vn(2:*),u(n),v(n),fact
!
!
!       Get out of this routine if n = 1.
!
        if(n .eq. 1) then
          v(1) = u(1)
          return
        endif
!
!
        if(ifrescal .eq. 1) then
!
!
!         Calculate |vn(2)|^2 + ... + |vn(n)|^2.
!
          sum = 0
          do k = 2,n
            sum = sum+vn(k)*conjg(vn(k))
          enddo ! k
!
!
!         Calculate scal.
!
          if(sum .eq. 0) scal = 0
          if(sum .ne. 0) scal = 2/(1+sum)
!
!
        endif
!
!
!       Calculate fact = scal * adjoint(vn) * u.
!
        fact = u(1)
!
        do k = 2,n
          fact = fact+conjg(vn(k))*u(k)
        enddo ! k
!
        fact = fact*scal
!
!
!       Subtract fact*vn from u, yielding v.
!      
        v(1) = u(1) - fact
!
        do k = 2,n
          v(k) = u(k) - fact*vn(k)
        enddo ! k
!
!
        return
        end
!
!
!
!
        subroutine idz_house(n,x,css,vn,scal)
!
!       constructs the vector vn with vn(1) = 1,
!       and the scalar scal, such that the obviously self-adjoint
!       H := identity_matrix - scal * vn * adjoint(vn) is unitary,
!       the absolute value of the first entry of Hx
!       is the root-sum-square of the entries of x,
!       and all other entries of Hx are zero
!       (H is the Householder matrix corresponding to x).
!
!       input:
!       n -- size of x and vn, though the indexing on vn goes
!            from 2 to n
!       x -- vector to reflect into its first component
!
!       output:
!       css -- root-sum-square of the entries of x * the phase of x(1)
!       vn -- entries 2 to n of the Householder vector vn;
!             vn(1) is assumed to be 1
!       scal -- scalar multiplying vn * adjoint(vn);
!
!               scal = 2/(1 + |vn(2)|^2 + ... + |vn(n)|^2)
!               when vn(2), ..., vn(n) don't all vanish;
!
!               scal = 0
!               when vn(2), ..., vn(n) do all vanish
!               (including when n = 1)
!
!       reference:
!       Golub and Van Loan, "Matrix Computations," 3rd edition,
!            Johns Hopkins University Press, 1996, Chapter 5.
!
        implicit none
        save
        integer n,k
        real*8 scal,test,rss,sum
        complex*16 x(n),v1,vn(2:*),x1,phase,css
!
!
        x1 = x(1)
!
!
!       Get out of this routine if n = 1.
!
        if(n .eq. 1) then
          css = x1
          scal = 0
          return
        endif
!
!
!       Calculate |x(2)|^2 + ... |x(n)|^2
!       and the root-sum-square value of the entries in x.
!
!
        sum = 0
        do k = 2,n
          sum = sum+x(k)*conjg(x(k))
        enddo ! k
!
!
!       Get out of this routine if sum = 0;
!       flag this case as such by setting v(2), ..., v(n) all to 0.
!
        if(sum .eq. 0) then
!
          css = x1
          do k = 2,n
            vn(k) = 0
          enddo ! k
          scal = 0
!
          return
!
        endif
!
!
        rss = x1*conjg(x1) + sum
        rss = sqrt(rss)
!
!
!       Determine the first component v1
!       of the unnormalized Householder vector
!       v = x - phase(x1) * rss * (1 0 0 ... 0 0)^T.
!
        if(x1 .eq. 0) phase = 1
        if(x1 .ne. 0) phase = x1/abs(x1)
        test = conjg(phase) * x1
        css = phase*rss
!
!       If test <= 0, then form x1-phase*rss directly,
!       since that expression cannot involve any cancellation.
!
        if(test .le. 0) v1 = x1-phase*rss
!
!       If test > 0, then use the fact that
!       x1-phase*rss = -phase*sum / ((phase)^* * x1 + rss),
!       in order to avoid potential cancellation.
!
        if(test .gt. 0) v1 = -phase*sum / (conjg(phase)*x1+rss)
!
!
!       Compute the vector vn and the scalar scal such that vn(1) = 1
!       in the Householder transformation
!       identity_matrix - scal * vn * adjoint(vn).
!
        do k = 2,n
          vn(k) = x(k)/v1
        enddo ! k
!
!       scal = 2
!            / ( |vn(1)|^2 + |vn(2)|^2 + ... + |vn(n)|^2 )
!
!            = 2
!            / ( 1 + |vn(2)|^2 + ... + |vn(n)|^2 )
!
!            = 2*|v(1)|^2
!            / ( |v(1)|^2 + |v(1)*vn(2)|^2 + ... + |v(1)*vn(n)|^2 )
!
!            = 2*|v(1)|^2
!            / ( |v(1)|^2 + (|v(2)|^2 + ... + |v(n)|^2) )
!
        scal = 2*v1*conjg(v1) / (v1*conjg(v1)+sum)
!
!
        rss = phase*rss
!
!
        return
        end
!
!




      subroutine ind_rearrange(n,krank,ind)

      ! Rearrange output ind from iddp_qrpiv or iddr_qrpiv to give list
      ! of krank columns selected by a pivoted QR process
      !
      ! Output overwrites first krank entries of ind
      !
      ! NOTE: This subroutine is not in the original ID library.

      implicit none
      integer n,krank,ind(n)

      integer k,iswap
      integer, allocatable :: tmp(:)

      allocate(tmp(n))

      do k = 1,n
        tmp(k) = k
      enddo
      
      do k = 1,krank
      
        iswap = tmp(k)
        tmp(k) = tmp(ind(k))
        tmp(ind(k)) = iswap
      
      enddo

      ind(1:krank) = tmp(1:krank)
     
      end subroutine ind_rearrange
