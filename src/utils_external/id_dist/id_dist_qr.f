      !
      !
      ! This file contains subroutines from the id_dist library used to
      ! implement pivoted QR
      !
      !
       
        
        
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



        subroutine idd_houseapp(n,vn,u,ifrescal,scal,v)
c
c       applies the Householder matrix
c       identity_matrix - scal * vn * transpose(vn)
c       to the vector u, yielding the vector v;
c
c       scal = 2/(1 + vn(2)^2 + ... + vn(n)^2)
c       when vn(2), ..., vn(n) don't all vanish;
c
c       scal = 0
c       when vn(2), ..., vn(n) do all vanish
c       (including when n = 1).
c
c       input:
c       n -- size of vn, u, and v, though the indexing on vn goes
c            from 2 to n
c       vn -- components 2 to n of the Householder vector vn;
c             vn(1) is assumed to be 1 
c       u -- vector to be transformed
c       ifrescal -- set to 1 to recompute scal from vn(2), ..., vn(n);
c                   set to 0 to use scal as input
c       scal -- see the entry for ifrescal in the decription
c               of the input
c
c       output:
c       scal -- see the entry for ifrescal in the decription
c               of the input
c       v -- result of applying the Householder matrix to u;
c            it's O.K. to have v be the same as u
c            in order to apply the matrix to the vector in place
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        save
        integer n,k,ifrescal
        real*8 vn(2:*),scal,u(n),v(n),fact,sum
c
c
c       Get out of this routine if n = 1.
c
        if(n .eq. 1) then
          v(1) = u(1)
          return
        endif
c
c
        if(ifrescal .eq. 1) then
c
c
c         Calculate (vn(2))^2 + ... + (vn(n))^2.
c
          sum = 0
          do k = 2,n
            sum = sum+vn(k)**2
          enddo ! k
c
c
c         Calculate scal.
c
          if(sum .eq. 0) scal = 0
          if(sum .ne. 0) scal = 2/(1+sum)
c
c
        endif
c
c
c       Calculate fact = scal * transpose(vn) * u.
c
        fact = u(1)
c
        do k = 2,n
          fact = fact+vn(k)*u(k)
        enddo ! k
c
        fact = fact*scal
c
c
c       Subtract fact*vn from u, yielding v.
c      
        v(1) = u(1) - fact
c
        do k = 2,n
          v(k) = u(k) - fact*vn(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_house(n,x,rss,vn,scal)
c
c       constructs the vector vn with vn(1) = 1
c       and the scalar scal such that
c       H := identity_matrix - scal * vn * transpose(vn) is orthogonal
c       and Hx = +/- e_1 * the root-sum-square of the entries of x
c       (H is the Householder matrix corresponding to x).
c
c       input:
c       n -- size of x and vn, though the indexing on vn goes
c            from 2 to n
c       x -- vector to reflect into its first component
c
c       output:
c       rss -- first entry of the vector resulting from the application
c              of the Householder matrix to x;
c              its absolute value is the root-sum-square
c              of the entries of x
c       vn -- entries 2 to n of the Householder vector vn;
c             vn(1) is assumed to be 1
c       scal -- scalar multiplying vn * transpose(vn);
c
c               scal = 2/(1 + vn(2)^2 + ... + vn(n)^2)
c               when vn(2), ..., vn(n) don't all vanish;
c
c               scal = 0
c               when vn(2), ..., vn(n) do all vanish
c               (including when n = 1)
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        save
        integer n,k
        real*8 x(n),rss,sum,v1,scal,vn(2:*),x1
c
c
        x1 = x(1)
c
c
c       Get out of this routine if n = 1.
c
        if(n .eq. 1) then
          rss = x1
          scal = 0
          return
        endif
c
c
c       Calculate (x(2))^2 + ... (x(n))^2
c       and the root-sum-square value of the entries in x.
c
c
        sum = 0
        do k = 2,n
          sum = sum+x(k)**2
        enddo ! k
c
c
c       Get out of this routine if sum = 0;
c       flag this case as such by setting v(2), ..., v(n) all to 0.
c
        if(sum .eq. 0) then
c
          rss = x1
          do k = 2,n
            vn(k) = 0
          enddo ! k
          scal = 0
c
          return
c
        endif
c
c
        rss = x1**2 + sum
        rss = sqrt(rss)
c
c
c       Determine the first component v1
c       of the unnormalized Householder vector
c       v = x - rss * (1 0 0 ... 0 0)^T.
c
c       If x1 <= 0, then form x1-rss directly,
c       since that expression cannot involve any cancellation.
c
        if(x1 .le. 0) v1 = x1-rss
c
c       If x1 > 0, then use the fact that
c       x1-rss = -sum / (x1+rss),
c       in order to avoid potential cancellation.
c
        if(x1 .gt. 0) v1 = -sum / (x1+rss)
c
c
c       Compute the vector vn and the scalar scal such that vn(1) = 1
c       in the Householder transformation
c       identity_matrix - scal * vn * transpose(vn).
c
        do k = 2,n
          vn(k) = x(k)/v1
        enddo ! k
c
c       scal = 2
c            / ( vn(1)^2 + vn(2)^2 + ... + vn(n)^2 )
c
c            = 2
c            / ( 1 + vn(2)^2 + ... + vn(n)^2 )
c
c            = 2*v(1)^2
c            / ( v(1)^2 + (v(1)*vn(2))^2 + ... + (v(1)*vn(n))^2 )
c
c            = 2*v(1)^2
c            / ( v(1)^2 + (v(2)^2 + ... + v(n)^2) )
c
        scal = 2*v1**2 / (v1**2+sum)
c
c
        return
        end



        subroutine idzr_qrpiv(m,n,a,krank,ind,ss)
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
c                  min(m,n,krank)
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
        integer n,m,ind(n),krank,k,j,kpiv,mm,nupdate,ifrescal,
     1          loops,loop
        real*8 ss(n),ssmax,scal,ssmaxin,rswap,feps
        complex*16 a(m,n),cswap
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
            ss(k) = ss(k)+a(j,k)*conjg(a(j,k))
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
            cswap = a(j,loop)
            a(j,loop) = a(j,kpiv)
            a(j,kpiv) = cswap
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
            call idz_house(mm,a(loop,loop),a(loop,loop),
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
                call idz_houseapp(mm,a(loop+1,loop),a(loop,k),
     1                            ifrescal,scal,a(loop,k))
              enddo ! k
            endif
c
c
c           Update the sums-of-squares array ss.
c
            do k = loop,n
              ss(k) = ss(k)-a(loop,k)*conjg(a(loop,k))
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
                    ss(k) = ss(k)+a(j,k)*conjg(a(j,k))
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



        subroutine idz_houseapp(n,vn,u,ifrescal,scal,v)
c
c       applies the Householder matrix
c       identity_matrix - scal * vn * adjoint(vn)
c       to the vector u, yielding the vector v;
c
c       scal = 2/(1 + |vn(2)|^2 + ... + |vn(n)|^2)
c       when vn(2), ..., vn(n) don't all vanish;
c
c       scal = 0
c       when vn(2), ..., vn(n) do all vanish
c       (including when n = 1).
c
c       input:
c       n -- size of vn, u, and v, though the indexing on vn goes
c            from 2 to n
c       vn -- components 2 to n of the Householder vector vn;
c             vn(1) is assumed to be 1 
c       u -- vector to be transformed
c       ifrescal -- set to 1 to recompute scal from vn(2), ..., vn(n);
c                   set to 0 to use scal as input
c       scal -- see the entry for ifrescal in the decription
c               of the input
c
c       output:
c       scal -- see the entry for ifrescal in the decription
c               of the input
c       v -- result of applying the Householder matrix to u;
c            it's O.K. to have v be the same as u
c            in order to apply the matrix to the vector in place
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        save
        integer n,k,ifrescal
        real*8 scal,sum
        complex*16 vn(2:*),u(n),v(n),fact
c
c
c       Get out of this routine if n = 1.
c
        if(n .eq. 1) then
          v(1) = u(1)
          return
        endif
c
c
        if(ifrescal .eq. 1) then
c
c
c         Calculate |vn(2)|^2 + ... + |vn(n)|^2.
c
          sum = 0
          do k = 2,n
            sum = sum+vn(k)*conjg(vn(k))
          enddo ! k
c
c
c         Calculate scal.
c
          if(sum .eq. 0) scal = 0
          if(sum .ne. 0) scal = 2/(1+sum)
c
c
        endif
c
c
c       Calculate fact = scal * adjoint(vn) * u.
c
        fact = u(1)
c
        do k = 2,n
          fact = fact+conjg(vn(k))*u(k)
        enddo ! k
c
        fact = fact*scal
c
c
c       Subtract fact*vn from u, yielding v.
c      
        v(1) = u(1) - fact
c
        do k = 2,n
          v(k) = u(k) - fact*vn(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idz_house(n,x,css,vn,scal)
c
c       constructs the vector vn with vn(1) = 1,
c       and the scalar scal, such that the obviously self-adjoint
c       H := identity_matrix - scal * vn * adjoint(vn) is unitary,
c       the absolute value of the first entry of Hx
c       is the root-sum-square of the entries of x,
c       and all other entries of Hx are zero
c       (H is the Householder matrix corresponding to x).
c
c       input:
c       n -- size of x and vn, though the indexing on vn goes
c            from 2 to n
c       x -- vector to reflect into its first component
c
c       output:
c       css -- root-sum-square of the entries of x * the phase of x(1)
c       vn -- entries 2 to n of the Householder vector vn;
c             vn(1) is assumed to be 1
c       scal -- scalar multiplying vn * adjoint(vn);
c
c               scal = 2/(1 + |vn(2)|^2 + ... + |vn(n)|^2)
c               when vn(2), ..., vn(n) don't all vanish;
c
c               scal = 0
c               when vn(2), ..., vn(n) do all vanish
c               (including when n = 1)
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        save
        integer n,k
        real*8 scal,test,rss,sum
        complex*16 x(n),v1,vn(2:*),x1,phase,css
c
c
        x1 = x(1)
c
c
c       Get out of this routine if n = 1.
c
        if(n .eq. 1) then
          css = x1
          scal = 0
          return
        endif
c
c
c       Calculate |x(2)|^2 + ... |x(n)|^2
c       and the root-sum-square value of the entries in x.
c
c
        sum = 0
        do k = 2,n
          sum = sum+x(k)*conjg(x(k))
        enddo ! k
c
c
c       Get out of this routine if sum = 0;
c       flag this case as such by setting v(2), ..., v(n) all to 0.
c
        if(sum .eq. 0) then
c
          css = x1
          do k = 2,n
            vn(k) = 0
          enddo ! k
          scal = 0
c
          return
c
        endif
c
c
        rss = x1*conjg(x1) + sum
        rss = sqrt(rss)
c
c
c       Determine the first component v1
c       of the unnormalized Householder vector
c       v = x - phase(x1) * rss * (1 0 0 ... 0 0)^T.
c
        if(x1 .eq. 0) phase = 1
        if(x1 .ne. 0) phase = x1/abs(x1)
        test = conjg(phase) * x1
        css = phase*rss
c
c       If test <= 0, then form x1-phase*rss directly,
c       since that expression cannot involve any cancellation.
c
        if(test .le. 0) v1 = x1-phase*rss
c
c       If test > 0, then use the fact that
c       x1-phase*rss = -phase*sum / ((phase)^* * x1 + rss),
c       in order to avoid potential cancellation.
c
        if(test .gt. 0) v1 = -phase*sum / (conjg(phase)*x1+rss)
c
c
c       Compute the vector vn and the scalar scal such that vn(1) = 1
c       in the Householder transformation
c       identity_matrix - scal * vn * adjoint(vn).
c
        do k = 2,n
          vn(k) = x(k)/v1
        enddo ! k
c
c       scal = 2
c            / ( |vn(1)|^2 + |vn(2)|^2 + ... + |vn(n)|^2 )
c
c            = 2
c            / ( 1 + |vn(2)|^2 + ... + |vn(n)|^2 )
c
c            = 2*|v(1)|^2
c            / ( |v(1)|^2 + |v(1)*vn(2)|^2 + ... + |v(1)*vn(n)|^2 )
c
c            = 2*|v(1)|^2
c            / ( |v(1)|^2 + (|v(2)|^2 + ... + |v(n)|^2) )
c
        scal = 2*v1*conjg(v1) / (v1*conjg(v1)+sum)
c
c
        rss = phase*rss
c
c
        return
        end
c
c
