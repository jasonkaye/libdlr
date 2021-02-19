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
