c
c
c       dependencies: prini, id_rand
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer m,n,krank,its,k
        real*8 snorm,diffsnorm,r1
        complex*16 a(len),dummy,u(len),v(len),b(len),rnd
        external matvec,matveca
c
        r1 = 1
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
c       Fill a with a matrix whose spectral norm is 2.
c
        call fill(krank,m,n,a)
c
c
c       Calculate the spectral norm of a.
c
        its = 100
c
        call idz_snorm(m,n,matveca,a,dummy,dummy,dummy,
     1                 matvec,a,dummy,dummy,dummy,its,snorm,v,u)
c
c
c       Divide snorm by 2 and display it.
c
        snorm = snorm/2
        call prin2('snorm (which should be 1) = *',snorm,1)
c
c
c       Add a little noise to a, obtaining b.
c
        do k = 1,m*n
          call id_srand(2,rnd)
          b(k) = a(k)+.1d-12*(2*rnd-1)
        enddo ! k
c
c
c       Calculate the spectral norm of a-b.
c
        its = 100
c
        call idz_diffsnorm(m,n,matveca,a,dummy,dummy,dummy,
     1                     matveca,b,dummy,dummy,dummy,
     2                     matvec,a,dummy,dummy,dummy,
     3                     matvec,b,dummy,dummy,dummy,its,diffsnorm,v)
c
c
c       Divide diffsnorm by .1d-12*sqrt(m*n) and display it.
c
        diffsnorm = diffsnorm/(.1d-12*sqrt(r1*m*n))
        call prin2('diffsnorm (which should be about 1) = *',
     1             diffsnorm,1)
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
c       krank -- one less than the rank of the matrix to be constructed
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
            l = krank+1
            sum = sum+exp(2*pi*ci*(j-r1)*(l-r1)/m)*sqrt(r1/m)
     1               *exp(2*pi*ci*(k-r1)*(l-r1)/n)*sqrt(r1/n)
     2               *1d-10
c
            a(j,k) = sum*2
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
        subroutine matveca(m,x,n,y,a,p2,p3,p4)
c
c       applies the adjoint of a to x, obtaining y.
c
c       input:
c       m -- first dimension of a, and length of x
c       x -- vector to which a^* is to be applied
c       n -- second dimension of a, and length of y
c       a -- matrix whose adjoint is to be applied to x
c            in order to create y
c       p2 -- dummy input
c       p3 -- dummy input
c       p4 -- dummy input
c
c       output:
c       y -- product of a^* and x
c
        implicit none
        integer m,n,j,k
        complex*16 a(m,n),p2,p3,p4,x(m),y(n),sum
c
c
        do k = 1,n
c
          sum = 0
c
          do j = 1,m
            sum = sum+conjg(a(j,k))*x(j)
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
        complex*16 a(m,n),p2,p3,p4,x(n),y(m),sum
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
