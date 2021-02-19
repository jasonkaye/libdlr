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
        real*8 snorm,diffsnorm,r1,a(len),dummy,
     1         u(len),v(len),b(len),rnd
        external matvec,matvect
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
        call fill(krank,m,n,a,u)
c
c
c       Calculate the spectral norm of a.
c
        its = 100
c
        call idd_snorm(m,n,matvect,a,dummy,dummy,dummy,
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
          call id_srand(1,rnd)
          b(k) = a(k)+.1d-12*rnd
        enddo ! k
c
c
c       Calculate the spectral norm of a-b.
c
        its = 100
c
        call idd_diffsnorm(m,n,matvect,a,dummy,dummy,dummy,
     1                     matvect,b,dummy,dummy,dummy,
     2                     matvec,a,dummy,dummy,dummy,
     3                     matvec,b,dummy,dummy,dummy,its,diffsnorm,v)
c
c
c       Divide diffsnorm by .5d-13*sqrt(m*n) and display it.
c
        diffsnorm = diffsnorm/(.5d-13*sqrt(r1*m*n))
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
