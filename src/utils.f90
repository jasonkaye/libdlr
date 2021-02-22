      !
      !
      ! This file contains common utility subroutines used in the code
      !
      !
     
      
      
      
      subroutine cgl_pt2coef(p,npan,nfun,legf,val,coef)

      ! ----- Convert point value representation on composite
      ! Gauss-Legendre grid to coefficient representation -----

      integer p,npan
      real *8 val(p*npan,nfun),coef(p*npan,nfun),legf(p,p)

      do i=1,npan
        coef((i-1)*p+1:i*p,:) = matmul(legf,val((i-1)*p+1:i*p,:))
      enddo

      end subroutine cgl_pt2coef


      subroutine ind_rearrange(n,krank,ind)

      ! Rearrange output ind from iddp_qrpiv or iddr_qrpiv to give list
      ! of krank columns selected by a pivoted QR process
      !
      ! Output overwrites first krank entries of ind

      implicit none
      integer n,krank,ind(n)

      integer k,iswap
      integer, allocatable :: tmp(:)

      allocate(tmp(n))

      do k = 1,n
        tmp(k) = k
      enddo
      
      do k = 1,krank
      
        ! Swap rnorms(k) and rnorms(list(k)).
      
        iswap = tmp(k)
        tmp(k) = tmp(ind(k))
        tmp(ind(k)) = iswap
      
      enddo

      ind(1:krank) = tmp(1:krank)
     
      end subroutine ind_rearrange
     

      subroutine evalexpfun(x,val)

      ! Evaluate the function (1-e^(-x))/x

      implicit none
      real *8 x,val

      integer i
      real *8 one,x1,c

      one = 1.0d0

      if (x>1.0d-1) then

        val = (one-exp(-x))/x

      else

        val = one
        c = one
        x1 = x
        do i=1,10
          c = c*(i+1)
          val = val + ((-1)**i)*x1/c
          x1 = x1*x
        enddo

      endif

      end subroutine evalexpfun


      subroutine barychebinit(n,x,w)

      ! Get Chebyshev nodes of first kind and corresponding barycentric
      ! Lagrange interpolation weights

      implicit none
      integer n
      real *8 x(n),w(n)

      integer j
      real *8 pi,c

      pi = 4*atan(1.0d0)

      do j=1,n
        c = (2.0d0*j-1)/(2*n)*pi
        x(n-j+1) = cos(c)
        w(n-j+1) = (1-2*mod(j-1,2))*sin(c)
        !w(j) = (-1)**(j-1)*sin(c)
      enddo

      end subroutine barychebinit


      subroutine barycheb(n,x,f,wc,xc,val)

      ! Barycentric Lagrange interpolation at Chebyshev nodes

      implicit none
      integer n
      real *8 x,f(n),wc(n),xc(n),val

      integer j
      real *8 dif,q,num,den

      do j=1,n
        if (x==xc(j)) then
          val = f(j)
          return
        endif
      enddo

      num = 0.0d0
      den = 0.0d0
      do j=1,n

        dif = x-xc(j)
        q = wc(j)/dif
        num = num + q*f(j)
        den = den + q

      enddo
      
      val = num/den

      end subroutine barycheb
