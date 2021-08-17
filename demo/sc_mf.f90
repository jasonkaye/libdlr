      program sc_mf

      ! Demonstration of discrete Lehmann representation for imaginary
      ! time Green's function with semi-circular density of states,
      ! rho(omega) = sqrt(1-omega^2), using Matsubara frequency
      ! sampling.
      !
      ! The DLR expansion is formed from samples at the DLR Matsubara
      ! frequency nodes, and then evaluated both in imaginary time and
      ! Matsubara frequency domains, where its accuracy is measured.
      !
      ! The Green's function is sampled by evaluating the Green's
      ! function using the Lehmann representation with the known
      ! spectral density. The integral in the Lehmann representation is
      ! computed by high-order numerical integration. This reference
      ! Green's function is also used to test the accuracy of the DLR
      ! expansion.


      implicit none
      integer ntst_it,ntst_mf
      real *8 lambda,eps,beta

      integer nmax

      ! Input parameters

      lambda = 1000.0d0 ! DLR high energy cutoff
      eps = 1.0d-14     ! DLR error tolerance
      nmax = lambda     ! DLR Matsubara frequency cutoff
      
      beta = 1000.0d0   ! Inverse temperature

      ntst_it = 2000    ! # imaginary time test points
      ntst_mf = 2000    ! Matsubara frequency test point cutoff
      

      ! Main subroutine

      call sc_mf_main(lambda,eps,nmax,ntst_it,ntst_mf,beta)


      end program sc_mf


      subroutine sc_mf_main(lambda,eps,nmax,ntst_it,ntst_mf,&
          beta)

      implicit none
      integer nmax,ntst_it,ntst_mf
      real *8 lambda,eps,beta

      integer npt,npo,i,j,r,pg,npg
      integer, allocatable :: mf2cfp(:),dlrmf(:),mf_tst(:)
      real *8 one
      real *8, allocatable :: it_tst(:),dlrrf(:)
      real *8, allocatable :: xgl(:),wgl(:),xgj(:),wgj(:),pbpg(:)
      real *8, allocatable :: gtst_it(:),gtrue_it(:),gc(:)
      complex *16, allocatable :: g(:),mf2cf(:,:)
      complex *16, allocatable :: gtst_mf(:),gtrue_mf(:)

      one = 1.0d0

      ! Get DLR frequencies, Matsubara frequency grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrmf(r))

      call dlr_buildmf(lambda,eps,nmax,-1,r,dlrrf,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients matrix (LU
      ! form)

      allocate(mf2cf(r,r),mf2cfp(r))

      call dlr_mf2cf(nmax,r,dlrrf,dlrmf,-1,mf2cf,mf2cfp)



      ! Initialize reference Green's function evaluator (note: this is
      ! not a DLR-related operation; rather, we are initializing a
      ! special-purpose subroutine to evaluate the Green's function with
      ! a semi-circular spectral density. The initialization and
      ! evaluation subroutines are defined at the bottom of this file.)

      pg = 24
      npg = max(ceiling(log(lambda)/log(2.0d0)),1)
      
      allocate(xgl(pg),wgl(pg),xgj(pg),wgj(pg),pbpg(2*npg+1))

      call gfun_init(pg,npg,pbpg,xgl,wgl,xgj,wgj)


      ! Sample G at Matsubara frequency nodes

      allocate(g(r))

      do i=1,r

        call gfun_mf(pg,npg,pbpg,xgl,wgl,xgj,wgj,beta,dlrmf(i),&
          g(i))

      enddo


      ! Get DLR coefficients of G

      allocate(gc(r))

      call dlr_mfexpnd(r,mf2cf,mf2cfp,g,gc)



      ! Get imaginary time test points in relative format

      allocate(it_tst(ntst_it))

      call eqpts_rel(ntst_it,it_tst)


      ! Evaluate DLR on imaginary time test grid

      allocate(gtst_it(ntst_it))

      do i=1,ntst_it

        call dlr_eval(r,dlrrf,gc,it_tst(i),gtst_it(i))

      enddo


      ! Get Matsubara frequency test points

      allocate(mf_tst(2*ntst_mf+1))

      do i=1,2*ntst_mf+1
        mf_tst(i) = -ntst_mf+i-1
      enddo


      ! Evaluate DLR on Matsubara frequency test grid

      allocate(gtst_mf(2*ntst_mf+1))

      do i=1,2*ntst_mf+1

        call dlr_mf_eval(r,dlrrf,-1,gc,mf_tst(i),gtst_mf(i))

      enddo



      ! Evaluate true Green's function on imaginary time test grid

      allocate(gtrue_it(ntst_it))

      do i=1,ntst_it

        call gfun_it(pg,npg,pbpg,xgl,wgl,xgj,wgj,beta,it_tst(i),&
          gtrue_it(i))

      enddo


      ! Evaluate true Green's function on Matsubara frequency test grid

      allocate(gtrue_mf(2*ntst_mf+1))

      do i=1,2*ntst_mf+1

        call gfun_mf(pg,npg,pbpg,xgl,wgl,xgj,wgj,beta,mf_tst(i),&
          gtrue_mf(i))

      enddo


      ! Output errors

      write(6,*) ''
      write(6,*) 'lambda = ',lambda
      write(6,*) 'epsilon = ',eps
      write(6,*) 'beta = ',beta
      write(6,*) ''
      write(6,*) 'DLR rank = ',r
      write(6,*) ''
      write(6,*) 'Imag time max error = ',maxval(abs(gtst_it-gtrue_it))
      write(6,*) 'Mats freq max error = ',maxval(abs(gtst_mf-gtrue_mf))
      write(6,*) ''

      end subroutine sc_mf_main



      subroutine gfun_init(n,np,pbp,xgl,wgl,xgj,wgj)

      ! Initialize subroutines gfun_it and gfun_mf, which evaluate
      ! Green's function with semi-circular density at a imaginary time
      ! and Matsubara frequency points, respectively

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n)

      integer i
      real *8 one

      one = 1.0d0

      ! Gauss-Legendre and Gauss-Jacobi quadrature nodes and weights

      call cdgqf(n,1,0.0d0,0.0d0,xgl,wgl)
      call cdgqf(n,4,0.5d0,0.0d0,xgj,wgj)


      ! Panel endpoints for composite quadrature rule

      pbp(np+1) = 0*one
      do i=1,np
        pbp(np+i+1) = one/2**(np-i)
      enddo
      pbp(1:np) = -pbp(2*np+1:np+2:-1)

      end subroutine gfun_init


      subroutine gfun_it(n,np,pbp,xgl,wgl,xgj,wgj,beta,t,val)

      ! Evaluate Green's function with semi-circular density at an
      ! imaginary time point using the Lehmann representation. Integral
      ! is computed by high-order composite Gauss-Legendre quadrature,
      ! with Gauss-Jacobi quadrature at endpoint panels to treat square
      ! root singularities.

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta,t,val

      integer ii,jj
      real *8 one,a,b,x,tt
      real *8, external :: kfunf

      one = 1.0d0

      ! Treat 0.5<tau<1, stored in relative format, by symmetry

      tt = abs(t)


      ! Composite Gauss quadrature

      val = 0.0d0
      do ii=2,2*np-1
        a = pbp(ii)
        b = pbp(ii+1)
        do jj=1,n
          x = a+(b-a)*(xgl(jj)+one)/2
          val = val + (b-a)/2*wgl(jj)*kfunf(tt,beta*x)*&
            sqrt(one-x**2)
        enddo
      enddo

      a = one/2
      b = one
      do jj=1,n
        x = a+(b-a)*(xgj(jj)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(jj)*&
          kfunf(tt,beta*x)*sqrt(one+x)
      enddo

      a = -one
      b = -one/2
      do jj=1,n
        x = a+(b-a)*(-xgj(n-jj+1)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(n-jj+1)*&
          kfunf(tt,beta*x)*sqrt(one-x)
      enddo


      end subroutine gfun_it


      subroutine gfun_mf(n,np,pbp,xgl,wgl,xgj,wgj,beta,m,val)

      ! Evaluate Green's function with semi-circular density at a
      ! Matsubara frequency point using the Lehmann representation. Integral
      ! is computed by high-order composite Gauss-Legendre quadrature,
      ! with Gauss-Jacobi quadrature at endpoint panels to treat square
      ! root singularities.

      implicit none
      integer n,np,m
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta
      complex *16 val

      integer ii,jj
      real *8 one,a,b,x
      complex *16, external :: kfunmf

      one = 1.0d0

      val = 0.0d0
      do ii=2,2*np-1
        a = pbp(ii)
        b = pbp(ii+1)
        do jj=1,n
          x = a+(b-a)*(xgl(jj)+one)/2
          val = val + (b-a)/2*wgl(jj)*kfunmf(2*m+1,beta*x)*&
            sqrt(one-x**2)
        enddo
      enddo

      a = one/2
      b = one
      do jj=1,n
        x = a+(b-a)*(xgj(jj)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(jj)*&
          kfunmf(2*m+1,beta*x)*sqrt(one+x)
      enddo

      a = -one
      b = -one/2
      do jj=1,n
        x = a+(b-a)*(-xgj(n-jj+1)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(n-jj+1)*&
          kfunmf(2*m+1,beta*x)*sqrt(one-x)
      enddo

      end subroutine gfun_mf
