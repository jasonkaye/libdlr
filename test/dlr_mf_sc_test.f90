      program dlr_sc_mf_test

      ! Test discrete Lehmann representation using Green's function
      ! generated from Lehmann representation with semi-circular
      ! density. Recover DLR coefficients from samples of Green's
      ! function at DLR Matsubara frequency points, and then measure the error of the
      ! resulting expansion on a test grid in imaginary time.
      
      implicit none
      integer ntst,nmax
      real *8 lambda,eps,beta
      character :: fb

      integer i

      ! --- Input parameters ---

      lambda = 100 ! Frequency cutoff
      eps = 1.0d-14 ! Desired accuracy
      nmax = lambda ! Matsubara frequency cutoff
      ntst = 1000 ! # test points to check representation of G
      beta = 100 ! Inverse temp: controls support of rho
      fb = 'f' ! Fermion or Boson? (this switch isn't working yet)


      ! --- Call main test subroutine ---

      call dlr_sc_mf_test_main(lambda,eps,nmax,ntst,beta,fb)


      !do i=1,7
      !  lambda = 20.0d0+20*(i-1)
      !  call dlr_sc_test_main(lambda,eps,ntst,beta,fb)
      !enddo

      end program dlr_sc_mf_test


      subroutine dlr_sc_mf_test_main(lambda,eps,nmax,ntst,beta,fb)

      ! Main driver routine for test of DLR basis on Green's function
      ! with semi-circular density

      implicit none
      integer ntst,nmax
      real *8 lambda,eps,beta
      character :: fb

      integer npt,npo,p,nt,no,i,j,rank,info,pg,npg
      integer, allocatable :: oidx(:),dlrmf(:),ipiv(:)
      real *8 one,gtrue,errl2,errlinf,kerr(2),gmax,gl2,gtest
      real *8, allocatable :: kmat(:,:),t(:),om(:),ttst(:),dlrrf(:)
      real *8, allocatable :: xgl(:),wgl(:),xgj(:),wgj(:),pbpg(:)
      complex *16, allocatable :: gdlr(:),mf2cf(:,:)

      one = 1.0d0

      write(6,*) ''
      write(6,*) '---------------- Input parameters ----------------'
      write(6,*) ''
      write(6,*) 'Cutoff lambda              = ',lambda
      write(6,*) 'Error tolerance eps        = ',eps
      write(6,*) 'Matsubara freq cutoff nmax = ',nmax
      write(6,*) 'Inverse temp beta          = ',beta
      write(6,*) '# test points              = ',ntst
      write(6,*) 'Fermion (f) or boson (b)   = ',fb


      ! --- Build DLR basis, grid, transform matrix ---

      ! Set parameters for the fine grid based on lambda

      call gridparams(lambda,p,npt,npo,nt,no)

      ! Get fine composite Chebyshev discretization of K(tau,omega)

      allocate(kmat(nt,no),t(nt),om(no))

      call kfine_cc(lambda,eps,fb,npt,npo,p,t,om,kmat,kerr)

      write(6,*) ''
      write(6,*) '-------------- Fine K discretization --------------'
      write(6,*) ''
      write(6,*) '# fine grid pts in tau     = ',nt
      write(6,*) '# fine grid pts in omega   = ',no
      write(6,*) 'Max rel L^inf err in tau   = ',kerr(1)
      write(6,*) 'Max rel L^inf err in omega = ',kerr(2)


      ! Select real frequency points for DLR basis

      rank = 500 ! Upper bound on possible rank

      allocate(dlrrf(rank),oidx(rank))

      call dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx)


      ! Get DLR Matsubara frequency grid

      allocate(dlrmf(rank))

      call dlr_mf(rank,dlrrf,nmax,fb,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients transform matrix in LU form

      allocate(mf2cf(rank,rank),ipiv(rank))

      call dlr_mf2cf(nmax,rank,dlrrf,dlrmf,fb,mf2cf,ipiv)


      ! --- Compute actual eps-rank of fine grid K matrix by SVD ---

      call epsrank(nt,no,kmat,eps,i)

      write(6,*) ''
      write(6,*) '-------------------- DLR basis --------------------'
      write(6,*) ''
      write(6,*) 'DLR rank                          = ',rank
      write(6,*) 'eps-rank of fine K discretization = ',i


      ! --- Sample Green's function and get DLR ---

      ! Initialize Green's function evaluator

      pg = 24
      npg = npo
      
      allocate(xgl(pg),wgl(pg),xgj(pg),wgj(pg),pbpg(2*no+1))

      call gfunsc_init(pg,npg,pbpg,xgl,wgl,xgj,wgj)


      ! Sample G(tau) at DLR grid points

      allocate(gdlr(rank))

      do i=1,rank

        call gfunsc_mf(pg,npg,pbpg,xgl,wgl,xgj,wgj,fb,beta,&
          dlrmf(i),gdlr(i))

      enddo


      ! Compute coefficients of DLR expansion from samples

      call dlr_mfexpnd(rank,mf2cf,ipiv,gdlr)


      ! --- Compare DLR with true Green's function ---

      allocate(ttst(ntst))

      ! Get test points at which to measure error of Green's function;
      ! note that values of t larger than 0.5 are computed as t-1 to
      ! full relative precision to avoid losing digits

      call testpts(ntst,ttst)

      errlinf = 0*one
      errl2 = 0*one
      gmax = 0*one
      gl2 = 0*one

      do i=1,ntst

        ! Evaluate Green's function

        call gfunsc(pg,npg,pbpg,xgl,wgl,xgj,wgj,fb,beta,ttst(i),gtrue)

        ! Evaluate DLR

        call dlr_eval(rank,fb,dlrrf,real(gdlr),ttst(i),gtest)

        ! Update L^inf and L^2 errors, norms

        errlinf = max(errlinf,abs(gtrue-gtest))
        errl2 = errl2 + (gtrue-gtest)**2

        gmax = max(gmax,abs(gtrue))
        gl2 = gl2 + gtrue**2

      enddo

      errl2 = sqrt((ttst(2)-ttst(1))*errl2)
      gl2 = sqrt((ttst(2)-ttst(1))*gl2)

      write(6,*) ''
      write(6,*) '-------------------- DLR error --------------------'
      write(6,*) ''
      write(6,*) 'Abs L^inf err = ',errlinf
      write(6,*) 'Abs L^2 err   = ',errl2
      write(6,*) 'Rel L^inf err = ',errlinf/gmax
      write(6,*) 'Rel L^2 err   = ',errl2/gl2
      write(6,*) ''


      end subroutine dlr_sc_mf_test_main



      subroutine gfunsc_init(n,np,pbp,xgl,wgl,xgj,wgj)

      ! --- Initialization routine for evaluation of Green's function
      ! with semi-circular density ---

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n)

      integer i
      real *8 one

      one = 1.0d0

      ! --- Gauss-Legendre and Gauss-Jacobi quadrature ---

      call cdgqf(n,1,0.0d0,0.0d0,xgl,wgl)
      call cdgqf(n,4,0.5d0,0.0d0,xgj,wgj)

      ! --- Panels endpoints for composite quadrature rule ---

      pbp(np+1) = 0*one
      do i=1,np
        pbp(np+i+1) = one/2**(np-i)
      enddo
      pbp(1:np) = -pbp(2*np+1:np+2:-1)


      end subroutine gfunsc_init


      subroutine gfunsc(n,np,pbp,xgl,wgl,xgj,wgj,fb,beta,t,val)

      ! Evaluate Green's function with semi-circular density
      !
      ! This is a wrapper for the main subroutine, gfunsc1

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta,t,val
      real *8, external :: kfunf,kfunb
      character :: fb

      if (fb.eq.'f') then
        call gfunsc1(n,np,pbp,xgl,wgl,xgj,wgj,kfunf,beta,t,val)
      elseif (fb.eq.'b') then
        call gfunsc1(n,np,pbp,xgl,wgl,xgj,wgj,kfunb,beta,t,val)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine gfunsc


      subroutine gfunsc1(n,np,pbp,xgl,wgl,xgj,wgj,kfun,beta,t,val)

      ! Evaluate Green's function with semi-circular density
      !
      ! Main subroutine

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta,t,val
      real *8, external :: kfun

      integer ii,jj
      real *8 one,a,b,x,tt

      one = 1.0d0

      ! Treat t near 1 by symmetry to maintain high relative precision
      ! in the value of t. Note t near 1 is store by the negative of
      ! its distance to 1.

      tt = abs(t)

      val = 0.0d0
      do ii=2,2*np-1
        a = pbp(ii)
        b = pbp(ii+1)
        do jj=1,n
          x = a+(b-a)*(xgl(jj)+one)/2
          val = val + (b-a)/2*wgl(jj)*kfun(tt,beta*x)*&
            sqrt(one-x**2)
        enddo
      enddo

      a = one/2
      b = one
      do jj=1,n
        x = a+(b-a)*(xgj(jj)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(jj)*&
          kfun(tt,beta*x)*sqrt(one+x)
      enddo

      a = -one
      b = -one/2
      do jj=1,n
        x = a+(b-a)*(-xgj(n-jj+1)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(n-jj+1)*&
          kfun(tt,beta*x)*sqrt(one-x)
      enddo
        

      end subroutine gfunsc1


      subroutine gfunsc_mf(n,np,pbp,xgl,wgl,xgj,wgj,fb,beta,m,val)

      ! Evaluate Green's function with semi-circular density in
      ! Matsubara frequency domain
      !
      ! This is a wrapper for the main subroutine, gfunsc1

      implicit none
      integer n,np,m
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta
      complex *16 val
      complex *16, external :: kfunf_mf
      character :: fb

      if (fb.eq.'f') then
        call gfunsc_mf1(n,np,pbp,xgl,wgl,xgj,wgj,kfunf_mf,beta,m,val)
      !elseif (fb.eq.'b') then
      !  call gfunsc_mf1(n,np,pbp,xgl,wgl,xgj,wgj,kfunb,beta,t,val)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine gfunsc_mf


      subroutine gfunsc_mf1(n,np,pbp,xgl,wgl,xgj,wgj,kfun,beta,m,val)

      ! Evaluate Green's function with semi-circular density in
      ! Matsubara frequency domain
      !
      ! Main subroutine

      implicit none
      integer n,np,m
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta
      complex *16 val
      complex *16, external :: kfun

      integer ii,jj
      real *8 one,a,b,x

      one = 1.0d0

      val = 0.0d0
      do ii=2,2*np-1
        a = pbp(ii)
        b = pbp(ii+1)
        do jj=1,n
          x = a+(b-a)*(xgl(jj)+one)/2
          val = val + (b-a)/2*wgl(jj)*kfun(m,beta*x)*&
            sqrt(one-x**2)
        enddo
      enddo

      a = one/2
      b = one
      do jj=1,n
        x = a+(b-a)*(xgj(jj)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(jj)*&
          kfun(m,beta*x)*sqrt(one+x)
      enddo

      a = -one
      b = -one/2
      do jj=1,n
        x = a+(b-a)*(-xgj(n-jj+1)+one)/2
        val = val + ((b-a)/2)**(1.5d0)*wgj(n-jj+1)*&
          kfun(m,beta*x)*sqrt(one-x)
      enddo
        

      end subroutine gfunsc_mf1


      subroutine testpts(ntst,ttst)

      ! Get equispaced points on [0,1] to test accuracy of Green's
      ! function

      ! To maintain high precision for extremely large lambda and small
      ! eps, all points larger than 0.5 are computed and stored as their negative
      ! distance from 1; t > 0.5 is stored as t* = t-1, but computed directly
      ! as t*.

      implicit none
      integer ntst
      real *8 ttst(ntst)

      integer i

      do i=1,ntst

        if (i.le.ntst/2) then
          ttst(i) = (i-1)*1.0d0/(ntst-1)
        else
          ttst(i) = -(ntst-i)*1.0d0/(ntst-1)
        endif

      enddo

      end subroutine testpts



      subroutine epsrank(m,n,amat,eps,r)

      ! --- Get relative epsilon-rank of a matrix by SVD ---

      implicit none
      integer m,n,r
      real *8 amat(m,n),eps

      integer lwork,info,i
      real *8, allocatable :: s(:),work(:),tmp(:,:)

      lwork = 10*max(m,n)
      allocate(s(min(m,n)),work(lwork),tmp(m,n))

      tmp = amat

      call dgesvd('n','n',m,n,tmp,m,s,1,1,1,1,work,lwork,info)

      do r=0,min(m,n)-1
        if (s(r+1).lt.s(1)*eps) then
          return
        endif
      enddo

      end subroutine epsrank






!      ! --- Sum of delta functions density ---
!
!      real *8 function gfun(t,beta,fb)
!
!      implicit none
!      real *8 t,beta
!      character fb
!
!      integer ier,maxrec,numint
!      real *8 par1,par2,rint
!
!      gfun = (exp(-t*beta) + exp(-(1-t)*beta))/(1+exp(-beta))/beta
!
!
!      end function gfun


!      ! --- Single delta function density ---
!
!      real *8 function gfun(t,beta,fb)
!
!      implicit none
!      real *8 t,beta
!      character fb
!
!      real *8 tt
!
!      ! t near 1 is store by the negative of its distance to 1.
!
!      if (t.lt.0.0d0) then
!        tt = 1+t
!      else
!        tt = t
!      endif
!
!      gfun = (exp(-tt*beta)/(1+exp(-beta)))/beta
!
!      end function gfun
