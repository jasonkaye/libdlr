      ! -------------------------------------------------------------
      !
      ! Copyright (C) 2021 The Simons Foundation
      ! 
      ! Author: Jason Kaye
      ! 
      ! -------------------------------------------------------------
      ! 
      ! libdlr is licensed under the Apache License, Version 2.0 (the
      ! "License"); you may not use this file except in compliance with
      ! the License.  You may obtain a copy of the License at
      ! 
      !     http://www.apache.org/licenses/LICENSE-2.0
      ! 
      ! Unless required by applicable law or agreed to in writing,
      ! software distributed under the License is distributed on an "AS
      ! IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
      ! express or implied.  See the License for the specific language
      ! governing permissions and limitations under the License.
      ! 
      ! -------------------------------------------------------------
      
      
      program dyson_sc

      ! Solve Dyson equation with self-energy Sigma = c^2*G both in
      ! imaginary time and Matsubara frequency
      !
      ! Since we know the solution has a semi-circular spectral
      ! function, we can begin with the correct self-energy, and make
      ! sure we recover the correct Green's function
      !
      ! The self-energy and correct reference Green's function are
      ! computed from the Lehmann representation in imaginary time using
      ! high-order numerical integration
      
      implicit none
      integer ntst,nmax
      real *8 lambda,eps,beta,c

      ! Input parameters

      lambda = 200.0d0  ! DLR high energy cutoff
      eps = 1.0d-14     ! DLR error tolerance
      nmax = 200.0d0    ! DLR Matsubara frequency cutoff

      beta = 20.0d0     ! Inverse temperature
      c = 1.0d0/2       ! Self-energy strength

      ntst = 1000       ! # output points


      ! Main test subroutine

      call dyson_sc_main(lambda,eps,nmax,ntst,beta,c)

      end program dyson_sc


      subroutine dyson_sc_main(lambda,eps,nmax,ntst,beta,c)
        
      ! Main driver routine for Dyson solver test 

      implicit none
      integer ntst,nmax
      real *8 lambda,eps,beta,c

      integer i,r,npg,pg
      integer, allocatable :: dlrmf(:),it2cfp(:),mf2cfp(:)
      real *8 err1,err2
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:),g1(:)
      real *8, allocatable :: phi(:,:),g0it(:),g1c(:),g2c(:)
      real *8, allocatable :: xgl(:),wgl(:),xgj(:),wgj(:),pbpg(:)
      real *8, allocatable :: sig(:),sigc(:),g0mat(:,:)
      real *8, allocatable :: it_tst(:),g1tst(:),g2tst(:),gtrue(:)
      complex *16, allocatable :: mf2cf(:,:),cf2mf(:,:),g0mf(:),sigmf(:)
      complex *16, allocatable :: gmf(:)
      real *8, external :: kfunf_rel


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Get Matsubara frequency grid

      allocate(dlrmf(r))

      call dlr_mf(nmax,r,dlrrf,-1,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients matrix (LU
      ! form)

      allocate(mf2cf(r,r),mf2cfp(r),cf2mf(r,r))

      call dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,-1,mf2cf,mf2cfp)


      ! Get DLR coefficients -> Matsubara frequency values matrix

      call dlr_cf2mf_init(r,dlrrf,dlrmf,-1,cf2mf)



      ! Get convolution tensor

      allocate(phi(r*r,r))

      call dlr_convtens(beta,-1,r,dlrrf,dlrit,it2cf,it2cfp,phi)


      ! Get free particle Green's function on imaginary time grid

      allocate(g0it(r))

      call getg0_it(beta,r,dlrit,0.0d0,g0it)


      ! Get free particle Green's function on Matsubara frequency grid

      allocate(g0mf(r))

      call getg0_mf(beta,r,dlrmf,0.0d0,g0mf)


      ! Get matrix of convolution by free particle Green's function

      allocate(g0mat(r,r))

      call dlr_convmat(r,1,it2cf,it2cfp,phi,g0it,g0mat)



      ! Initialize evaluator for Green's function with semi-circular
      ! spectral density

      pg = 24
      npg = max(ceiling(log(lambda)/log(2.0d0)),1)
      
      allocate(xgl(pg),wgl(pg),xgj(pg),wgj(pg),pbpg(2*npg+1))

      call gfun_init(pg,npg,pbpg,xgl,wgl,xgj,wgj)


      ! Get self-energy at imaginary time nodes

      allocate(sig(r))

      do i=1,r

        call gfun_it(pg,npg,pbpg,xgl,wgl,xgj,wgj,beta,dlrit(i),&
          sig(i))

        sig(i) = c*c*sig(i)

      enddo


      ! Solve Dyson equation in imaginary time
      
      allocate(g1(r),g1c(r))

      call dyson_it(r,1,it2cf,it2cfp,phi,g0it,g0mat,sig,g1)


      ! Get DLR coefficients of solution

      call dlr_it2cf(r,it2cf,it2cfp,g1,g1c)




      ! Transform self-energy to Matsubara frequency and solve Dyson
      ! equation

      allocate(sigc(r),sigmf(r),gmf(r),g2c(r))

      ! Get DLR coefficients of self-energy

      call dlr_it2cf(r,it2cf,it2cfp,sig,sigc)


      ! Get self-energy on Matsubara frequency grid

      call dlr_cf2mf(r,cf2mf,sigc,sigmf)


      ! Solve Dyson equation by diagonal inversion

      call dyson_mf(beta,r,g0mf,sigmf,gmf)


      ! Get DLR coefficients of solution

      call dlr_mf2cf(r,mf2cf,mf2cfp,gmf,g2c)



      ! Get test points in relative format

      allocate(it_tst(ntst))

      call eqpts_rel(ntst,it_tst)


      ! Evaluate solutions on test grid

      allocate(g1tst(ntst),g2tst(ntst),gtrue(ntst))

      do i=1,ntst

        call dlr_it_eval(r,dlrrf,g1c,it_tst(i),g1tst(i))
        
        call dlr_it_eval(r,dlrrf,g2c,it_tst(i),g2tst(i))

        call gfun_it(pg,npg,pbpg,xgl,wgl,xgj,wgj,beta,it_tst(i),&
          gtrue(i))

      enddo


      ! Compute L^inf errors

      err1 = maxval(abs(g1tst-gtrue))
      err2 = maxval(abs(g2tst-gtrue))

      write(6,*) ''
      write(6,*) '-------------------- DLR error --------------------'
      write(6,*) ''
      write(6,*) 'Imag time L^inf err = ',err1
      write(6,*) 'Mats freq L^inf err = ',err2
      write(6,*) ''


      ! Return failed status if error is not sufficiently small

      if (err1.gt.1.0d-12.or.err2.gt.1.0d-12) then
        call exit(1)
      endif

      end subroutine dyson_sc_main


      subroutine getg0_it(beta,r,dlrit,mu,g0)

      ! Evaluate free particle Green's function in imaginary time

      implicit none
      integer r
      real *8 beta,dlrit(r),mu,g0(r)

      integer i
      real *8, external :: kfunf_rel

      do i=1,r
        g0(i) = -kfunf_rel(dlrit(i),beta*mu)
      enddo

      end subroutine getg0_it


      subroutine getg0_mf(beta,r,dlrmf,mu,g0)

      ! Evaluate free particle Green's function in Matsubara frequency

      implicit none
      integer r,dlrmf(r)
      real *8 beta,mu
      complex *16 g0(r)

      integer i
      complex *16, external :: kfunmf

      do i=1,r
        g0(i) = -kfunmf(2*dlrmf(i)+1,beta*mu)
      enddo

      end subroutine getg0_mf


      subroutine gfun_init(n,np,pbp,xgl,wgl,xgj,wgj)

      ! Initialize subroutine gfun_it, which evaluates Green's function
      ! with semi-circular density at an imaginary time point

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n)

      integer i
      real *8 one

      one = 1.0d0

      ! Gauss-Legendre and Gauss-Jacobi quadrature

      call cdgqf(n,1,0.0d0,0.0d0,xgl,wgl)
      call cdgqf(n,4,0.5d0,0.0d0,xgj,wgj)

      ! Panels endpoints for composite quadrature rule

      pbp(np+1) = 0*one
      do i=1,np
        pbp(np+i+1) = one/2**(np-i)
      enddo
      pbp(1:np) = -pbp(2*np+1:np+2:-1)


      end subroutine gfun_init


      subroutine gfun_it(n,np,pbp,xgl,wgl,xgj,wgj,beta,t,val)

      ! Evaluate Green's function with semi-circular density
      !
      ! Uses composite Gauss-Legendre quadrature with dyadically-refined
      ! panels, Gauss-Jacobi quadrature at endpoints

      implicit none
      integer n,np
      real *8 pbp(2*np+1),xgl(n),wgl(n),xgj(n),wgj(n),beta,t,val

      integer ii,jj
      real *8 one,pi,a,b,x,tt
      real *8, external :: kfunf

      one = 1.0d0
      pi = 4*atan(1.0d0)

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

      val = -2/pi*val

      end subroutine gfun_it
