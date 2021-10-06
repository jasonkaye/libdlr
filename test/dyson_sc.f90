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

      ! Solve nonlinear Dyson equation with self-energy Sigma = c^2*G
      ! self-consistently by weighted fixed point iteration using both
      ! pure imaginary time method and Matsubara frequency/imaginary
      ! time method
      !
      ! We solve the Dyson equation self-consistently by a weighted
      ! fixed point iteration, with weight w assigned to the new iterate
      ! and weight 1-w assigned to the previous iterate.
      !
      ! The solution has a semi-circular spectral function, and we
      ! measure the error against a reference solution computed from the
      ! Lehmann representation using high-order numerical integration
      
      implicit none
      integer ntst,maxit,nmax
      real *8 lambda,eps,fptol,w,beta,c

      ! Input parameters

      lambda = 200.0d0  ! DLR high energy cutoff
      eps = 1.0d-14     ! DLR error tolerance
      nmax = 200.0d0    ! DLR Matsubara frequency cutoff

      beta = 20.0d0     ! Inverse temperature
      c = 1.0d0/2       ! Self-energy strength

      maxit = 1000      ! Max # fixed point iterations
      fptol = 1.0d-12   ! Fixed point tolerance
      w = 0.8d0         ! Fixed point iteration weighting
      
      ntst = 1000       ! # output points


      ! Main test subroutine

      call dyson_sc_main(lambda,eps,nmax,ntst,fptol,maxit,w,&
        beta,c)

      end program dyson_sc


      subroutine dyson_sc_main(lambda,eps,nmax,ntst,fptol,&
          maxit,w,beta,c)
        
      ! Main driver routine for nonlinear Dyson solver test 

      implicit none
      integer ntst,maxit,nmax
      real *8 lambda,eps,fptol,w,beta,c

      integer i,j,r,info,numit,npg,npo,pg
      integer, allocatable :: dlrmf(:),it2cfp(:),mf2cfp(:)
      real *8 err1,err2
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:),g1(:),g2(:)
      real *8, allocatable :: cf2it(:,:),phi(:,:),g0it(:),g1c(:),g2c(:)
      real *8, allocatable :: xgl(:),wgl(:),xgj(:),wgj(:),pbpg(:)
      real *8, allocatable :: it_tst(:),g1tst(:),g2tst(:),gtrue(:)
      complex *16, allocatable :: mf2cf(:,:),cf2mf(:,:),g0mf(:)
      real *8, external :: kfunf_rel


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_buildit(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Get DLR coefficients -> imaginary time values matrix

      allocate(cf2it(r,r))

      call dlr_cf2it(r,dlrrf,dlrit,cf2it)


      ! Get Matsubara frequency grid

      allocate(dlrmf(r))

      call dlr_mf(nmax,r,dlrrf,-1,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients matrix (LU
      ! form)

      allocate(mf2cf(r,r),mf2cfp(r),cf2mf(r,r))

      call dlr_mf2cf(nmax,r,dlrrf,dlrmf,-1,mf2cf,mf2cfp)


      ! Get DLR coefficients -> Matsubara frequency values matrix

      call dlr_cf2mf(r,dlrrf,dlrmf,-1,cf2mf)



      ! Get convolution tensor

      allocate(phi(r*r,r))

      call dlr_convtens(beta,-1,r,dlrrf,dlrit,it2cf,it2cfp,phi)


      ! Get free particle Green's function on imaginary time grid

      allocate(g0it(r))

      call getg0_it(beta,r,dlrit,0.0d0,g0it)


      ! Get free particle Green's function on Matsubara frequency grid

      allocate(g0mf(r))

      call getg0_mf(beta,r,dlrmf,0.0d0,g0mf)


      ! Solve nonlinear Dyson equation by imaginary time method

      allocate(g1(r))

      numit = maxit

      g1 = g0it

      call dlr_dyson_it(beta,r,dlrit,it2cf,it2cfp,cf2it,phi,&
        sigfun,w,fptol,numit,g0it,g1,info)


      ! Solve nonlinear Dyson equation by Matsubara frequency method

      allocate(g2(r))

      numit = maxit

      g2 = g0it

      call dlr_dyson_mf(beta,r,dlrit,it2cf,it2cfp,cf2it,&
        dlrmf,mf2cf,mf2cfp,cf2mf,sigfun,w,fptol,numit,g0mf,g2,info)


      ! Get DLR coefficients of solutions

      allocate(g1c(r),g2c(r))

      call dlr_it_expnd(r,it2cf,it2cfp,g1,g1c)
      call dlr_it_expnd(r,it2cf,it2cfp,g2,g2c)


      ! Initialize Green's function evaluator (semi-circular spectral
      ! density)

      pg = 24
      npg = max(ceiling(log(lambda)/log(2.0d0)),1)
      
      allocate(xgl(pg),wgl(pg),xgj(pg),wgj(pg),pbpg(2*npg+1))

      call gfun_init(pg,npg,pbpg,xgl,wgl,xgj,wgj)


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


      contains

        subroutine sigfun(r,g,sig)

        ! Self-energy evaluator

        implicit none
        integer r
        real *8 g(r),sig(r)

        sig = c**2*g

        end subroutine sigfun
      
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

