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
      
           
      program ha_mf

      ! Test discrete Lehmann representation using Green's function
      ! generated from Lehmann representation with density which is a
      ! sum of two delta functions. Recover DLR coefficients from
      ! values of Green's function on fermionic Matsubara frequency
      ! grid, and measure the error of the resulting DLR expansion on a
      ! test grid in imaginary time.

      implicit none
      integer ntst,nmax
      real *8 lambda,eps,beta

      ! Input parameters

      lambda = 1000.0d0   ! DLR high energy cutoff
      eps = 1.0d-14       ! DLR error tolerance
      nmax = 1000         ! DLR Matsubara frequency cutoff

      ntst = 10000        ! # imaginary time test points
      beta = 1000.0d0     ! Inverse temperature


      ! Main test subroutine

      call ha_mf_main(lambda,eps,nmax,ntst,beta)

      end program ha_mf


      subroutine ha_mf_main(lambda,eps,nmax,ntst,beta)

      implicit none
      integer ntst,nmax
      real *8 lambda,eps,beta

      integer i,j,r
      integer, allocatable :: dlrmf(:),mf2cfp(:)
      real *8 one,gtrue,errl2,errlinf,gmax,gl2,gtest
      real *8, allocatable :: ttst(:),dlrrf(:)
      real *8, allocatable :: xgl(:),wgl(:),xgj(:),wgj(:),pbpg(:),gc(:)
      complex *16, allocatable :: g(:),mf2cf(:,:)

      one = 1.0d0


      ! Get DLR frequencies, Matsubara frequency grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrmf(r))

      call dlr_mf_build(lambda,eps,nmax,-1,r,dlrrf,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients matrix (LU
      ! form)

      allocate(mf2cf(r,r),mf2cfp(r))

      call dlr_mf2cf(nmax,r,dlrrf,dlrmf,-1,mf2cf,mf2cfp)


      ! Sample G at Matsubara frequency nodes

      allocate(g(r))

      do i=1,r

        call gfun_mf(beta,dlrmf(i),g(i))

      enddo


      ! Get DLR coefficients of G

      allocate(gc(r))

      call dlr_mf_expnd(r,mf2cf,mf2cfp,g,gc)


      ! Get test points in relative format

      allocate(ttst(ntst))

      call eqpts_rel(ntst,ttst)


      ! Measure L^inf and L^2 errors

      errlinf = 0*one
      errl2 = 0*one
      gmax = 0*one
      gl2 = 0*one

      do i=1,ntst

        ! Evaluate Green's function

        call gfun_it(beta,ttst(i),gtrue)

        ! Evaluate DLR

        call dlr_it_eval(r,dlrrf,gc,ttst(i),gtest)

        ! Update L^inf and L^2 errors, norms

        errlinf = max(errlinf,abs(gtrue-gtest))
        errl2 = errl2 + (gtrue-gtest)**2

        gmax = max(gmax,abs(gtrue))
        gl2 = gl2 + gtrue**2

      enddo

      errl2 = sqrt((ttst(2)-ttst(1))*errl2)
      gl2 = sqrt((ttst(2)-ttst(1))*gl2)


      write(6,*) ''
      write(6,*) 'DLR rank = ',r
      write(6,*) 'Abs L^inf err = ',errlinf
      write(6,*) 'Abs L^2 err   = ',errl2
      write(6,*) 'Rel L^inf err = ',errlinf/gmax
      write(6,*) 'Rel L^2 err   = ',errl2/gl2
      write(6,*) ''


      ! Return failed status if error is not sufficiently small

      if (errlinf.gt.1.0d-12) then
        call exit(1)
      endif


      end subroutine ha_mf_main


      subroutine gfun_mf(beta,n,g)

      ! Evaluate Green's function with sum-of-delta-functions spectral
      ! density in fermionic Matsubara frequency domain

      implicit none
      integer n
      real *8 beta
      complex *16 g
      complex *16, external :: kfunmf

      real *8 a1,a2,a3,a4,a5

      a1 = -0.804d0
      a2 = -0.443d0
      a3 =  0.093d0
      a4 =  0.915d0
      a5 =  0.929d0

      g = kfunmf(2*n+1,beta*a1) + kfunmf(2*n+1,beta*a2) +&
        kfunmf(2*n+1,beta*a3) + kfunmf(2*n+1,beta*a4) +&
        kfunmf(2*n+1,beta*a5)
        
      end subroutine gfun_mf



      subroutine gfun_it(beta,t,g)

      ! Evaluate Green's function with sum-of-delta-functions spectral
      ! density in imaginary time domain

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a1,a2,a3,a4,a5

      a1 = -0.804d0
      a2 = -0.443d0
      a3 =  0.093d0
      a4 =  0.915d0
      a5 =  0.929d0

      g = kfunf_rel(t,beta*a1) + kfunf_rel(t,beta*a2) &
        + kfunf_rel(t,beta*a3) + kfunf_rel(t,beta*a4) &
        + kfunf_rel(t,beta*a5)

      end subroutine gfun_it
