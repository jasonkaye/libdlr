      ! -------------------------------------------------------------
      !
      ! Copyright (C) 2025 The Simons Foundation
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
      
           
      program timetofreqtotime

      ! For a single exponential Green's function, form the DLR expansion from
      ! samples in imaginary time, transform to Matsubara frequency, and then
      ! transform back, measuring errors at each step on a dense grid
        
      implicit none
      integer ntst_it,ntst_mf,nmax
      real *8 lambda,eps,beta

      ! Input parameters

      lambda = 100.0d0  ! DLR high energy cutoff
      eps = 1.0d-10     ! DLR error tolerance
      nmax = 100        ! DLR Matsubara frequency cutoff

      beta = 100.0d0    ! Inverse temperature

      ntst_it = 1000    ! # imaginary time test points
      ntst_mf = 1000    ! Matsubara frequency test point cutoff

      call timetofreqtotime_main(lambda,eps,nmax,ntst_it,ntst_mf,beta)

      end program timetofreqtotime


      subroutine timetofreqtotime_main(lambda,eps,nmax,ntst_it,ntst_mf,beta)

      implicit none
      integer ntst_it,ntst_mf,nmax
      real *8 lambda,eps,beta

      integer i,r
      integer, allocatable :: it2cfp(:),mf_tst(:),dlrmf(:),mf2cfp(:)
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:),g(:),gc(:)
      real *8, allocatable :: it_tst(:),gtst_it(:),gtrue_it(:),gc2(:)
      complex *16, allocatable :: gtst_mf(:),gtrue_mf(:),mf2cf(:,:),gmf(:)

      ! Get DLR frequencies, imaginary time grid
      r = 500 ! Upper bound on DLR rank
      allocate(dlrrf(r),dlrit(r))
      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)

      ! Get imaginary time values -> DLR coefficients matrix (LU form)
      allocate(it2cf(r,r),it2cfp(r))
      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Initialize reference Green's function evaluator (note: this is
      ! not a DLR-related operation; rather, we are initializing a
      ! special-purpose subroutine to evaluate the Green's function with
      ! a semi-circular spectral density. The initialization and
      ! evaluation subroutines are defined at the bottom of this file.)

      ! Sample G at imaginary time nodes
      allocate(g(r))
      do i=1,r
        call gfun_it(beta,dlrit(i),g(i))
      enddo

      ! Get DLR coefficients of G
      allocate(gc(r))
      call dlr_it2cf(r,1,it2cf,it2cfp,g,gc)

      ! Get imaginary time test points in relative format
      allocate(it_tst(ntst_it))
      call eqpts_rel(ntst_it,it_tst)

      ! Evaluate DLR and true Green's function on imaginary time test grid
      allocate(gtst_it(ntst_it),gtrue_it(ntst_it))
      do i=1,ntst_it
        call dlr_it_eval(r,1,dlrrf,gc,it_tst(i),gtst_it(i))
        call gfun_it(beta,it_tst(i),gtrue_it(i))
      enddo

      ! Measure error
      write(6,*) 'lambda = ',lambda
      write(6,*) 'epsilon = ',eps
      write(6,*) 'beta = ',beta
      write(6,*) ''
      write(6,*) 'DLR rank = ',r
      write(6,*) ''
      write(6,*) 'Testing error of initial imaginary time DLR expansion...'
      write(6,*) 'Imag time L^inf error = ',maxval(abs(gtst_it-gtrue_it))
      write(6,*) ''

      
      ! Get Matsubara frequency test points
      allocate(mf_tst(2*ntst_mf+1))
      do i=1,2*ntst_mf+1
        mf_tst(i) = -ntst_mf+i-1
      enddo

      ! Evaluate DLR and true Green's function on Matsubara frequency test grid
      allocate(gtst_mf(2*ntst_mf+1),gtrue_mf(2*ntst_mf+1))
      i=1
      do i=1,2*ntst_mf+1
        call dlr_mf_eval(r,1,dlrrf(1:r),-1,gc,mf_tst(i),beta,gtst_mf(i))
        call gfun_mf(beta,mf_tst(i),gtrue_mf(i))
      enddo

      ! Measure error
      write(6,*) 'Testing error of Matsubara frequency DLR expansion...'
      write(6,*) 'Mats freq L^inf error = ',maxval(abs(gtst_mf-gtrue_mf))
      write(6,*) ''

      ! Get DLR Matsubara frequency nodes and values -> coefficients matrix
      allocate(dlrmf(r),mf2cf(r,r),mf2cfp(r))
      call dlr_mf(nmax,r,dlrrf,-1,dlrmf)
      call dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,-1,mf2cf,mf2cfp)

      ! Evaluate DLR expansion at Matsubara frequency nodes
      allocate(gmf(r))
      do i=1,r
        call dlr_mf_eval(r,1,dlrrf,-1,gc,dlrmf(i),beta,gmf(i))
      enddo

       ! Form DLR expansion from Matsubara frequency values
       allocate(gc2(r))
       call dlr_mf2cf(r,1,mf2cf,mf2cfp,beta,gmf,gc2)
 
       ! Evaluate DLR expansion on imaginary time test grid
       do i=1,ntst_it
         call dlr_it_eval(r,1,dlrrf,gc2,it_tst(i),gtst_it(i))
       enddo
       write(6,*) 'Testing error of imaginary time DLR expansion formed from DLR Matsubara grid values...'
       write(6,*) 'Imag time L^inf error = ',maxval(abs(gtst_it-gtrue_it))
       write(6,*) ''

      end subroutine timetofreqtotime_main


      subroutine gfun_it(beta,t,g)

      ! Evaluate Green's function with sum-of-delta-functions spectral
      ! density in imaginary time domain

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 om0

      om0 = -0.804d0
      g = kfunf_rel(t,beta*om0)


      end subroutine gfun_it




      subroutine gfun_mf(beta,n,g)

      ! Evaluate Green's function with sum-of-delta-functions spectral
      ! density in fermionic Matsubara frequency domain

      implicit none
      integer n
      real *8 beta
      complex *16 g
      complex *16, external :: kfunmf

      real *8 om0

      om0 = -0.804d0
      g = beta*(kfunmf(2*n+1,beta*om0))
        
      end subroutine gfun_mf