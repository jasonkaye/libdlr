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
      
      
      program ip_exp
        
      ! Test L2 inner product of two DLR expansions for Green's
      ! functions which are each a single exponential. Compare result
      ! with analytically-known inner product.
      
      implicit none
      real *8 lambda,eps,beta

      ! Input parameters

      lambda = 10.0d0   ! DLR high energy cutoff
      eps = 1.0d-14       ! DLR error tolerance

      beta = 10.0d0     ! Inverse temperature


      ! Main test subroutine

      call ip_exp_main(lambda,eps,beta)


      end program ip_exp


      subroutine ip_exp_main(lambda,eps,beta)

      implicit none
      real *8 lambda,eps,beta

      integer j,k,r
      integer, allocatable :: it2cfp(:)
      real *8 ip,err
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:),g1(:),g2(:)
      real *8, allocatable :: ipmat(:,:)
      real *8, external :: iptrue


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Sample G1 and G2 at imaginary time nodes

      allocate(g1(r),g2(r))

      do j=1,r

        call gfun1(beta,dlrit(j),g1(j))
        call gfun2(beta,dlrit(j),g2(j))

      enddo


      ! Get inner product weight matrix ipmat

      allocate(ipmat(r,r))

      call dlr_ipmat(beta,r,dlrit,dlrrf,it2cf,it2cfp,ipmat)

      
      ! Compute inner product by g1^T * ipmat * g2

      g2 = matmul(ipmat,g2)
      ip = sum(g1*g2)


      ! Compute error

      err = abs(ip-iptrue(beta))

      write(6,*) ''
      write(6,*) 'DLR rank = ',r
      write(6,*) 'Abs err = ',err
      write(6,*) 'Rel err   = ',err/abs(iptrue(beta))
      write(6,*) ''


      ! Return failed status if error is not sufficiently small

      if (err.gt.1.0d-13) then
        call exit(1)
      endif

      end subroutine ip_exp_main



      subroutine gfun1(beta,t,g)

      ! Evaluate single exponential Green's function

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a1

      a1 = 0.804d0

      g = kfunf_rel(t,beta*a1)

      end subroutine gfun1

      subroutine gfun2(beta,t,g)

      ! Evaluate single exponential Green's function

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a2

      a2 = -0.443d0

      g = kfunf_rel(t,beta*a2)

      end subroutine gfun2

      function iptrue(beta)

      ! Evaluate Green's function corresponding to convolution of two 
      ! single exponentials 

      implicit none
      real *8 iptrue,beta

      real *8 a1,a2
      real *8, external :: kfunf

      a1 = 0.804d0
      a2 = -0.443d0

      a1 = beta*a1
      a2 = beta*a2

      iptrue = beta*kfunf(0.0d0,a1)*kfunf(0.0d0,a2)&
        *(1.0d0-exp(-(a1+a2)))/(a1+a2)

      end function iptrue
