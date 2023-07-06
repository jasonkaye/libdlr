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
      
      
      program conv_exp
        
      ! Test convolution of two DLR expansions for Green's functions
      ! which are each a single exponential. Compare result with
      ! analytically-known convolution. Test two methods of
      ! convolution; directly forming matrix of convolution, and fast
      ! convolution in the DLR coefficient domain.
      
      implicit none
      integer ntst
      real *8 lambda,eps,beta

      ! Input parameters

      lambda = 1000.0d0   ! DLR high energy cutoff
      eps = 1.0d-14       ! DLR error tolerance

      ntst = 10000        ! # imaginary time test points
      beta = 1000.0d0     ! Inverse temperature


      ! Main test subroutine

      call conv_exp_main(lambda,eps,ntst,beta)


      end program conv_exp


      subroutine conv_exp_main(lambda,eps,ntst,beta)

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,j,r
      integer, allocatable :: it2cfp(:)
      real *8 one,gtrue,gtst1,gtst2
      real *8 err1l2,err2l2,err1linf,err2linf,gmax,gl2
      real *8, allocatable :: dlrit(:),dlrrf(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),fstconv(:,:)
      real *8, allocatable :: ttst(:),g1(:),g2(:),g31(:),g32(:)

      real *8, allocatable :: phi(:,:),gmat(:,:)

      one = 1.0d0


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get DLR coefficients -> imaginary time values matrix

      allocate(cf2it(r,r))

      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Sample G1 and G2 at imaginary time nodes

      allocate(g1(r),g2(r),g31(r),g32(r))

      do i=1,r

        call gfun1(beta,dlrit(i),g1(i))
        call gfun2(beta,dlrit(i),g2(i))

      enddo


      ! Get convolution tensor

      allocate(phi(r*r,r))

      call dlr_convtens(beta,-1,r,dlrrf,dlrit,it2cf,it2cfp,phi)


      ! Get matrix of convolution by G1

      allocate(gmat(r,r))

      call dlr_convmat(r,1,it2cf,it2cfp,phi,g1,gmat)


      ! Get convolution G3 of G1 with G2

      call dlr_conv(r,1,gmat,g2,g31)

      
      ! Get DLR coefficients of G3

      call dlr_it2cf(r,1,it2cf,it2cfp,g31,g31)



      ! Initialize fast convolution routine

      allocate(fstconv(r,2*r))

      call dlr_fstconv_init(beta,r,dlrrf,dlrit,cf2it,fstconv)


      ! Get convolution G3 of G1 with G2 by fast convolution

      call dlr_fstconv(r,1,cf2it,it2cf,it2cfp,fstconv,g1,g2,g32)


      ! Get DLR coefficients of G3

      call dlr_it2cf(r,1,it2cf,it2cfp,g32,g32)


      ! Get test points in relative format

      allocate(ttst(ntst))

      call eqpts_rel(ntst,ttst)


      ! Measure L^inf and L^2 errors of convolution

      err1linf = 0*one
      err1l2 = 0*one
      err2linf = 0*one
      err2l2 = 0*one
      gmax = 0*one
      gl2 = 0*one

      do i=1,ntst

        ! Evaluate true convolution

        call gfun3(beta,ttst(i),gtrue)

        ! Evaluate DLR

        call dlr_it_eval(r,1,dlrrf,g31,ttst(i),gtst1)
        call dlr_it_eval(r,1,dlrrf,g32,ttst(i),gtst2)

        ! Update L^inf and L^2 errors, norms

        err1linf = max(err1linf,abs(gtrue-gtst1))
        err1l2 = err1l2 + (gtrue-gtst1)**2
        err2linf = max(err2linf,abs(gtrue-gtst2))
        err2l2 = err2l2 + (gtrue-gtst2)**2

        gmax = max(gmax,abs(gtrue))
        gl2 = gl2 + gtrue**2

      enddo

      err1l2 = sqrt((ttst(2)-ttst(1))*err1l2)
      err2l2 = sqrt((ttst(2)-ttst(1))*err2l2)
      gl2 = sqrt((ttst(2)-ttst(1))*gl2)

      write(6,*) ''
      write(6,*) 'DLR rank = ',r
      write(6,*) 'Errors for first method -- convolution tensor:'
      write(6,*) 'Abs L^inf err = ',err1linf
      write(6,*) 'Abs L^2 err   = ',err1l2
      write(6,*) 'Rel L^inf err = ',err1linf/gmax
      write(6,*) 'Rel L^2 err   = ',err1l2/gl2

      write(6,*) 'Errors for second method -- fast convolution:'
      write(6,*) 'Abs L^inf err = ',err2linf
      write(6,*) 'Abs L^2 err   = ',err2l2
      write(6,*) 'Rel L^inf err = ',err2linf/gmax
      write(6,*) 'Rel L^2 err   = ',err2l2/gl2
      write(6,*) ''


      ! Return failed status if error is not sufficiently small

      if (err1linf.lt.1.0d-13.and.err2linf.lt.1.0d-13) then
        call exit(0)
      else
        call exit(1)
      endif

      end subroutine conv_exp_main



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

      subroutine gfun3(beta,t,g)

      ! Evaluate Green's function corresponding to convolution of two 
      ! single exponentials 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a1,a2

      a1 = 0.804d0
      a2 = -0.443d0

      g = (kfunf_rel(t,beta*a1)-kfunf_rel(t,beta*a2))/(a1-a2)

      end subroutine gfun3
