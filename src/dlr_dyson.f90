      ! -------------------------------------------------------------
      !
      ! This file contains subroutines for solving the Dyson equation
      ! using the DLR
      !
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



      !> Solve the Dyson equation in imaginary time.
      !!
      !! Given a fixed self-energy, this routine forms the linear Dyson
      !! equation in imaginary time and solves it using Gaussian
      !! elimination.
      !!
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   it2cf   imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU factors
      !! @param[in]   it2cfp  imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU pivots
      !! @param[in]   phi     tensor taking DLR coefficients of g to
      !!                        matrix of convolution by g.
      !! @param[in]   g0      values of the right hand side G0 on
      !!                        the imaginary time grid
      !! @param[in]   g0mat   matrix of convolution by G0
      !! @param[in]   sig     values of the self-energy on the
      !!                          imaginary time grid
      !! @param[out]  g       solution of the linear Dyson equation
      !!                          on imaginary time grid

      subroutine dyson_it(r,it2cf,it2cfp,phi,g0,g0mat,sig,g)

      implicit none
      integer r,it2cfp(r)
      real *8 it2cf(r,r)
      real *8 g0(r),g0mat(r,r),sig(r),g(r)
      real *8 phi(r*r,r)

      integer j,info1
      integer, allocatable :: ipiv(:)
      real *8 one
      real *8, allocatable :: sigc(:),sigmat(:,:),sysmat(:,:)

      one = 1.0d0

      allocate(sigc(r),sigmat(r,r),sysmat(r,r),ipiv(r))

      ! Get matrix of convolution by self-energy

      call dlr_convmat(r,it2cf,it2cfp,phi,sig,sigmat)

      ! Form system matrix for linear Dyson equation

      call dgemm('N','N',r,r,r,-one,g0mat,r,sigmat,r,0*one,sysmat,r)

      do j=1,r
        sysmat(j,j) = one + sysmat(j,j)
      enddo

      ! Solve linear equation by LU factorization + backsolve

      call dgetrf(r,r,sysmat,r,ipiv,info1)

      g = g0

      call dgetrs('N',r,1,sysmat,r,ipiv,g,r,info1)


      end subroutine dyson_it





      !> Solve the Dyson equation in Matsubara frequency.
      !!
      !! Given a fixed self-energy, this routine solves the Dyson
      !! equation in the Matsubara frequency domain by diagonal
      !! inversion.
      !!
      !! @param[in]   beta    inverse temperature
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   g0      values of the right hand side G0 on
      !!                        the Matsubara frequency grid
      !! @param[in]   sigmf   values of the self-energy on the
      !!                        Matsubara frequency grid
      !! @param[out]  gmf     solution of the Dyson equation on the
      !!                        Matsubara frequency grid

      subroutine dyson_mf(beta,r,g0,sigmf,gmf)

      implicit none
      integer r
      real *8 beta
      complex *16 g0(r),sigmf(r),gmf(r)

      gmf = g0/(1.0d0-beta**2*g0*sigmf)

      end subroutine dyson_mf
