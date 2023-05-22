      ! -------------------------------------------------------------
      !
      ! Copyright (C) 2022 The Simons Foundation
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
      
      
      program dyson_mat

      ! Solve 2x2 matrix-valued Dyson equation in imaginary time
      !
      ! To obtain an exact solution for testing, we begin with a 3x3
      ! Hamiltonian h and zero self-energy. The corresponding Green's
      ! function G_ij can be obtained by diagonalization. It can then be
      ! shown that the upper left 2x2 block of G_ij can also be obtained
      ! by solving the Dyson equation with single-particle Hamiltonian
      ! given by h(1:2,1:2), and self-energy given by
      ! h(1:2,3)*G0*h(3,1:2), where G0 is the free-particle Green's
      ! function with energy h(3,3). By comparing the Green's function
      ! obtained by this approach with the known exact solution, we can
      ! test a matrix-valued Dyson solver.
      
      implicit none
      integer ntst
      real *8 lambda,eps,beta,h(3,3)

      ! Input parameters

      lambda = 200.0d0  ! DLR high energy cutoff
      eps = 1.0d-14     ! DLR error tolerance

      beta = 20.0d0     ! Inverse temperature

      h(1,1) = -3.0d0   ! 3x3 Hamiltonian matrix
      h(2,2) = -1.0d0
      h(3,3) = 1.0d0
      h(2,1) = 0.3d0
      h(1,2) = 0.3d0
      h(3,1) = -0.2d0
      h(1,3) = -0.2d0
      h(3,2) = 0.1d0
      h(2,3) = 0.1d0

      ntst = 101        ! # output points


      ! Main test subroutine

      call dyson_mat_main(lambda,eps,ntst,beta,h)

      end program dyson_mat


      subroutine dyson_mat_main(lambda,eps,ntst,beta,h)
        
      ! Main driver routine for matrix-valued Dyson solver test 

      implicit none
      integer ntst
      real *8 lambda,eps,beta,h(3,3)

      integer i,j,k,r
      integer, allocatable :: dlrmf(:),it2cfp(:)
      real *8 err
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:)
      real *8, allocatable :: phi(:,:),g0it(:,:,:),g0mat(:,:)
      real *8, allocatable :: g(:,:,:),gc(:,:,:),g02it(:),sig(:,:,:)
      real *8, allocatable :: gtrue(:,:,:),gctrue(:,:,:)
      real *8, allocatable :: it_tst(:),gtst(:,:,:),gtruetst(:,:,:)


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Get convolution tensor

      allocate(phi(r*r,r))

      call dlr_convtens(beta,-1,r,dlrrf,dlrit,it2cf,it2cfp,phi)



      ! Get free particle Green's function for upper-left 2x2 block of
      ! Hamiltonian

      allocate(g0it(r,2,2))

      call getg0_it(beta,r,2,dlrit,h(1:2,1:2),0.0d0,g0it)


      ! Get matrix of convolution by free particle Green's function

      allocate(g0mat(2*r,2*r))

      call dlr_convmat(r,2,it2cf,it2cfp,phi,g0it,g0mat)

     

      ! Get self-energy from free particle Green's function for
      ! lower-right entry of Hamiltonian

      allocate(g02it(r),sig(r,2,2))

      call getg0_it(beta,r,1,dlrit,h(3,3),0.0d0,g02it)

      sig(:,1,1) = g02it*h(1,3)*h(3,1)
      sig(:,2,1) = g02it*h(2,3)*h(3,1)
      sig(:,1,2) = g02it*h(1,3)*h(3,2)
      sig(:,2,2) = g02it*h(2,3)*h(3,2)


      ! Solve Dyson equation
      
      allocate(g(r,2,2),gc(r,2,2))

      call dyson_it(r,2,it2cf,it2cfp,phi,g0it,g0mat,sig,g)


      ! Get DLR coefficients of solution

      call dlr_it2cf(r,2,it2cf,it2cfp,g,gc)


      ! Get exact 3x3 Green's function on imaginary time grid by
      ! diagonalization

      allocate(gtrue(r,3,3),gctrue(r,2,2))

      call getg0_it(beta,r,3,dlrit,h,0.0d0,gtrue)


      ! Get DLR coefficients of upper left 2x2 block of exact Green's
      ! function

      call dlr_it2cf(r,2,it2cf,it2cfp,gtrue(:,1:2,1:2),gctrue)



      ! Get test points in relative format

      allocate(it_tst(ntst))

      call eqpts_rel(ntst,it_tst)


      ! Evaluate solutions on test grid

      allocate(gtst(ntst,2,2),gtruetst(ntst,2,2))

      do k=1,ntst

        call dlr_it_eval(r,2,dlrrf,gc,it_tst(k),gtst(k,:,:))
        call dlr_it_eval(r,2,dlrrf,gctrue,it_tst(k),gtruetst(k,:,:))

      enddo

      err = maxval(abs(gtst-gtruetst))

      write(6,*) 'Max L^inf error = ',err


      ! Return failed status if error is not sufficiently small

      if (err.lt.1.0d-13) then
        call exit(0)
      else
        call exit(1)
      endif

      end subroutine dyson_mat_main





      subroutine getg0_it(beta,r,n,dlrit,h,mu,g0)

      ! Obtain nxn single-particle Green's function with Hamiltonian h
      ! on DLR imaginary time grid

      implicit none
      integer r,n
      real *8 beta,dlrit(r),h(n,n),mu,g0(r,n,n)

      integer i,j,lwork,info
      real *8, allocatable :: eval(:),evec(:,:),d(:),e(:),tau(:),work(:)
      real *8, external :: kfunf_rel

      ! Get eigenvalues and eigenvectors of Hamiltonian matrix h

      lwork = n*n
      allocate(eval(n),evec(n,n),e(n-1),tau(n-1),work(lwork))

      evec = h

      call dsytrd('U',n,evec,n,eval,e,tau,work,lwork,info)
      call dorgtr('U',n,evec,n,tau,work,lwork,info)
      call dsteqr('V',n,eval,e,evec,n,work,info)

      ! Get free particle Green's function

      g0 = 0
      do i=1,r
        do j=1,n
          g0(i,j,j) = -kfunf_rel(dlrit(i),beta*(eval(j)-mu))
        enddo
          g0(i,:,:) = matmul(evec,matmul(g0(i,:,:),transpose(evec)))
      enddo

      end subroutine getg0_it
