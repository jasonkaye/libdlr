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
      
           
      program syk_it

      ! Demonstration of solution of SYK equation in imaginary time
      ! using the discrete Lehmann representation
      !
      ! The SYK equation is the nonlinear Dyson equation corresponding
      ! to self-energy
      !
      ! Sigma(tau) = c^2 * G^2(tau) * G(beta-tau).
      !
      ! We solve the Dyson equation self-consistently by a weighted
      ! fixed point iteration, with weight w assigned to the new iterate
      ! and weight 1-w assigned to the previous iterate.
      !
      ! To solve the equation with a desired single particle energy mu,
      ! we pick a number nmu>1, and solve a sequence intermediate
      ! problems to obtain a good initial guess. First we solve the
      ! equation with single particle energy mu_0 = 0, then use this
      ! solution as an initial guess for the fixed point iteration with
      ! mu_1 = mu/nmu, then use this the solution as an initial guess
      ! for the fixed point iteration with mu = 2*mu/nmu, and so on,
      ! until we reach m_{nmu} = mu.
      !
      ! The solution is output at nout equispaced points in imaginary
      ! time, in a file "gfun". The first column of gfun contains these
      ! imaginary time points, and the second column contains the
      ! corresponding values of the Green's function.

      
      implicit none
      integer nout,maxit,nmu
      real *8 lambda,eps,fptol,w,beta,mu,c

      ! Input parameters

      lambda = 500.0d0  ! DLR high energy cutoff
      eps = 1.0d-14     ! DLR error tolerance

      beta = 50.0d0     ! Inverse temperature
      mu = 0.1d0        ! Chemical potential
      c = 1.0d0         ! Self-energy strength

      maxit = 1000      ! Max # fixed point iterations
      fptol = 1.0d-12   ! Fixed point tolerance
      w = 0.5d0         ! Fixed point iteration weighting
      
      nmu = 1           ! # intermediate problems to solve
      
      nout = 1000       ! # output points


      ! Main subroutine

      call syk_it_main(lambda,eps,nout,fptol,maxit,w,beta,mu,c,nmu)


      end program syk_it


      subroutine syk_it_main(lambda,eps,nout,fptol,maxit,w,beta,&
          mu,c,nmu)

      implicit none
      integer nout,maxit,nmu
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer i,j,r,info,numit
      integer, allocatable :: it2cfp(:)
      real *8 gtest,gtest2
      real *8, allocatable :: ttst(:),it2cf(:,:),dlrit(:),dlrrf(:),g(:)
      real *8, allocatable :: it2itr(:,:),phi(:,:),g0(:)
      real *8, external :: kfunf_rel


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Get DLR coefficients -> reflected imaginary time values matrix
      ! (need this because self-energy involves a factor of the form
      ! G(beta-tau)

      allocate(it2itr(r,r))

      call dlr_it2itr_init(r,dlrrf,dlrit,it2cf,it2cfp,it2itr)


      ! Get convolution tensor

      allocate(phi(r*r,r))

      call dlr_convtens(beta,-1,r,dlrrf,dlrit,it2cf,it2cfp,phi)


      ! Get free particle Green's function on imaginary time grid

      allocate(g0(r))

      call getg0_it(beta,r,dlrit,0.0d0,g0)



      ! Solve nonlinear Dyson equation for mu = 0

      allocate(g(r))

      numit = maxit

      g = g0 ! Initial guess

      call dlr_dyson_it(beta,r,dlrit,it2cf,it2cfp,phi,&
        sigfun,w,fptol,numit,g0,g,info)

      write(6,*) 'mu = ',0.0d0

      if (info.ne.0) then
        write(6,*) 'Did not converge after ',numit,' iterations.'
      else
        write(6,*) 'Converged in ',numit,' iterations.'
      endif


      ! Solve nonlinear Dyson equation for mu = i*mu/nmu, i=1,...,nmu,
      ! using solution from previous mu as initial guess for next mu

      do i=1,nmu

        numit = maxit

        call getg0_it(beta,r,dlrit,i*mu/nmu,g0)

        call dlr_dyson_it(beta,r,dlrit,it2cf,it2cfp,phi,&
          sigfun,w,fptol,numit,g0,g,info)

        write(6,*) 'mu = ',i*mu/nmu

        if (info.ne.0) then
          write(6,*) 'Did not converge after ',numit,' iterations.' 
        else
          write(6,*) 'Converged in ',numit,' iterations.'
        endif

      enddo


      ! Get DLR coefficients of solution

      call dlr_it2cf(r,it2cf,it2cfp,g,g)


      ! Get output points in relative format

      allocate(ttst(nout))

      call eqpts_rel(nout,ttst)

      
      ! Generate output file

      open(1,file='gfun')

      do i=1,nout

        call dlr_it_eval(r,dlrrf,g,ttst(i),gtest)

        call rel2abs(1,ttst(i),ttst(i))

        write(1,*) ttst(i),gtest

      enddo

      close(1)


      contains

        subroutine sigfun(r,g,sig)

        ! Evaluate SYK self-energy Sigma = G(tau)*G(tau)*G(beta-tau)

        implicit none
        integer r
        real *8 g(r),sig(r)

        call dlr_it2itr(r,it2itr,g,sig)

        sig = c*c*g*g*sig

        end subroutine sigfun
      
      end subroutine syk_it_main


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

