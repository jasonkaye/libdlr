      program syk_mf

      ! Demonstration of solution of SYK equation by Matsubara
      ! frequency/imaginary time alternation method using the discrete
      ! Lehmann representation
      !
      ! The SYK equation is the nonlinear Dyson equation corresponding
      ! to self-energy
      !
      ! Sigma(tau) = c^2 * G^2(tau) * G(beta-tau).
      !
      ! We solve the Dyson equation self-consistently by a weighted
      ! fixed point iteration, with weight w assigned to the new iterate
      ! and weight 1-w assigned to the previous iterate. The self-energy
      ! is evaluated in the imaginary time domain, and each linear Dyson
      ! equation, corresponding to fixed self-energy, is solved in the
      ! Matsubara frequency domain, where it is diagonal.
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
      integer nout,maxit,nmu,nmax
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer i

      ! --- Input parameters ---

      lambda = 500.0d0  ! DLR high energy cutoff
      eps = 1.0d-14     ! DLR error tolerance
      nmax = 500        ! DLR Matsubara frequency cutoff

      beta = 50.0d0     ! Inverse temperature
      mu = 0.1d0        ! Chemical potential
      c = 1.0d0         ! Self-energy strength

      maxit = 1000      ! Max # fixed point iterations
      fptol = 1.0d-12   ! Fixed point tolerance
      w = 0.5d0         ! Fixed point iteration weighting
      
      nmu = 1           ! # intermediate problems to solve
      
      nout = 1000       ! # output points


      ! Main subroutine

      call syk_mf_main(lambda,eps,nmax,nout,fptol,maxit,w,beta,mu,&
        c,nmu)


      end program syk_mf


      subroutine syk_mf_main(lambda,eps,nmax,nout,fptol,maxit,w,&
          beta,mu,c,nmu)

      implicit none
      integer nout,maxit,nmu,nmax
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer i,j,r,info,numit
      integer, allocatable :: dlrmf(:),it2cfp(:),mf2cfp(:)
      real *8 gtest,gtest2
      real *8, allocatable :: ttst(:),it2cf(:,:),dlrit(:),dlrrf(:),g(:)
      real *8, allocatable :: cf2it(:,:),it2itr(:,:)
      complex *16, allocatable :: mf2cf(:,:),cf2mf(:,:),g0(:)


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_buildit(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Get DLR Matsubara frequency grid

      allocate(dlrmf(r))

      call dlr_mf(nmax,r,dlrrf,-1,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients matrix (LU
      ! form)

      allocate(mf2cf(r,r),mf2cfp(r))

      call dlr_mf2cf(nmax,r,dlrrf,dlrmf,-1,mf2cf,mf2cfp)


      ! Get DLR coefficients -> imaginary time values matrix

      allocate(cf2it(r,r))

      call dlr_cf2it(r,dlrrf,dlrit,cf2it)


      ! Get DLR coefficients -> reflected imaginary time values matrix
      ! (need this because self-energy involves a factor of the form
      ! G(beta-tau)

      allocate(it2itr(r,r))

      call dlr_it2itr(r,dlrrf,dlrit,it2cf,it2cfp,it2itr)


      ! Get DLR coefficients -> Matsubara frequency values matrix

      allocate(cf2mf(r,r))

      call dlr_cf2mf(r,dlrrf,dlrmf,-1,cf2mf)


      ! Get free particle Green's function on Matsubara frequency grid

      allocate(g0(r))

      call getg0_mf(beta,r,dlrmf,0.0d0,g0)



      ! Solve nonlinear Dyson equation for mu = 0

      allocate(g(r))

      numit = maxit

      ! Set initial guess to g0

      call getg0_it(beta,r,dlrit,0.0d0,g)

      call dlr_dyson_mf(beta,r,dlrit,it2cf,it2cfp,cf2it,&
        dlrmf,mf2cf,mf2cfp,cf2mf,sigfun,w,fptol,numit,g0,g,info)

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

        call getg0_mf(beta,r,dlrmf,i*mu/nmu,g0)

        call dlr_dyson_mf(beta,r,dlrit,it2cf,it2cfp,cf2it,&
          dlrmf,mf2cf,mf2cfp,cf2mf,sigfun,w,fptol,numit,g0,&
          g,info)

        write(6,*) 'mu = ',i*mu/nmu

        if (info.ne.0) then
          write(6,*) 'Did not converge after ',numit,' iterations.' 
        else
          write(6,*) 'Converged in ',numit,' iterations.'
        endif

      enddo


      ! Get DLR coefficients of solution

      call dlr_expnd(r,it2cf,it2cfp,g,g)


      ! Get output points in relative format

      allocate(ttst(nout))

      call eqpts_rel(nout,ttst)

      
      ! Generate output file

      open(1,file='gfun')

      do i=1,nout

        call dlr_eval(r,dlrrf,g,ttst(i),gtest)

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

        sig = c**2*g**2*matmul(it2itr,g)

        end subroutine sigfun

      end subroutine syk_mf_main


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

