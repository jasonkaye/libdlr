      program dlr_syk_it

      ! Solve SYK equation self-consistently by weighted fixed point
      ! iteration in imaginary time
      !
      ! The SYK equation is the Dyson equation corresponding to
      ! self-energy
      !
      ! Sigma(tau) = c^2 * G^2(tau) * G(beta-tau).
      !
      ! We solve the Dyson equation self-consistently by a weighted
      ! fixed point iteration, with weight w assigned to the new iterate
      ! and weight 1-w assigned to the previous iterate.
      !
      ! Furthermore, to solve the equation with a desired single
      ! particle energy mu, we pick a number nmu>1, and solve a sequence
      ! intermediate problems to obtain a good initial guess. First we
      ! solve the equation with single particle energy mu_0 = 0, then
      ! use this solution as an initial guess for the fixed point
      ! iteration with mu_1 = mu/nmu, then use this the solution as an
      ! initial guess for the fixed point iteration with mu = 2*mu/nmu,
      ! and so on, until we reach m_{nmu} = mu.
      !
      ! The solution will be given at nout equispaced points in
      ! imaginary time, in a file "gfun". The first column of gfun
      ! contains the imaginary time points, and the second column
      ! contains the corresponding values of the Green's function.

      
      implicit none
      integer nout,maxit,nmu
      real *8 lambda,eps,fptol,w,beta,mu,c

      ! --- Input parameters ---

      lambda = 500.0d0 ! DLR cutoff
      eps = 1.0d-14 ! Desired accuracy

      beta = 50.0d0 ! Inverse temperature
      mu = 0.1d0 ! Single particle energy
      c = 1.0d0 ! Self-energy strength

      maxit = 1000 ! Max # fixed point iterations
      fptol = 1.0d-12 ! Fixed point tolerance
      w = 0.5d0 ! Fixed point iteration weighting
      
      nmu = 1 ! # intermediate problems to solve
      
      nout = 1000 ! # points at which to output solution


      ! --- Call main test subroutine ---

      call dlr_syk_it_main(lambda,eps,nout,fptol,maxit,w,beta,mu,c,nmu)


      end program dlr_syk_it


      subroutine dlr_syk_it_main(lambda,eps,nout,fptol,maxit,w,beta,&
          mu,c,nmu)

      ! Main driver routine for imaginary time SYK solver 

      implicit none
      integer nout,maxit,nmu
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer i,j,rank,info,numit
      integer, allocatable :: it2cfpiv(:)
      real *8 one,gtest,gtest2
      real *8, allocatable :: ttst(:),it2cf(:,:),dlrit(:),dlrrf(:),g(:)
      real *8, allocatable :: cf2it(:,:),it2itr(:,:),phi(:,:),g0(:)
      real *8, external :: kfunf_rel

      one = 1.0d0

      ! Build DLR basis, grid

      rank = 500 ! Upper bound on rank

      allocate(dlrrf(rank),dlrit(rank))

      call dlr_buildit(lambda,eps,rank,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients transform matrix in LU form

      allocate(it2cf(rank,rank),it2cfpiv(rank))

      call dlr_it2cf(rank,dlrrf,dlrit,it2cf,it2cfpiv)


      ! Get DLR coefficients -> imaginary time values transform matrix,
      ! and DLR coefficients -> reflected imaginary time values
      ! transform matrix

      allocate(cf2it(rank,rank),it2itr(rank,rank))

      call dlr_cf2it(rank,dlrrf,dlrit,cf2it)

      call dlr_it2itr(rank,dlrrf,dlrit,it2cf,it2cfpiv,it2itr)


      ! --- Solve SYK equation ---

      ! Get tensor used to form matrix of convolution by a Green's
      ! function

      allocate(phi(rank*rank,rank))

      call dlr_convtens(beta,rank,dlrrf,dlrit,it2cf,it2cfpiv,phi)


      ! Get free particle Green's function

      allocate(g0(rank))

      call getg0(beta,rank,dlrit,0.0d0,g0)


      ! Solve Dyson equation by slowly increasing mu

      allocate(g(rank))

      numit = maxit

      ! Set initial guess to g0

      g = g0

      call dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,phi,&
        sigeval,w,fptol,numit,g0,g,info)

      write(6,*) 'mu = ',0.0d0

      if (info.ne.0) then
        write(6,*) 'Did not converge after ',numit,' iterations.'
      else
        write(6,*) 'Converged in ',numit,' iterations.'
      endif

      do i=1,nmu

        numit = maxit

        call getg0(beta,rank,dlrit,i*mu/nmu,g0)

        call dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,phi,&
          sigeval,w,fptol,numit,g0,g,info)

        write(6,*) 'mu = ',i*mu/nmu

        if (info.ne.0) then
          write(6,*) 'Did not converge after ',numit,' iterations.' 
        else
          write(6,*) 'Converged in ',numit,' iterations.'
        endif

      enddo


      ! Evaluate at output points and generate output file 

      allocate(ttst(nout))

      call eqpts_rel(nout,ttst)

      call dlr_expnd(rank,it2cf,it2cfpiv,g,g)

      open(1,file='gfun')

      do i=1,nout

        call dlr_eval(rank,dlrrf,g,ttst(i),gtest)

        call rel2abs(1,ttst(i),ttst(i))

        write(1,*) ttst(i),gtest

      enddo

      close(1)


      contains

        subroutine sigeval(rank,g,sig)

        ! Evaluator for SYK self-energy

        implicit none
        integer rank
        real *8 g(rank),sig(rank)

        sig = c**2*g**2*matmul(it2itr,g)

        end subroutine sigeval
      
      end subroutine dlr_syk_it_main


      subroutine getg0(beta,rank,dlrit,mu,g0)

      implicit none
      integer rank
      real *8 beta,dlrit(rank),mu,g0(rank)

      integer i
      real *8, external :: kfunf_rel

      do i=1,rank
        g0(i) = -kfunf_rel(dlrit(i),beta*mu)
      enddo

      end subroutine getg0
