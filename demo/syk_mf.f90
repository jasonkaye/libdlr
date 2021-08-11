      program dlr_syk_mf

      ! Solve SYK equation self-consistently by weighted fixed point
      ! iteration, transforming back and forth between imaginary time
      ! and Matsubara frequency
      !
      ! The SYK equation is the Dyson equation corresponding to
      ! self-energy
      !
      ! Sigma(tau) = c^2 * G^2(tau) * G(beta-tau).
      !
      ! We solve the Dyson equation self-consistently by a weighted
      ! fixed point iteration, with weight w assigned to the new iterate
      ! and weight 1-w assigned to the previous iterate. The Dyson
      ! equation is inverted in the Matsubara frequency domain, where it
      ! is diagonal, and the self-energy is computed in the imaginary
      ! time domain.
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
      integer nout,maxit,nmu,nmax
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer i

      ! --- Input parameters ---

      lambda = 500.0d0 ! DLR cutoff
      eps = 1.0d-14 ! Desired accuracy
      nmax = ceiling(lambda) ! Matsubara frequency cutoff

      beta = 50.0d0 ! Inverse temperature
      mu = 0.1d0 ! Single particle energy
      c = 1.0d0 ! Self-energy strength

      maxit = 1000 ! Max # fixed point iterations
      fptol = 1.0d-12 ! Fixed point tolerance
      w = 0.5d0 ! Fixed point iteration weighting
      
      nmu = 1 ! # intermediate problems to solve
      
      nout = 1000 ! # points at which to output solution


      ! --- Call main test subroutine ---

      call dlr_syk_mf_main(lambda,eps,nmax,nout,fptol,maxit,w,beta,mu,&
        c,nmu)


      end program dlr_syk_mf


      subroutine dlr_syk_mf_main(lambda,eps,nmax,nout,fptol,maxit,w,&
          beta,mu,c,nmu)

      ! Main driver routine for imaginary time SYK solver 

      implicit none
      integer nout,maxit,nmu,nmax
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer i,j,r,info,numit
      integer, allocatable :: dlrmf(:),it2cfp(:),mf2cfp(:)
      real *8 one,gtest,gtest2
      real *8, allocatable :: ttst(:),it2cf(:,:),dlrit(:),dlrrf(:),g(:)
      real *8, allocatable :: cf2it(:,:),it2itr(:,:)
      complex *16, allocatable :: mf2cf(:,:),cf2mf(:,:),g0(:)

      one = 1.0d0


      ! Build DLR basis, grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_buildit(lambda,eps,r,dlrrf,dlrit)


      ! Get imaginary time values -> DLR coefficients transform matrix in LU form

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Get DLR Matsubara frequency grid

      allocate(dlrmf(r))

      call dlr_mf(nmax,r,dlrrf,dlrmf)


      ! Get Matsubara frequency values -> DLR coefficients transform matrix in LU form

      allocate(mf2cf(r,r),mf2cfp(r))

      call dlr_mf2cf(nmax,r,dlrrf,dlrmf,mf2cf,mf2cfp)



      ! Get DLR coefficients -> imaginary time values transform matrix,
      ! DLR coefficients -> reflected imaginary time values
      ! transform matrix, and DLR coefficients -> Matsubara frequency
      ! values transform matrix

      allocate(cf2it(r,r),it2itr(r,r),cf2mf(r,r))

      call dlr_cf2it(r,dlrrf,dlrit,cf2it)

      call dlr_it2itr(r,dlrrf,dlrit,it2cf,it2cfp,it2itr)

      call dlr_cf2mf(r,dlrrf,dlrmf,cf2mf)
      

      ! --- Solve SYK equation ---


      ! Solve Dyson equation by marching in mu

      allocate(g(r))

      numit = maxit

      ! Get free particle Green's function

      allocate(g0(r))

      call getg0(beta,r,dlrmf,0.0d0,g0)

      ! Set initial guess to g0

      call dlr_mfexpnd(r,mf2cf,mf2cfp,g0,g)

      g = matmul(cf2it,g)

      call dlr_dyson_mf(beta,r,dlrit,it2cf,it2cfp,cf2it,&
        dlrmf,mf2cf,mf2cfp,cf2mf,sigfun,w,fptol,numit,g0,g,info)

      write(6,*) 'mu = ',0.0d0

      if (info.ne.0) then
        write(6,*) 'Did not converge after ',numit,' iterations.'
      else
        write(6,*) 'Converged in ',numit,' iterations.'
      endif

      do i=1,nmu

        numit = maxit

        call getg0(beta,r,dlrmf,i*mu/nmu,g0)

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


      ! Evaluate at output points and generate output file 

      allocate(ttst(nout))

      call eqpts_rel(nout,ttst)

      call dlr_expnd(r,it2cf,it2cfp,g,g)

      open(1,file='gfun')

      do i=1,nout

        call dlr_eval(r,dlrrf,g,ttst(i),gtest)

        call rel2abs(1,ttst(i),ttst(i))

        write(1,*) ttst(i),gtest

      enddo

      close(1)


      contains

        subroutine sigfun(r,g,sig)

        ! Evaluator for SYK self-energy

        implicit none
        integer r
        real *8 g(r),sig(r)

        sig = c**2*g**2*matmul(it2itr,g)

        end subroutine sigfun
      
      end subroutine dlr_syk_mf_main


      subroutine getg0(beta,r,dlrmf,mu,g0)

      implicit none
      integer r,dlrmf(r)
      real *8 beta,mu
      complex *16 g0(r)

      integer i
      complex *16, external :: kfunf_mf

      do i=1,r
        g0(i) = -kfunf_mf(dlrmf(i),beta*mu)
      enddo

      end subroutine getg0
