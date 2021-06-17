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
      
      nout = 10000 ! # points at which to output solution


      ! --- Call main test subroutine ---

      call dlr_syk_it_main(lambda,eps,nout,fptol,maxit,w,beta,mu,c,nmu)


      end program dlr_syk_it


      subroutine dlr_syk_it_main(lambda,eps,nout,fptol,maxit,w,beta,&
          mu,c,nmu)

      ! Main driver routine for imaginary time SYK solver 

      implicit none
      integer nout,maxit,nmu
      real *8 lambda,eps,fptol,w,beta,mu,c

      integer npt,npo,p,nt,no,i,j,rank,info,pg,npg,numit
      integer, allocatable :: it2cfpiv(:),tidx(:),oidx(:)
      real *8 one,gtest,kerr(2),gtest2
      real *8, allocatable :: kmat(:,:),t(:),om(:),ttst(:)
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:),g(:)
      real *8, allocatable :: cf2it(:,:),cf2itr(:,:),phi(:,:)

      one = 1.0d0

      ! --- Build DLR basis, grid, transform matrix ---

      ! Set parameters for the fine grid based on lambda

      call gridparams(lambda,p,npt,npo,nt,no)


      ! Get fine composite Chebyshev discretization of K(tau,omega)

      allocate(kmat(nt,no),t(nt),om(no))

      call kfine_cc(lambda,p,npt,npo,t,om,kmat,kerr)


      ! Select real frequency points for DLR basis

      rank = 500 ! Upper bound on possible rank

      allocate(dlrrf(rank),oidx(rank))

      call dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx)


      ! Get DLR imaginary time grid

      allocate(dlrit(rank),tidx(rank))

      call dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit,tidx)


      ! Get imaginary time values -> DLR coefficients transform matrix in LU form

      allocate(it2cf(rank,rank),it2cfpiv(rank))

      call dlr_it2cf(nt,no,kmat,rank,oidx,tidx,it2cf,it2cfpiv)


      ! Get DLR coefficients -> imaginary time values transform matrix,
      ! and DLR coefficients -> reflected imaginary time values
      ! transform matrix

      allocate(cf2it(rank,rank),cf2itr(rank,rank))

      call dlr_cf2it(nt,no,kmat,rank,oidx,tidx,cf2it)

      call dlr_cf2itr(rank,dlrrf,dlrit,cf2itr)


      ! --- Solve SYK equation ---

      ! Get tensor used to form matrix of convolution by a Green's
      ! function

      allocate(phi(rank*rank,rank))

      call dlr_convtens(rank,dlrrf,dlrit,phi)
      phi = phi*beta


      ! Solve Dyson equation by marching in mu

      allocate(g(rank))

      numit = maxit

      call dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,phi,0*one,&
        sigeval,w,fptol,numit,1,g,info)

      write(6,*) 'mu = ',0.0d0

      if (info.ne.0) then
        write(6,*) 'Did not converge after ',numit,' iterations.'
      else
        write(6,*) 'Converged in ',numit,' iterations.'
      endif

      do i=1,nmu

        numit = maxit

        call dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,phi,&
          i*mu/nmu,sigeval,w,fptol,numit,0,g,info)

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

        sig = c**2*matmul(cf2it,g)**2*matmul(cf2itr,g)

        end subroutine sigeval
      
      end subroutine dlr_syk_it_main


      subroutine dlr_cf2itr(rank,dlrrf,dlrit,cf2itr)

      implicit none
      integer rank
      real *8 dlrrf(rank),dlrit(rank),cf2itr(rank,rank)

      integer i,j
      real *8, external :: kfunf_rel

      ! Get matrix taking DLR coefficients to values of DLR expansion at
      ! imaginary time points reflected about tau = 1/2.

      do j=1,rank
        do i=1,rank
          cf2itr(i,j) = kfunf_rel(-dlrit(i),dlrrf(j))
        enddo
      enddo

      end subroutine dlr_cf2itr

