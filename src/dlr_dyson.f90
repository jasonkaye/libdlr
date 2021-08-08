      !
      !
      ! This file contains subroutines for solving the Dyson equation
      ! using the discrete Lehmann representation 
      !
      !




      !> Solve the nonlinear Dyson equation in imaginary time,
      !! with a given expression for the self-energy in terms of the
      !! Green's function.
      !!
      !! This solver uses weighted fixed point iteration to achieve
      !! self-consistency; this is an iteration
      !! G_{n+1} = w*G_n + (1-w)*G_{n-1}, where w is a given weight. The
      !! iteration terminates when the difference between G_{n+1} and
      !! G_n at all of the imaginary time grid points is less than a
      !! given tolerance in absolute value, or when a maximum number of
      !! iterations is reached.
      !!
      !! @param[in]     beta      inverse temperature
      !! @param[in]     rank      number of DLR basis functions
      !! @param[in]     dlrit     DLR imaginary time nodes
      !! @param[in]     it2cf     imaginary time grid values ->
      !!                            DLR coefficients transform matrix,
      !!                            stored in LAPACK LU factored format;
      !!                            LU factors
      !! @param[in]     it2cfpiv  imaginary time grid values ->
      !!                            DLR coefficients transform matrix,
      !!                            stored in LAPACK LU factored format;
      !!                            LU pivots
      !! @param[in]     cf2it     DLR coefficients -> imaginary time grid
      !!                            values transform matrix
      !! @param[in]     phi       tensor taking DLR coefficients of g to
      !!                            matrix of convolution by g.
      !! @param[in]     sigeval   subroutine with calling sequence
      !!                            sigeval(rank,g,sig), which takes in
      !!                            value of an imaginary time Green's
      !!                            function at the imaginary time grid
      !!                            points, and returns values of the
      !!                            self-energy Sigma at those grid points
      !! @param[in]     w         weighting parameter for fixed point
      !!                            iteration
      !! @param[in]     fptol     fixed point iteration tolerance
      !! @param[in,out] numit     on input: max number of fixed point
      !!                            iterations; on output: number of
      !!                            fixed point iterations taken
      !! @param[in]     g0        right hand side of Dyson equation, on
      !!                            imaginary time grid
      !! @param[in,out] g         on input, initial guess for fixed
      !!                            point iteration; on output, solution
      !!                            of the Dyson equation on imaginary
      !!                            time grid
      !! @param[out]    info      =0 if iteration converged to tolerance
      !!                            fptol; =-1 if iteration did not
      !!                            converge

      subroutine dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,&
          phi,sigeval,w,fptol,numit,g0,g,info)

      implicit none
      integer rank,numit,it2cfpiv(rank),info
      real *8 beta,dlrit(rank),it2cf(rank,rank)
      real *8 cf2it(rank,rank),g0(rank),g(rank),w,fptol
      real *8 phi(rank*rank,rank)

      integer i,info1
      real *8 one
      real *8, allocatable :: g0mat(:,:),sig(:),gnew(:)

      one = 1.0d0

      ! --- Get matrix of convolution by G0 ---

      allocate(g0mat(rank,rank))

      call dlr_convmat(rank,it2cf,it2cfpiv,phi,g0,g0mat)


      ! --- Fixed point iteration for G ---

      allocate(sig(rank),gnew(rank))

      do i=1,numit

        ! Get Sigma in imaginary time grid representation

        call sigeval(rank,g,sig)

        ! Solve linear Dyson equation

        call dyson_it_lin(rank,it2cf,it2cfpiv,phi,g0,g0mat,sig,gnew)

        ! Check self-consistency

        if (maxval(abs(gnew-g))<fptol) then

          g = gnew
          numit = i
          info = 0

          return

        else

          ! Next G is weighted linear combination of previous and
          ! current iterates

          g = w*gnew + (one-w)*g

        endif
      enddo

      info = -1

      end subroutine dlr_dyson_it





      !> Solve the linear Dyson equation in imaginary time, with a fixed
      !! self-energy.
      !!
      !! This solver forms the Dyson equation in imaginary time using
      !! the DLR basis and solves it using Gaussian elimination.
      !!
      !! @param[in]   rank      number of DLR basis functions
      !! @param[in]   it2cf     imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]   it2cfpiv  imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[in]   phi       tensor taking DLR coefficients of g to
      !!                          matrix of convolution by g.
      !! @param[in]   g0        values of the right hand side G0 on
      !!                          the imaginary time grid
      !! @param[in]   g0mat     matrix of convolution by G0
      !! @param[in]   sig       values of the self-energy on the
      !!                            imaginary time grid
      !! @param[out]  g         solution of the linear Dyson equation
      !!                            on imaginary time grid

      subroutine dyson_it_lin(rank,it2cf,it2cfpiv,phi,g0,g0mat,sig,g)

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 it2cf(rank,rank)
      real *8 g0(rank),g0mat(rank,rank),sig(rank),g(rank)
      real *8 phi(rank*rank,rank)

      integer j,info1
      integer, allocatable :: ipiv(:)
      real *8 one
      real *8, allocatable :: sigc(:),sysmat(:,:)

      one = 1.0d0

      allocate(sigc(rank),sysmat(rank,rank),ipiv(rank))

      ! Matrix of convolution by Sigma

      call dlr_convmat(rank,it2cf,it2cfpiv,phi,sig,sysmat)

      ! System matrix for linear Dyson equation

      sysmat = -matmul(g0mat,sysmat)

      do j=1,rank
        sysmat(j,j) = one + sysmat(j,j)
      enddo

      ! Solve linear equation by LU factorization + backsolve

      call dgetrf(rank,rank,sysmat,rank,ipiv,info1)

      g = g0

      call dgetrs('N',rank,1,sysmat,rank,ipiv,g,rank,info1)


      end subroutine dyson_it_lin





      !> Solve the nonlinear Dyson equation in Matsubara frequency,
      !! with a given expression for the self-energy in terms of the
      !! Green's function, evaluated in imaginary time.
      !!
      !! This solver uses weighted fixed point iteration to achieve
      !! self-consistency; this is an iteration
      !! G_{n+1} = w*G_n + (1-w)*G_{n-1}, where w is a given weight. The
      !! iteration terminates when the difference between G_{n+1} and
      !! G_n at all of the imaginary time grid points is less than a
      !! given tolerance in absolute value, or when a maximum number of
      !! iterations is reached.
      !!
      !! The solver transforms back and forth between the imaginary time
      !! and Matsubara frequency domains, evaluating the self-energy in
      !! imaginary time, and solving the Dyson equation in Matsubara
      !! frequency.
      !!
      !! @param[in]     beta      inverse temperature
      !! @param[in]     rank      number of DLR basis functions
      !! @param[in]     dlrit     DLR imaginary time nodes
      !! @param[in]     it2cf     imaginary time grid values ->
      !!                            DLR coefficients transform matrix,
      !!                            stored in LAPACK LU factored format;
      !!                            LU factors
      !! @param[in]     it2cfpiv  imaginary time grid values ->
      !!                            DLR coefficients transform matrix,
      !!                            stored in LAPACK LU factored format;
      !!                            LU pivots
      !! @param[in]     cf2it     DLR coefficients -> imaginary time grid
      !!                            values transform matrix
      !! @param[in]     dlrmf     DLR Matsubara frequency nodes
      !! @param[in]     mf2cf     Matsubara frequency grid values ->
      !!                            DLR coefficients transform matrix,
      !!                            stored in LAPACK LU factored format;
      !!                            LU factors
      !! @param[in]     mf2cfpiv  Matsubra frequency grid values ->
      !!                            DLR coefficients transform matrix,
      !!                            stored in LAPACK LU factored format;
      !!                            LU pivots
      !! @param[in]     cf2mf     DLR coeffs -> Matsubara freq grid
      !!                            values transform matrix
      !! @param[in]     sigeval   subroutine with calling sequence
      !!                            sigeval(rank,g,sig), which takes in
      !!                            value of an imaginary time Green's
      !!                            function at the imaginary time grid
      !!                            points, and returns values of the
      !!                            self-energy Sigma at those grid
      !!                            points
      !! @param[in]     w         weighting parameter for fixed point
      !!                            iteration
      !! @param[in]     fptol     fixed point iteration tolerance
      !! @param[in,out] numit     on input: max number of fixed point
      !!                            iterations; on output: number of
      !!                            fixed point iterations taken
      !! @param[in]     g0        right hand side of Dyson equation, on
      !!                            Matsubara frequency grid
      !! @param[in,out] g         on input, initial guess for fixed
      !!                            point iteration; on output, solution
      !!                            of the Dyson equation on imaginary
      !!                            time grid
      !! @param[out]    info      =0 if iteration converged to tolerance
      !!                            fptol; =-1 if iteration did not
      !!                            converge

      subroutine dlr_dyson_mf(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,&
          dlrmf,mf2cf,mf2cfpiv,cf2mf,sigeval,w,fptol,numit,g0,g,&
          info)

      implicit none
      integer rank,numit,it2cfpiv(rank),mf2cfpiv(rank),dlrmf(rank)
      integer info
      real *8 beta,dlrit(rank),it2cf(rank,rank)
      real *8 cf2it(rank,rank),g(rank),w,fptol
      complex *16 mf2cf(rank,rank),cf2mf(rank,rank),g0(rank)

      integer i
      real *8 one
      real *8, allocatable :: sig(:),gnew(:)

      one = 1.0d0

      ! --- Fixed point iteration for G ---

      allocate(sig(rank),gnew(rank))

      do i=1,numit

        ! Get Sigma in imaginary time grid representation from DLR
        ! coefficient representation of previous G

        call sigeval(rank,g,sig)

        ! Solve linear Dyson equation

        call dyson_mf_lin(beta,rank,it2cf,it2cfpiv,cf2it,&
          mf2cf,mf2cfpiv,cf2mf,g0,sig,gnew)

        ! Check self-consistency

        if (maxval(abs(gnew-g))<fptol) then

          g = gnew
          numit = i
          info = 0

          return

        else

          g = w*gnew + (one-w)*g

        endif
      enddo

      info = -1

      end subroutine dlr_dyson_mf





      !> Solve the linear Dyson equation in Matsubara frequency, with a
      !! fixed self-energy.
      !!
      !! This solver takes in the self-energy and returns the imaginary
      !! time Green's function on the imaginary time grid, but performs
      !! the solve in the Matsubara frequency domain by diagonal
      !! inversion.
      !!
      !! @param[in]   beta      inverse temperature
      !! @param[in]   rank      number of DLR basis functions
      !! @param[in]   it2cf     imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]   it2cfpiv  imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[in]   cf2it     DLR coefficients -> imaginary time grid
      !!                          values transform matrix
      !! @param[in]   mf2cf     Matsubara frequency grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]   mf2cfpiv  Matsubra frequency grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[in]   cf2mf     DLR coeffs -> Matsubara freq grid
      !!                          values transform matrix
      !! @param[in]   g0        values of the right hand side G0 on
      !!                          the Matsubara frequency grid
      !! @param[in]   sig       values of the self-energy on the
      !!                            imaginary time grid
      !! @param[out]  g         solution of the linear Dyson equation
      !!                            on imaginary time grid

      subroutine dyson_mf_lin(beta,rank,it2cf,it2cfpiv,cf2it,&
          mf2cf,mf2cfpiv,cf2mf,g0,sig,g)

      implicit none
      integer rank,it2cfpiv(rank),mf2cfpiv(rank)
      real *8 beta,it2cf(rank,rank),cf2it(rank,rank)
      real *8 sig(rank),g(rank)
      complex *16 mf2cf(rank,rank),cf2mf(rank,rank),g0(rank)

      integer j,info1
      integer, allocatable :: ipiv(:)
      real *8 one
      real *8, allocatable :: sigc(:)
      complex *16, allocatable :: sigmf(:),gmf(:)

      one = 1.0d0

      allocate(sigc(rank),gmf(rank))

      ! DLR coefficient representation of Sigma

      call dlr_expnd(rank,it2cf,it2cfpiv,sig,sigc)

      ! Matsubara frequency representation of Sigma

      sigmf = matmul(cf2mf,sigc)

      ! Invert linear Dyson equation

      gmf = g0/(one-beta**2*g0*sigmf)

      ! DLR coefficient representation of solution

      call dlr_mfexpnd(rank,mf2cf,mf2cfpiv,gmf,g)

      ! Evaluate on imaginary time grid

      g = matmul(cf2it,g)

      end subroutine dyson_mf_lin
