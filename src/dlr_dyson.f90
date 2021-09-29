      module dlr_dyson_mod
        use dlr_it_mod
        use dlr_mf_mod
        use dlr_conv_mod
        implicit none
      contains

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
      !! @param[in]     beta    inverse temperature
      !! @param[in]     r       number of DLR basis functions
      !! @param[in]     dlrit   DLR imaginary time nodes
      !! @param[in]     it2cf   imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]     it2cfp  imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[in]     cf2it   DLR coefficients -> imaginary time grid
      !!                          values transform matrix
      !! @param[in]     phi     tensor taking DLR coefficients of g to
      !!                          matrix of convolution by g.
      !! @param[in]     sigfun  subroutine with calling sequence
      !!                          sigfun(r,g,sig), which takes in
      !!                          value of an imaginary time Green's
      !!                          function at the imaginary time grid
      !!                          points, and returns values of the
      !!                          self-energy Sigma at those grid points
      !! @param[in]     w       weighting parameter for fixed point
      !!                          iteration
      !! @param[in]     fptol   fixed point iteration tolerance
      !! @param[in,out] numit   on input: max number of fixed point
      !!                          iterations; on output: number of
      !!                          fixed point iterations taken
      !! @param[in]     g0      right hand side of Dyson equation, on
      !!                          imaginary time grid
      !! @param[in,out] g       on input, initial guess for fixed
      !!                          point iteration; on output, solution
      !!                          of the Dyson equation on imaginary
      !!                          time grid
      !! @param[out]    info    =0 if iteration converged to tolerance
      !!                          fptol; =-1 if iteration did not
      !!                          converge

      subroutine dlr_dyson_it(beta,r,dlrit,it2cf,it2cfp,cf2it,&
          phi,sigfun,w,fptol,numit,g0,g,info)

      implicit none
      integer r,numit,it2cfp(r),info
      real *8 beta,dlrit(r),it2cf(r,r)
      real *8 cf2it(r,r),g0(r),g(r),w,fptol
      real *8 phi(r*r,r)

      integer i,info1
      real *8 one
      real *8, allocatable :: g0mat(:,:),sig(:),gnew(:)

      one = 1.0d0

      ! Get matrix of convolution by G0

      allocate(g0mat(r,r))

      call dlr_convmat(r,it2cf,it2cfp,phi,g0,g0mat)


      ! Weighted fixed point iteration

      allocate(sig(r),gnew(r))

      do i=1,numit

        ! Evaluate self-energy

        call sigfun(r,g,sig)

        ! Solve linear Dyson equation

        call dyson_it_lin(r,it2cf,it2cfp,phi,g0,g0mat,sig,gnew)

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

      subroutine dyson_it_lin(r,it2cf,it2cfp,phi,g0,g0mat,sig,g)

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
      !! @param[in]     beta    inverse temperature
      !! @param[in]     r       number of DLR basis functions
      !! @param[in]     dlrit   DLR imaginary time nodes
      !! @param[in]     it2cf   imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]     it2cfp  imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[in]     cf2it   DLR coefficients -> imaginary time grid
      !!                          values transform matrix
      !! @param[in]     dlrmf   DLR Matsubara frequency nodes
      !! @param[in]     mf2cf   Matsubara frequency grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]     mf2cfp  Matsubra frequency grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[in]     cf2mf   DLR coeffs -> Matsubara freq grid
      !!                          values transform matrix
      !! @param[in]     sigfun  subroutine with calling sequence
      !!                          sigfun(r,g,sig), which takes in
      !!                          value of an imaginary time Green's
      !!                          function at the imaginary time grid
      !!                          points, and returns values of the
      !!                          self-energy Sigma at those grid
      !!                          points
      !! @param[in]     w       weighting parameter for fixed point
      !!                          iteration
      !! @param[in]     fptol   fixed point iteration tolerance
      !! @param[in,out] numit   on input: max number of fixed point
      !!                          iterations; on output: number of
      !!                          fixed point iterations taken
      !! @param[in]     g0      right hand side of Dyson equation, on
      !!                          Matsubara frequency grid
      !! @param[in,out] g       on input, initial guess for fixed
      !!                          point iteration; on output, solution
      !!                          of the Dyson equation on imaginary
      !!                          time grid
      !! @param[out]    info    =0 if iteration converged to tolerance
      !!                          fptol; =-1 if iteration did not
      !!                          converge

      subroutine dlr_dyson_mf(beta,r,dlrit,it2cf,it2cfp,cf2it,&
          dlrmf,mf2cf,mf2cfp,cf2mf,sigfun,w,fptol,numit,g0,g,&
          info)

      implicit none
      integer r,numit,it2cfp(r),mf2cfp(r),dlrmf(r)
      integer info
      real *8 beta,dlrit(r),it2cf(r,r)
      real *8 cf2it(r,r),g(r),w,fptol
      complex *16 mf2cf(r,r),cf2mf(r,r),g0(r)

      integer i
      real *8 one
      real *8, allocatable :: sig(:),gnew(:)

      one = 1.0d0

      ! Weighted fixed point iteration

      allocate(sig(r),gnew(r))

      do i=1,numit

        ! Evaluate self-energy

        call sigfun(r,g,sig)

        ! Solve linear Dyson equation

        call dyson_mf_lin(beta,r,it2cf,it2cfp,cf2it,&
          mf2cf,mf2cfp,cf2mf,g0,sig,gnew)

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

      end subroutine dlr_dyson_mf





      !> Solve the linear Dyson equation in Matsubara frequency, with a
      !! fixed self-energy.
      !!
      !! This solver takes in the self-energy and returns the imaginary
      !! time Green's function on the imaginary time grid, but performs
      !! the solve in the Matsubara frequency domain by diagonal
      !! inversion.
      !!
      !! @param[in]   beta    inverse temperature
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   it2cf   imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU factors
      !! @param[in]   it2cfp  imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU pivots
      !! @param[in]   cf2it   DLR coefficients -> imaginary time grid
      !!                        values transform matrix
      !! @param[in]   mf2cf   Matsubara frequency grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU factors
      !! @param[in]   mf2cfp  Matsubra frequency grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU pivots
      !! @param[in]   cf2mf   DLR coeffs -> Matsubara freq grid
      !!                        values transform matrix
      !! @param[in]   g0      values of the right hand side G0 on
      !!                        the Matsubara frequency grid
      !! @param[in]   sig     values of the self-energy on the
      !!                          imaginary time grid
      !! @param[out]  g       solution of the linear Dyson equation
      !!                          on imaginary time grid

      subroutine dyson_mf_lin(beta,r,it2cf,it2cfp,cf2it,&
          mf2cf,mf2cfp,cf2mf,g0,sig,g)

      implicit none
      integer r,it2cfp(r),mf2cfp(r)
      real *8 beta,it2cf(r,r),cf2it(r,r)
      real *8 sig(r),g(r)
      complex *16 mf2cf(r,r),cf2mf(r,r),g0(r)

      real *8 one
      real *8, allocatable :: sigc(:),tmp2(:)
      complex *16, allocatable :: sigmf(:),gmf(:),tmp(:)

      one = 1.0d0

      allocate(sigc(r),gmf(r),sigmf(r),tmp(r),tmp2(r))

      ! Get DLR coefficients of self-energy

      call dlr_it_expnd(r,it2cf,it2cfp,sig,sigc)

      ! Get self-energy on Matsubara frequency grid

      tmp = sigc
      call zgemv('N',r,r,one,cf2mf,r,tmp,1,0*one,sigmf,1)

      ! Solve Dyson equation by diagonal inversion

      gmf = g0/(one-beta**2*g0*sigmf)

      ! Get DLR coefficients of solution

      call dlr_mf_expnd(r,mf2cf,mf2cfp,gmf,g)

      ! Evaluate solution on imaginary time grid

      tmp2 = g
      call dgemv('N',r,r,one,cf2it,r,tmp2,1,0*one,g,1)

      end subroutine dyson_mf_lin

      end module dlr_dyson_mod
      
