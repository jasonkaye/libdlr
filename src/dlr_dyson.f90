      !
      !
      ! This file contains subroutines for solving the Dyson equation
      ! using the discrete Lehmann representation 
      !
      !


      subroutine dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,&
          phi,mu,sigeval,w,fptol,numit,useg0,g,info)

      ! Solve the Dyson equation by weighted fix point iteration using
      ! the DLR in imaginary time
      !
      ! Input:
      !
      ! beta      - inverse temperature
      ! rank      - rank of DLR (# basis functions)
      ! dlrit     - selected imaginary time nodes (tau points)
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! cf2it     - DLR coefficients -> imaginary time grid values
      !               transform
      ! phi       - tensor taking a set of DLR coefficients to matrix of
      !               convolution by corresponding DLR expansion
      ! mu        - single particle energy parameter
      ! sigeval   - subroutine with calling sequence
      !               sigeval(rank,g,sig), which takes in DLR
      !               coefficients of Green's function G and returns
      !               values of Sigma at DLR imaginary time grid points.
      ! w         - weighting parameter for weighted fixed point
      !               iteration
      ! fptol     - fixed point iteration tolerance
      ! numit     - max number of fixed point iterations
      ! useg0     - useg0==1, use G0(tau) as initial guess in fixed
      !               point iteration
      ! g         - if useg0/=0, DLR coefficients of initial guess for
      !               fixed point iteration
      !
      ! Output:
      !
      ! numit     - Number of fixed point iterations taken
      ! g         - DLR coefficients of self-consistent solution of
      !               Dyson equation
      ! info      - info=0: iteration converged to tolerance fptol;
      !               info=-1 iteration did not converge

      implicit none
      integer rank,numit,useg0,it2cfpiv(rank),info
      real *8 beta,dlrit(rank),it2cf(rank,rank),mu
      real *8 cf2it(rank,rank),g(rank),w,fptol,phi(rank*rank,rank)

      integer i,info1
      real *8 one
      real *8, allocatable :: g0(:),g0mat(:,:),sig(:),gnew(:)

      one = 1.0d0

      ! --- Get DLR of G0, matrix of convolution by G0 ---

      allocate(g0(rank),g0mat(rank,rank))

      call g0init_it(beta,rank,dlrit,it2cf,it2cfpiv,phi,mu,g0,g0mat)


      ! --- Fixed point iteration for G ---

      ! Either use input initial guess, or use G0 as initial guess

      if (useg0==1) then

        g = g0

      endif

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


      subroutine dyson_it_lin(rank,it2cf,it2cfpiv,phi,&
          g0,g0mat,sig,g)

      ! Solve the linear Dyson equation, with given self-energy, using
      ! DLR in imaginary time
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! phi       - tensor taking a set of DLR coefficients to matrix of
      !               convolution by corresponding DLR expansion
      ! g0        - DLR coefficients of g0, produced by subroutine
      !               g0init_it
      ! g0mat     - matrix of convolution by G0, produced by subroutine
      !               g0init_it
      ! sig       - values of self-energy Sigma at DLR imaginary time
      !               nodes.
      !
      ! Output:
      !
      ! g         - DLR coefficients of solution of linear Dyson
      !               equation 

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

      ! DLR coefficient representation of Sigma
      
      !call dlr_expnd(rank,it2cf,it2cfpiv,sig,sigc)

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


      subroutine g0init_it(beta,rank,dlrit,it2cf,it2cfpiv,phi,mu,g0,&
          g0mat)

      ! Get DLR coefficient representation of G0, and matrix of
      ! convolution by G0, for use in Dyson solver
      !
      ! Input:
      !
      ! beta      - inverse temperature
      ! rank      - rank of DLR (# basis functions)
      ! dlrit     - selected imaginary time nodes (tau points)
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! phi       - tensor taking a set of DLR coefficients to matrix of
      !               convolution by corresponding DLR expansion
      ! mu        - single particle energy parameter
      !
      ! Output:
      !
      ! g0        - values of G0 on DLR imaginary time grid
      ! g0mat     - matrix of convolution by G0

      integer rank,it2cfpiv(rank)
      real *8 beta,dlrit(rank),it2cf(rank,rank),phi(rank*rank,rank),mu
      real *8 g0(rank),g0mat(rank,rank)

      integer i
      real *8, allocatable :: g0c(:)
      real *8, external :: kfunf_rel

      ! Imaginary time grid representation of G0

      do i=1,rank

        g0(i) = -kfunf_rel(dlrit(i),beta*mu)

      enddo

      ! Matrix of convolution by G0

      call dlr_convmat(rank,it2cf,it2cfpiv,phi,g0,g0mat)

      end subroutine g0init_it




      subroutine dlr_dyson_mf(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,&
          dlrmf,mf2cf,mf2cfpiv,cf2mf,mu,sigeval,w,fptol,numit,useg0,g,&
          info)

      ! Solve the Dyson equation by weighted fix point iteration using
      ! the DLR, computing the self-energy in the imaginary time domain
      ! and inverting the Dyson equation in the Matsubara frequency
      ! domain
      !
      ! Input:
      !
      ! beta      - inverse temperature
      ! rank      - rank of DLR (# basis functions)
      ! dlrit     - selected imaginary time nodes (tau points)
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! cf2it     - DLR coefficients -> imaginary time grid values
      !               transform
      ! dlrmf     - selected Matsubara frequency nodes
      ! mf2cf     - Matsubara frequency grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! mf2cfpiv  - pivot matrix for mf2cf in lapack LU storage format
      ! cf2mf     - DLR coefficients -> Matsubara frequency grid values
      ! mu        - single particle energy parameter
      ! sigeval   - subroutine with calling sequence
      !               sigeval(rank,g,sig), which takes in DLR
      !               coefficients of Green's function G and returns
      !               values of Sigma at DLR imaginary time grid points.
      ! w         - weighting parameter for weighted fixed point
      !               iteration
      ! fptol     - fixed point iteration tolerance
      ! numit     - max number of fixed point iterations
      ! useg0     - useg0==1, use G0(tau) as initial guess in fixed
      !               point iteration
      ! g         - if useg0/=0, DLR coefficients of initial guess for
      !               fixed point iteration
      !
      ! Output:
      !
      ! numit     - Number of fixed point iterations taken
      ! g         - DLR coefficients of self-consistent solution of
      !               Dyson equation
      ! info      - info=0: iteration converged to tolerance fptol;
      !               info=-1 iteration did not converge

      implicit none
      integer rank,numit,useg0,it2cfpiv(rank),mf2cfpiv(rank),dlrmf(rank)
      integer info
      real *8 beta,dlrit(rank),mu,it2cf(rank,rank)
      real *8 cf2it(rank,rank),g(rank),w,fptol
      complex *16 mf2cf(rank,rank),cf2mf(rank,rank)

      integer i
      real *8 one
      real *8, allocatable :: sig(:),gnew(:)
      complex *16, allocatable :: g0(:)

      one = 1.0d0

      ! --- Get Matsubara frequency representation of G0 ---

      allocate(g0(rank))

      call g0init_mf(beta,rank,dlrmf,mu,g0)


      ! --- Fixed point iteration for G ---

      ! Either use input initial guess, or use G0 as initial guess

      if (useg0==1) then
        
        call dlr_mfexpnd(rank,mf2cf,mf2cfpiv,g0,g)

      endif

      allocate(sig(rank),gnew(rank))

      do i=1,numit

        ! Get Sigma in imaginary time grid representation from DLR
        ! coefficient representation of previous G

        call sigeval(rank,g,sig)

        ! Solve linear Dyson equation

        call dyson_mf_lin(beta,rank,it2cf,it2cfpiv,mf2cf,mf2cfpiv,&
          cf2mf,g0,sig,gnew)

        ! Check self-consistency

        if (maxval(abs(matmul(cf2it,gnew-g)))<fptol) then

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


      subroutine g0init_mf(beta,rank,dlrmf,mu,g0)

      ! Get G0 on Matsubara frequency DLR grid, for use in Dyson solver
      !
      ! Input:
      !
      ! beta      - inverse temperature
      ! rank      - rank of DLR (# basis functions)
      ! dlrit     - selected imaginary time nodes (tau points)
      ! dlrmf     - selected Matsubara frequency nodes
      ! mu        - single particle energy parameter
      !
      ! Output:
      !
      ! g0        - G0 on Matsubara frequency DLR grid

      integer rank,dlrmf(rank)
      real *8 beta,mu
      complex *16 g0(rank)

      integer i
      complex *16, external :: kfunf_mf

      ! Evaluate G0 on Matsubara frequency DLR grid 

      do i=1,rank

        g0(i) = -kfunf_mf(dlrmf(i),beta*mu)

      enddo

      end subroutine g0init_mf


      subroutine dyson_mf_lin(beta,rank,it2cf,it2cfpiv,mf2cf,mf2cfpiv,&
          cf2mf,g0,sig,g)

      ! Solve the linear Dyson equation, with given self-energy, using
      ! DLR in imaginary time
      !
      ! Input:
      !
      ! beta      - inverse temperature
      ! rank      - rank of DLR (# basis functions)
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! mf2cf     - Matsubara frequency grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! mf2cfpiv  - pivot matrix for mf2cf in lapack LU storage format
      ! cf2mf     - DLR coefficients -> Matsubara frequency grid values
      ! g0        - Values of G0 on Matsubara frequency DLR grid,
      !               produced by subroutine g0init_mf
      ! sig       - values of self-energy Sigma at DLR imaginary time
      !               nodes.
      !
      ! Output:
      !
      ! g         - DLR coefficients of solution of linear Dyson
      !               equation 

      implicit none
      integer rank,it2cfpiv(rank),mf2cfpiv(rank)
      real *8 beta,it2cf(rank,rank)
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


      end subroutine dyson_mf_lin





!      subroutine dlr_dyson_it(beta,rank,dlrit,it2cf,it2cfpiv,cf2it,&
!          phi,mu,sigeval,w,fptol,numit,useg0,g,info)
!
!      ! Solve the Dyson equation by weighted fix point iteration using
!      ! the DLR in imaginary time
!      !
!      ! Input:
!      !
!      ! beta      - inverse temperature
!      ! rank      - rank of DLR (# basis functions)
!      ! dlrit     - selected imaginary time nodes (tau points)
!      ! it2cf     - imaginary time grid values -> DLR coefficients
!      !               transform matrix in lapack LU storage format
!      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
!      !               format
!      ! cf2it     - DLR coefficients -> imaginary time grid values
!      !               transform
!      ! phi       - tensor taking a set of DLR coefficients to matrix of
!      !               convolution by corresponding DLR expansion
!      ! mu        - single particle energy parameter
!      ! sigeval   - subroutine with calling sequence
!      !               sigeval(rank,g,sig), which takes in DLR
!      !               coefficients of Green's function G and returns
!      !               values of Sigma at DLR imaginary time grid points.
!      ! w         - weighting parameter for weighted fixed point
!      !               iteration
!      ! fptol     - fixed point iteration tolerance
!      ! numit     - max number of fixed point iterations
!      ! useg0     - useg0==1, use G0(tau) as initial guess in fixed
!      !               point iteration
!      ! g         - if useg0/=0, DLR coefficients of initial guess for
!      !               fixed point iteration
!      !
!      ! Output:
!      !
!      ! numit     - Number of fixed point iterations taken
!      ! g         - DLR coefficients of self-consistent solution of
!      !               Dyson equation
!      ! info      - info=0: iteration converged to tolerance fptol;
!      !               info=-1 iteration did not converge
!
!      implicit none
!      integer rank,numit,useg0,it2cfpiv(rank),info
!      real *8 beta,dlrit(rank),it2cf(rank,rank),mu
!      real *8 cf2it(rank,rank),g(rank),w,fptol,phi(rank*rank,rank)
!
!      integer i,info1
!      real *8 one
!      real *8, allocatable :: g0(:),g0mat(:,:),sig(:),gnew(:)
!
!      one = 1.0d0
!
!      ! --- Get DLR of G0, matrix of convolution by G0 ---
!
!      allocate(g0(rank),g0mat(rank,rank))
!
!      call g0init_it(beta,rank,dlrit,it2cf,it2cfpiv,phi,mu,g0,g0mat)
!
!
!      ! --- Fixed point iteration for G ---
!
!      ! Either use input initial guess, or use G0 as initial guess
!
!      if (useg0==1) then
!
!        g = g0
!
!      endif
!
!      allocate(sig(rank),gnew(rank))
!
!      do i=1,numit
!
!        ! Get Sigma in imaginary time grid representation from DLR
!        ! coefficient representation of previous G
!
!        call sigeval(rank,g,sig)
!
!        ! Solve linear Dyson equation
!
!        call dyson_it_lin(rank,it2cf,it2cfpiv,phi,g0,g0mat,sig,gnew)
!
!        ! Check self-consistency
!
!        if (maxval(abs(matmul(cf2it,gnew-g)))<fptol) then
!
!          g = gnew
!          numit = i
!          info = 0
!
!          return
!
!        else
!
!          ! Next G is weighted linear combination of previous and
!          ! current iterates
!
!          g = w*gnew + (one-w)*g
!
!        endif
!      enddo
!
!      info = -1
!
!      end subroutine dlr_dyson_it
!
!
!
!
!      subroutine dyson_it_lin(rank,it2cf,it2cfpiv,phi,&
!          g0,g0mat,sig,g)
!
!      ! Solve the linear Dyson equation, with given self-energy, using
!      ! DLR in imaginary time
!      !
!      ! Input:
!      !
!      ! rank      - rank of DLR (# basis functions)
!      ! it2cf     - imaginary time grid values -> DLR coefficients
!      !               transform matrix in lapack LU storage format
!      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
!      !               format
!      ! phi       - tensor taking a set of DLR coefficients to matrix of
!      !               convolution by corresponding DLR expansion
!      ! g0        - DLR coefficients of g0, produced by subroutine
!      !               g0init_it
!      ! g0mat     - matrix of convolution by G0, produced by subroutine
!      !               g0init_it
!      ! sig       - values of self-energy Sigma at DLR imaginary time
!      !               nodes.
!      !
!      ! Output:
!      !
!      ! g         - DLR coefficients of solution of linear Dyson
!      !               equation 
!
!      implicit none
!      integer rank,it2cfpiv(rank)
!      real *8 it2cf(rank,rank)
!      real *8 g0(rank),g0mat(rank,rank),sig(rank),g(rank)
!      real *8 phi(rank*rank,rank)
!
!      integer j,info1
!      integer, allocatable :: ipiv(:)
!      real *8 one
!      real *8, allocatable :: sigc(:),sysmat(:,:)
!
!      one = 1.0d0
!
!      allocate(sigc(rank),sysmat(rank,rank),ipiv(rank))
!
!      ! DLR coefficient representation of Sigma
!      
!      call dlr_expnd(rank,it2cf,it2cfpiv,sig,sigc)
!
!      ! Matrix of convolution by Sigma
!
!      call dlr_convmat(rank,phi,it2cf,it2cfpiv,sigc,sysmat)
!
!      ! System matrix for linear Dyson equation
!
!      sysmat = -matmul(g0mat,sysmat)
!
!      do j=1,rank
!        sysmat(j,j) = one + sysmat(j,j)
!      enddo
!
!      ! Solve linear equation by LU factorization + backsolve
!
!      call dgetrf(rank,rank,sysmat,rank,ipiv,info1)
!
!      g = g0
!
!      call dgetrs('N',rank,1,sysmat,rank,ipiv,g,rank,info1)
!
!
!      end subroutine dyson_it_lin
!
!
!      subroutine g0init_it(beta,rank,dlrit,it2cf,it2cfpiv,phi,mu,g0,&
!          g0mat)
!
!      ! Get DLR coefficient representation of G0, and matrix of
!      ! convolution by G0, for use in Dyson solver
!      !
!      ! Input:
!      !
!      ! beta      - inverse temperature
!      ! rank      - rank of DLR (# basis functions)
!      ! dlrit     - selected imaginary time nodes (tau points)
!      ! it2cf     - imaginary time grid values -> DLR coefficients
!      !               transform matrix in lapack LU storage format
!      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
!      !               format
!      ! phi       - tensor taking a set of DLR coefficients to matrix of
!      !               convolution by corresponding DLR expansion
!      ! mu        - single particle energy parameter
!      !
!      ! Output:
!      !
!      ! g0        - DLR coefficients of G0
!      ! g0mat     - matrix of convolution by G0
!
!      integer rank,it2cfpiv(rank)
!      real *8 beta,dlrit(rank),it2cf(rank,rank),phi(rank*rank,rank),mu
!      real *8 g0(rank),g0mat(rank,rank)
!
!      integer i
!      real *8, external :: kfunf_rel
!
!      ! Imaginary time grid representation of G0
!
!      do i=1,rank
!
!        g0(i) = -kfunf_rel(dlrit(i),beta*mu)
!
!      enddo
!
!      ! DLR coefficient representation of G0
!
!      call dlr_expnd(rank,it2cf,it2cfpiv,g0,g0)
!
!      ! Matrix of convolution by G0
!
!      call dlr_convmat(rank,phi,it2cf,it2cfpiv,g0,g0mat)
!
!      end subroutine g0init_it

