      program dlr_conv_test

      ! Test convolution of two DLR expansions
      
      implicit none
      integer ntst
      real *8 lambda,eps,beta
      character :: fb

      ! --- Input parameters ---

      lambda = 1000 ! Frequency cutoff
      eps = 1.0d-14 ! Desired accuracy
      ntst = 10000 ! # test points to check representation of G
      beta = 1000 ! Inverse temp: controls support of rho
      fb = 'f' ! Fermion or Boson? (this switch isn't working yet)


      ! --- Call main test subroutine ---

      call dlr_conv_test_main(lambda,eps,ntst,beta,fb)


      end program dlr_conv_test


      subroutine dlr_conv_test_main(lambda,eps,ntst,beta,fb)

      ! Main driver routine for test

      implicit none
      integer ntst
      real *8 lambda,eps,beta
      character :: fb

      integer npt,npo,p,nt,no,i,j,rank,info,pg,npg
      integer, allocatable :: ipiv(:),tidx(:),oidx(:)
      real *8 one,gtrue,gtest,errl2,errlinf,kerr(2),gmax,gl2,tabs
      real *8, allocatable :: kmat(:,:),t(:),om(:),ttst(:)
      real *8, allocatable :: it2cf(:,:),dlrit(:),dlrrf(:)
      real *8, allocatable :: g1(:),g2(:),g3(:) 

      real *8, allocatable :: phi(:,:),gmat(:,:)

      one = 1.0d0

      write(6,*) ''
      write(6,*) '---------------- Input parameters ----------------'
      write(6,*) ''
      write(6,*) 'Cutoff lambda            = ',lambda
      write(6,*) 'Error tolerance eps      = ',eps
      write(6,*) 'Inverse temp beta        = ',beta
      write(6,*) '# test points            = ',ntst
      write(6,*) 'Fermion (f) or boson (b) = ',fb


      ! --- Build DLR basis, grid, transform matrix ---

      ! Set parameters for the fine grid based on lambda

      call gridparams(lambda,p,npt,npo,nt,no)

      ! Get fine composite Chebyshev discretization of K(tau,omega)

      allocate(kmat(nt,no),t(nt),om(no))

      call kfine_cc(fb,lambda,p,npt,npo,t,om,kmat,kerr)

      write(6,*) ''
      write(6,*) '-------------- Fine K discretization --------------'
      write(6,*) ''
      write(6,*) '# fine grid pts in tau     = ',nt
      write(6,*) '# fine grid pts in omega   = ',no
      write(6,*) 'Max rel L^inf err in tau   = ',kerr(1)
      write(6,*) 'Max rel L^inf err in omega = ',kerr(2)


      ! Select real frequency points for DLR basis

      rank = 500 ! Upper bound on possible rank

      allocate(dlrrf(rank),oidx(rank))

      call dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx)


      ! Get DLR imaginary time grid

      allocate(dlrit(rank),tidx(rank))

      call dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit,tidx)


      ! Get imaginary time values -> DLR coefficients transform matrix in LU form

      allocate(it2cf(rank,rank),ipiv(rank))

      call dlr_it2cf(nt,no,kmat,rank,oidx,tidx,it2cf,ipiv)



      ! --- Compute actual eps-rank of fine grid K matrix by SVD ---

      write(6,*) ''
      write(6,*) '-------------------- DLR basis --------------------'
      write(6,*) ''
      write(6,*) 'DLR rank                          = ',rank


      ! --- Sample two Green's function and get DLR expansions ---


      ! Sample G1 and G2 at DLR grid points

      allocate(g1(rank),g2(rank),g3(rank))

      do i=1,rank

        call gfun1(beta,dlrit(i),g1(i))
        call gfun2(beta,dlrit(i),g2(i))

      enddo


      ! Compute coefficients of DLR expansions from samples

      call dlr_expnd(rank,it2cf,ipiv,g1)
      call dlr_expnd(rank,it2cf,ipiv,g2)


      ! --- Compute convolution and measure error ---

      ! Get convolution tensor

      allocate(phi(rank*rank,rank))

      call dlr_convtens(rank,dlrrf,dlrit,phi)
      phi = beta*phi


      ! Form matrix of convolution by G1

      allocate(gmat(rank,rank))

      call dlr_convmat(rank,phi,it2cf,ipiv,g1,gmat)


      ! Apply matrix to obtain convolution G3

      g3 = matmul(gmat,g2)

      ! Get test points and compare computed G3 against exact
      ! convolution

      allocate(ttst(ntst))

      call eqpts_rel(ntst,ttst)

      errlinf = 0*one
      errl2 = 0*one
      gmax = 0*one
      gl2 = 0*one

      do i=1,ntst

        ! Evaluate true convolution

        call gfun3(beta,ttst(i),gtrue)

        ! Evaluate DLR

        call dlr_eval(fb,rank,dlrrf,g3,ttst(i),gtest)

        ! Update L^inf and L^2 errors, norms

        errlinf = max(errlinf,abs(gtrue-gtest))
        errl2 = errl2 + (gtrue-gtest)**2

        gmax = max(gmax,abs(gtrue))
        gl2 = gl2 + gtrue**2

      enddo

      errl2 = sqrt((ttst(2)-ttst(1))*errl2)
      gl2 = sqrt((ttst(2)-ttst(1))*gl2)

      write(6,*) ''
      write(6,*) '-------------------- DLR error --------------------'
      write(6,*) ''
      write(6,*) 'Abs L^inf err = ',errlinf
      write(6,*) 'Abs L^2 err   = ',errl2
      write(6,*) 'Rel L^inf err = ',errlinf/gmax
      write(6,*) 'Rel L^2 err   = ',errl2/gl2
      write(6,*) ''

      ! Return failed status if error is not sufficiently small

      if (errlinf.gt.1.0d-13) then
        call exit(1)
      endif

      end subroutine dlr_conv_test_main



      subroutine gfun1(beta,t,g)

      ! Evaluate Green's function corresponding to
      ! delta function spectral density 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf2

      real *8 a1

      a1 = 0.804d0

      g = kfunf2(t,beta*a1)

      end subroutine gfun1

      subroutine gfun2(beta,t,g)

      ! Evaluate Green's function corresponding to
      ! delta function spectral density 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf2

      real *8 a2

      a2 = -0.443d0
      !a2 = 0.804d0

      g = kfunf2(t,beta*a2)

      end subroutine gfun2

      subroutine gfun3(beta,t,g)

      ! Evaluate Green's function corresponding to convolution of two 
      ! delta function spectral densities 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf2

      real *8 a1,a2

      a1 = 0.804d0
      a2 = -0.443d0

      g = (kfunf2(t,beta*a2)-kfunf2(t,beta*a1))/(a1-a2)
      
      !g = beta*(t-kfunf2(beta*1.0d0,a1))*kfunf2(t,beta*a1)

      end subroutine gfun3
