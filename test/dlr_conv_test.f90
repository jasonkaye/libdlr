      program dlr_conv_test

      ! Test convolution of two DLR expansions
      
      implicit none
      integer ntst
      real *8 lambda,eps,beta

      ! --- Input parameters ---

      lambda = 1000 ! Frequency cutoff
      eps = 1.0d-14 ! Desired accuracy
      ntst = 10000 ! # test points to check representation of G
      beta = 1000 ! Inverse temp: controls support of rho


      ! --- Call main test subroutine ---

      call dlr_conv_test_main(lambda,eps,ntst,beta)


      end program dlr_conv_test


      subroutine dlr_conv_test_main(lambda,eps,ntst,beta)

      ! Main driver routine for test

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,j,r
      integer, allocatable :: it2cfp(:)
      real *8 one,gtrue,gtest,errl2,errlinf,gmax,gl2
      real *8, allocatable :: ttst(:),it2cf(:,:),dlrit(:),dlrrf(:)
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


      ! Build DLR basis, grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_buildit(lambda,eps,r,dlrrf,dlrit)

      ! Get imaginary time values -> DLR coefficients transform matrix in LU form

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf(r,dlrrf,dlrit,it2cf,it2cfp)


      write(6,*) ''
      write(6,*) '-------------------- DLR basis --------------------'
      write(6,*) ''
      write(6,*) 'DLR rank                          = ',r


      ! --- Sample two Green's function and get DLR expansions ---


      ! Sample G1 and G2 at DLR grid points

      allocate(g1(r),g2(r),g3(r))

      do i=1,r

        call gfun1(beta,dlrit(i),g1(i))
        call gfun2(beta,dlrit(i),g2(i))

      enddo


      ! --- Compute convolution and measure error ---

      ! Get convolution tensor

      allocate(phi(r*r,r))

      call dlr_convtens(beta,r,dlrrf,dlrit,it2cf,it2cfp,phi)


      ! Form matrix of convolution by G1

      allocate(gmat(r,r))

      call dlr_convmat(r,it2cf,it2cfp,phi,g1,gmat)


      ! Apply matrix to obtain convolution G3

      g3 = matmul(gmat,g2)

      
      ! Get DLR coefficients of g3

      call dlr_expnd(r,it2cf,it2cfp,g3,g3)


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

        call dlr_eval(r,dlrrf,g3,ttst(i),gtest)

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
      real *8, external :: kfunf_rel

      real *8 a1

      a1 = 0.804d0

      g = kfunf_rel(t,beta*a1)

      end subroutine gfun1

      subroutine gfun2(beta,t,g)

      ! Evaluate Green's function corresponding to
      ! delta function spectral density 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a2

      a2 = -0.443d0
      !a2 = 0.804d0

      g = kfunf_rel(t,beta*a2)

      end subroutine gfun2

      subroutine gfun3(beta,t,g)

      ! Evaluate Green's function corresponding to convolution of two 
      ! delta function spectral densities 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a1,a2

      a1 = 0.804d0
      a2 = -0.443d0

      g = (kfunf_rel(t,beta*a2)-kfunf_rel(t,beta*a1))/(a1-a2)
      
      !g = beta*(t-kfunf_rel(beta*1.0d0,a1))*kfunf_rel(t,beta*a1)

      end subroutine gfun3
