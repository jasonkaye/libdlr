      program dlr_ha_test

      ! Test discrete Lehmann representation using Green's function
      ! generated from Lehmann representation with density which is a
      ! sum of two delta functions. Recover DLR coefficients from
      ! samples of Green's function at DLR grid points, and then measure
      ! the error of the resulting expansion on a test grid.
      
      implicit none
      integer ntst
      real *8 lambda,eps,beta

      ! --- Input parameters ---

      lambda = 1000 ! Frequency cutoff
      eps = 1.0d-14 ! Desired accuracy
      ntst = 10000 ! # test points to check representation of G
      beta = 1000 ! Inverse temp: controls support of rho


      ! --- Call main test subroutine ---

      call dlr_ha_test_main(lambda,eps,ntst,beta)


      end program dlr_ha_test


      subroutine dlr_ha_test_main(lambda,eps,ntst,beta)

      ! Main driver routine for test of DLR basis on Green's function
      ! with two delta function density

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,j,rank
      integer, allocatable :: it2cfpiv(:)
      real *8 one,gtrue,gtest,errl2,errlinf,gmax,gl2
      real *8, allocatable :: ttst(:),it2cf(:,:),dlrit(:),dlrrf(:)
      real *8, allocatable :: g(:),gc(:)

      one = 1.0d0

      write(6,*) ''
      write(6,*) '---------------- Input parameters ----------------'
      write(6,*) ''
      write(6,*) 'Cutoff lambda            = ',lambda
      write(6,*) 'Error tolerance eps      = ',eps
      write(6,*) 'Inverse temp beta        = ',beta
      write(6,*) '# test points            = ',ntst


      ! Build DLR basis, grid

      rank = 500 ! Upper bound on rank

      allocate(dlrrf(rank),dlrit(rank))

      call dlr_buildit(lambda,eps,rank,dlrrf,dlrit)


      write(6,*) ''
      write(6,*) '-------------------- DLR basis --------------------'
      write(6,*) ''
      write(6,*) 'DLR rank                          = ',rank


      ! Get imaginary time values -> DLR coefficients transform matrix in LU form

      allocate(it2cf(rank,rank),it2cfpiv(rank))

      call dlr_it2cf(rank,dlrrf,dlrit,it2cf,it2cfpiv)


      ! Sample G(tau) at DLR grid points

      allocate(g(rank),gc(rank))

      do i=1,rank

        call gfun(beta,dlrit(i),g(i))

      enddo


      ! Compute coefficients of DLR expansion from samples

      call dlr_expnd(rank,it2cf,it2cfpiv,g,gc)


      ! Get test points at which to measure error of Green's function;
      ! test points given in relative format

      allocate(ttst(ntst))

      call eqpts_rel(ntst,ttst)

      errlinf = 0*one
      errl2 = 0*one
      gmax = 0*one
      gl2 = 0*one

      do i=1,ntst

        ! Evaluate Green's function

        call gfun(beta,ttst(i),gtrue)

        ! Evaluate DLR

        call dlr_eval(rank,dlrrf,gc,ttst(i),gtest)

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

      end subroutine dlr_ha_test_main



      subroutine gfun(beta,t,g)

      ! Evaluate Green's function corresponding to
      ! sum-of-delta-functions spectral density 

      implicit none
      real *8 beta,t,g
      real *8, external :: kfunf_rel

      real *8 a1,a2,a3,a4,a5

      a1 = -0.804d0
      a2 = -0.443d0
      a3 =  0.093d0
      a4 =  0.915d0
      a5 =  0.929d0

      g = kfunf_rel(t,beta*a1) + kfunf_rel(t,beta*a2) &
        + kfunf_rel(t,beta*a3) + kfunf_rel(t,beta*a4) &
        + kfunf_rel(t,beta*a5)

      end subroutine gfun
