      !
      !
      ! This file contains subroutines for imaginary time convolution of
      ! Green's functions in the discrete Lehmann representation
      !
      !
      
      !> Get tensor taking a set of DLR coefficients to the matrix of
      !! convolution by the corresponding Green's function.
      !!
      !! To obtain the matrix C of convolution by a Green's function G,
      !! use the output of this subroutine with the subroutine
      !! dlr_convtens.
      !!
      !! @param[in]   beta      inverse temperature
      !! @param[in]   rank      number of DLR basis functions
      !! @param[in]   dlrrf     DLR frequency nodes
      !! @param[in]   dlrit     DLR imaginary time nodes
      !! @param[in]   it2cf     imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU factors
      !! @param[in]   it2cfpiv  imaginary time grid values ->
      !!                          DLR coefficients transform matrix,
      !!                          stored in LAPACK LU factored format;
      !!                          LU pivots
      !! @param[out]  phi       tensor taking DLR coefficients of g to
      !!                          matrix C of convolution by g. C takes
      !!                          DLR coefficients of a function f ->
      !!                          DLR grid values of the convolution
      !!                          g * f.

      subroutine dlr_convtens(beta,rank,dlrrf,dlrit,it2cf,it2cfpiv,phi)

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 beta,dlrrf(rank),dlrit(rank),it2cf(rank,rank)
      real *8 phi(rank*rank,rank)
      real *8, external :: kfun

      integer j,k,l,ier,maxrec,numint,info
      real *8 one,rint1,rint2
      real *8, allocatable :: phitmp(:,:,:),phitmp2(:,:)
      real *8, external :: kfunf,kfunf_rel

      one = 1.0d0

      allocate(phitmp(rank,rank,rank),phitmp2(rank,rank*rank))

      do l=1,rank
        do k=1,rank
          do j=1,rank

            if (k.ne.l) then

              phitmp(j,k,l) = (kfunf_rel(dlrit(j),dlrrf(l)) -&
                kfunf_rel(dlrit(j),dlrrf(k)))/(dlrrf(k)-dlrrf(l))

            else

              if (dlrit(j).gt.0.0d0) then

                phitmp(j,k,l) = (dlrit(j)-kfunf(1.0d0,dlrrf(k)))*&
                  kfunf_rel(dlrit(j),dlrrf(k))

              else

                phitmp(j,k,l) = (dlrit(j)+kfunf(0.0d0,dlrrf(k)))*&
                  kfunf_rel(dlrit(j),dlrrf(k))

              endif
            endif

          enddo
        enddo
      enddo



      do l=1,rank
        do k=1,rank
          do j=1,rank
            phitmp2(k,(l-1)*rank+j) = phitmp(j,k,l)
          enddo
        enddo
      enddo
            
      call dgetrs('T',rank,rank*rank,it2cf,rank,it2cfpiv,phitmp2,rank,&
        info)

      do l=1,rank
        do k=1,rank
          do j=1,rank
            phitmp(j,k,l) = phitmp2(k,(l-1)*rank+j)
          enddo
        enddo
      enddo


      do l=1,rank
        do k=1,rank
          do j=1,rank
            phi((k-1)*rank+j,l) = phitmp(j,k,l)
          enddo
        enddo
      enddo


      phi = beta*phi

      end subroutine dlr_convtens





      !> Get the matrix of convolution by a Green's function.
      !!
      !! The output of this subroutine is the matrix gmat of convolution
      !! by G. gmat is applied to a vector of imaginary time grid
      !! values of a Green's function F, and returns the values of
      !! G * F at the imaginary time grid points.
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
      !! @param[in]   phi       tensor produced by subroutine
      !!                          dlr_convtens taking DLR coefficients
      !!                          of g to matrix of convolution by g.
      !! @param[in]   g         values of Green's function at imaginary
      !!                          time grid points
      !! @param[out]  gmat      matrix of convolution by g

      subroutine dlr_convmat(rank,it2cf,it2cfpiv,phi,g,gmat)

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 phi(rank*rank,rank),it2cf(rank,rank),g(rank)
      real *8 gmat(rank,rank)

      integer info
      real *8, allocatable :: gc(:)

      ! Get DLR coefficients of G

      allocate(gc(rank))

      call dlr_expnd(rank,it2cf,it2cfpiv,g,gc)

      ! Get convolution matrix taking coefficients -> values

      call dgemv('N',rank*rank,rank,1.0d0,phi,rank*rank,gc,1,0.0d0,&
        gmat,1)


      end subroutine dlr_convmat





      !> Get tensor taking a set of DLR coefficients to the matrix of
      !! convolution by the corresponding Green's function, acting from
      !! DLR coefficients to DLR grid values.
      !!
      !! Given the DLR coefficients of a Green's function g, applying
      !! this tensor produces a matrix A of convolution by g which takes
      !! the DLR coefficients of a function f to the DLR grid values of
      !! the convolution f * g. Since we typically want a convolution
      !! matrix operating on grid values, we typically use the subroutine
      !! dlr_convtens instead of this one. Alternatively,
      !! one can use the tensor produced by this subroutine with the
      !! subroutine dlr_convmat_vcc to produce a convolution matrix
      !! operating on grid values.
      !!
      !! @param[in]   beta    inverse temperature
      !! @param[in]   rank    number of DLR basis functions
      !! @param[in]   dlrrf   DLR frequency nodes
      !! @param[in]   dlrit   DLR imaginary time nodes
      !! @param[out]  phivcc  tensor taking DLR coefficients of g to
      !!                        matrix C of convolution by g. C takes
      !!                        DLR coefficients of a function f -> DLR
      !!                        grid values of the convolution g * f.
      !!                        "vcc" indicates that for the tensor
      !!                        phi_{ijk}, i indexes grid values, and
      !!                        j,k index DLR coefficients.
      
      subroutine dlr_convtens_vcc(beta,rank,dlrrf,dlrit,phivcc)

      implicit none
      integer rank
      real *8 beta,dlrrf(rank),dlrit(rank)
      real *8 phivcc(rank*rank,rank)
      real *8, external :: kfun

      integer j,k,l,ier,maxrec,numint
      real *8 one,rint1,rint2
      real *8, external :: kfunf,kfunf_rel

      one = 1.0d0

      do l=1,rank
        do k=1,rank
          do j=1,rank

            if (k.ne.l) then

              phivcc((k-1)*rank+j,l) = (kfunf_rel(dlrit(j),dlrrf(l)) -&
                kfunf_rel(dlrit(j),dlrrf(k)))/(dlrrf(k)-dlrrf(l))

            else

              if (dlrit(j).gt.0.0d0) then

                phivcc((k-1)*rank+j,l) = (dlrit(j)-kfunf(1.0d0,dlrrf(k)))*&
                  kfunf_rel(dlrit(j),dlrrf(k))

              else

                phivcc((k-1)*rank+j,l) = (dlrit(j)+kfunf(0.0d0,dlrrf(k)))*&
                  kfunf_rel(dlrit(j),dlrrf(k))

              endif
            endif

          enddo
        enddo
      enddo

      phivcc = beta*phivcc

      end subroutine dlr_convtens_vcc





      !> Get the matrix of convolution by a Green's function from the
      !! tensor produced by the subroutine dlr_convtens_vcc
      !!
      !! The output of this subroutine is the matrix gmat of convolution
      !! by G. gmat is applied to a vector of imaginary time grid
      !! values of a Green's function F, and returns the values of
      !! G * F at the imaginary time grid points.
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
      !! @param[in]   phivcc    tensor produced by subroutine
      !!                          dlr_convtens_vcc taking DLR coefficients
      !!                          of g to matrix of convolution by g,
      !!                          acting from DLR coefficients to DLR
      !!                          grid values.
      !! @param[in]   g         values of Green's function at imaginary
      !!                          time grid points
      !! @param[out]  gmat      matrix of convolution by g

      subroutine dlr_convmat_vcc(rank,it2cf,it2cfpiv,phivcc,g,gmat)

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 phivcc(rank*rank,rank),it2cf(rank,rank),g(rank)
      real *8 gmat(rank,rank)

      integer i,j,info
      real *8, allocatable :: gc(:)

      ! Get DLR coefficients of G

      allocate(gc(rank))

      call dlr_expnd(rank,it2cf,it2cfpiv,g,gc)

      ! Get convolution matrix taking coefficients -> values

      call dgemv('N',rank*rank,rank,1.0d0,phivcc,rank*rank,gc,1,0.0d0,&
        gmat,1)

      ! Precompose with matrix taking values -> coefficients

      gmat = transpose(gmat)

      call dgetrs('T',rank,rank,it2cf,rank,it2cfpiv,gmat,rank,info)

      gmat = transpose(gmat)

      end subroutine dlr_convmat_vcc
