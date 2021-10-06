
      !
      !
      ! This file contains subroutines for imaginary time convolution of
      ! Green's functions in the discrete Lehmann representation
      !
      !


      !> Get tensor taking a set of DLR coefficients to the matrix of
      !! convolution by the corresponding Green's function
      !!
      !! To obtain the matrix C of convolution by a Green's function G,
      !! use the output of this subroutine with the subroutine
      !! dlr_convtens.
      !!
      !! @param[in]   beta    inverse temperature
      !! @param[in]   xi      xi=-1 for fermionic case; xi=1 for bosonic
      !!                        case
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   dlrrf   DLR frequency nodes
      !! @param[in]   dlrit   DLR imaginary time nodes
      !! @param[in]   it2cf   imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU factors
      !! @param[in]   it2cfp  imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU pivots
      !! @param[out]  phi     tensor taking DLR coefficients of g to
      !!                        matrix C of convolution by g. C takes
      !!                        DLR coefficients of a function f ->
      !!                        DLR grid values of the convolution
      !!                        g * f.

      subroutine dlr_convtens(beta,xi,r,dlrrf,dlrit,it2cf,it2cfp,phi)

      implicit none
      integer r,xi,it2cfp(r)
      real *8 beta,dlrrf(r),dlrit(r),it2cf(r,r)
      real *8 phi(r*r,r)
      real *8, external :: kfun

      integer j,k,l,ier,maxrec,numint,info
      real *8 one,rint1,rint2
      real *8, allocatable :: phitmp(:,:,:),phitmp2(:,:)
      real *8, external :: kfunf,kfunf_rel,expfun

      one = 1.0d0

      ! Get tensor taking DLR coefficients to matrix of convolution from
      ! DLR coefficients -> imaginary time grid values

      allocate(phitmp(r,r,r),phitmp2(r,r*r))

      do l=1,r
        do k=1,r
          do j=1,r

            if (k.ne.l) then

              phitmp(j,k,l) =&
                (kfunf_rel(dlrit(j),dlrrf(l))*expfun(dlrrf(k),xi) -&
                kfunf_rel(dlrit(j),dlrrf(k))*expfun(dlrrf(l),xi))/&
                (dlrrf(k)-dlrrf(l))

            else

              if (dlrit(j).gt.0.0d0) then

                phitmp(j,k,l) = (dlrit(j)*expfun(dlrrf(k),xi)+&
                  xi*kfunf(1.0d0,dlrrf(k)))*kfunf_rel(dlrit(j),dlrrf(k))

              else

                phitmp(j,k,l) = (dlrit(j)*expfun(dlrrf(k),xi)+&
                  kfunf(0.0d0,dlrrf(k)))*kfunf_rel(dlrit(j),dlrrf(k))

              endif
            endif

          enddo
        enddo
      enddo


      ! Precompose with imaginary time grid values -> DLR coefficients
      ! transformation to obtain final tensor

      do l=1,r
        do k=1,r
          do j=1,r
            phitmp2(k,(l-1)*r+j) = phitmp(j,k,l)
          enddo
        enddo
      enddo
            
      call dgetrs('T',r,r*r,it2cf,r,it2cfp,phitmp2,r,info)

      do l=1,r
        do k=1,r
          do j=1,r
            phitmp(j,k,l) = phitmp2(k,(l-1)*r+j)
          enddo
        enddo
      enddo


      ! Reshape

      do l=1,r
        do k=1,r
          do j=1,r
            phi((k-1)*r+j,l) = phitmp(j,k,l)
          enddo
        enddo
      enddo


      ! Scale by beta for convolutions on [0,beta]

      phi = beta*phi

      end subroutine dlr_convtens





      !> Get the matrix of convolution by a Green's function.
      !!
      !! The output of this subroutine is the matrix gmat of convolution
      !! by G. gmat is applied to a vector of imaginary time grid
      !! values of a Green's function F, and returns the values of
      !! G * F at the imaginary time grid points.
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
      !! @param[in]   phi     tensor produced by subroutine
      !!                        dlr_convtens taking DLR coefficients
      !!                        of g to matrix of convolution by g.
      !! @param[in]   g       values of Green's function at imaginary
      !!                        time grid points
      !! @param[out]  gmat    matrix of convolution by g

      subroutine dlr_convmat(r,it2cf,it2cfp,phi,g,gmat)

      implicit none
      integer r,it2cfp(r)
      real *8 phi(r*r,r),it2cf(r,r),g(r)
      real *8 gmat(r,r)

      integer info
      real *8, allocatable :: gc(:)

      ! Get DLR coefficients of G

      allocate(gc(r))

      call dlr_it_expnd(r,it2cf,it2cfp,g,gc)


      ! Get convolution matrix

      call dgemv('N',r*r,r,1.0d0,phi,r*r,gc,1,0.0d0,gmat,1)

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
      !! @param[in]   xi      xi=-1 for fermionic case; xi=1 for bosonic
      !!                        case
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   dlrrf   DLR frequency nodes
      !! @param[in]   dlrit   DLR imaginary time nodes
      !! @param[out]  phivcc  tensor taking DLR coefficients of g to
      !!                        matrix C of convolution by g. C takes
      !!                        DLR coefficients of a function f -> DLR
      !!                        grid values of the convolution g * f.
      !!                        "vcc" indicates that for the tensor
      !!                        phi_{ijk}, i indexes grid values, and
      !!                        j,k index DLR coefficients.
      
      subroutine dlr_convtens_vcc(beta,xi,r,dlrrf,dlrit,phivcc)

      implicit none
      integer r,xi
      real *8 beta,dlrrf(r),dlrit(r)
      real *8 phivcc(r*r,r)
      real *8, external :: kfun

      integer j,k,l,ier,maxrec,numint
      real *8 one,rint1,rint2
      real *8, external :: kfunf,kfunf_rel,expfun

      one = 1.0d0

      ! Get tensor taking DLR coefficients to matrix of convolution from
      ! DLR coefficients -> imaginary time grid values

      do l=1,r
        do k=1,r
          do j=1,r

            if (k.ne.l) then

              phivcc((k-1)*r+j,l) = &
                (kfunf_rel(dlrit(j),dlrrf(l))*expfun(dlrrf(k),xi) -&
                kfunf_rel(dlrit(j),dlrrf(k))*expfun(dlrrf(l),xi))/&
                (dlrrf(k)-dlrrf(l))

            else

              if (dlrit(j).gt.0.0d0) then

                phivcc((k-1)*r+j,l) = (dlrit(j)*expfun(dlrrf(k),xi)+&
                  xi*kfunf(1.0d0,dlrrf(k)))*kfunf_rel(dlrit(j),dlrrf(k))

              else

                phivcc((k-1)*r+j,l) = (dlrit(j)*expfun(dlrrf(k),xi)+&
                  kfunf(0.0d0,dlrrf(k)))*kfunf_rel(dlrit(j),dlrrf(k))

              endif
            endif

          enddo
        enddo
      enddo

      
      ! Scale by beta for convolutions on [0,beta]

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
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   it2cf   imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU factors
      !! @param[in]   it2cfp  imaginary time grid values ->
      !!                        DLR coefficients transform matrix,
      !!                        stored in LAPACK LU factored format;
      !!                        LU pivots
      !! @param[in]   phivcc  tensor produced by subroutine
      !!                        dlr_convtens_vcc taking DLR coefficients
      !!                        of g to matrix of convolution by g,
      !!                        acting from DLR coefficients to DLR
      !!                        grid values.
      !! @param[in]   g       values of Green's function at imaginary
      !!                        time grid points
      !! @param[out]  gmat    matrix of convolution by g

      subroutine dlr_convmat_vcc(r,it2cf,it2cfp,phivcc,g,gmat)

      implicit none
      integer r,it2cfp(r)
      real *8 phivcc(r*r,r),it2cf(r,r),g(r)
      real *8 gmat(r,r)

      integer i,j,info
      real *8, allocatable :: gc(:)

      ! Get DLR coefficients of G

      allocate(gc(r))

      call dlr_it_expnd(r,it2cf,it2cfp,g,gc)


      ! Get convolution matrix taking DLR coefficients -> imaginary time
      ! grid values

      call dgemv('N',r*r,r,1.0d0,phivcc,r*r,gc,1,0.0d0,&
        gmat,1)


      ! Precompose with matrix taking imaginary time grid values -> DLR
      ! coefficients

      gmat = transpose(gmat)

      call dgetrs('T',r,r,it2cf,r,it2cfp,gmat,r,info)

      gmat = transpose(gmat)

      end subroutine dlr_convmat_vcc


      !> Evaluate the function \f$f(\omega) = (1-\xi \cdot \exp(-\omega))/(1+\exp(-\omega))\f$,
      !! for \f$\xi = \pm 1\f$
      !!
      !! @param[in]   om      frequency \f$\omega\f$
      !! @param[in]   xi      xi=-1 for fermionic case; xi=1 for bosonic
      !!                        case
      !!
      !! @param[out]  expfun  value of the function \f$f(\omega)\f$

      function expfun(om,xi)

      ! Evaluate the function f(om) = (1-xi*exp(-om))/(1+exp(-om)), for
      ! xi = +-1

      implicit none
      real *8 expfun
      real *8, intent(in) :: om
      integer, intent(in) :: xi

      if (om.ge.0.0d0) then
        expfun = (1.0d0-xi*exp(-om))/(1.0d0+exp(-om))
      else
        expfun = (exp(om)-1.0d0*xi)/(exp(om)+1.0d0)
      endif

      end function expfun
