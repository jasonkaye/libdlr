      !
      !
      ! This file contains subroutines for imaginary time convolution of
      ! Green's functions in the discrete Lehmann representation
      !
      !
      
     

      
      subroutine dlr_convtens(beta,rank,dlrrf,dlrit,phi)

      ! Get tensor phi_{jkl} used to take a set of DLR coefficients to
      ! the matrix A of convolution by the corresponding DLR expansion.
      !
      ! A is applied to a vector of DLR coefficients, and returns the
      ! convolution at the DLR imaginary time nodes.
      !
      ! Given phi, the matrix A of convolution by a function with DLR
      ! coefficients rho_l is given by
      !
      ! A_jk = sum_l phi_jkl rho_l.
      !
      ! Input:
      !
      ! beta  - inverse temperature
      ! rank  - rank of DLR (# basis functions)
      ! dlrrf - selected real frequency nodes (omega points)
      ! dlrit - selected imaginary time nodes (tau points)
      !
      ! Output:
      !
      ! phi   - convolution tensor


      implicit none
      integer rank
      real *8 beta,dlrrf(rank),dlrit(rank)
      real *8 phi(rank*rank,rank)
      real *8, external :: kfun

      integer j,k,l,ier,maxrec,numint
      real *8 one,rint1,rint2
      real *8, external :: kfunf,kfunf_rel

      one = 1.0d0

      do l=1,rank
        do k=1,rank
          do j=1,rank

            if (k.ne.l) then

              phi((k-1)*rank+j,l) = (kfunf_rel(dlrit(j),dlrrf(l)) -&
                kfunf_rel(dlrit(j),dlrrf(k)))/(dlrrf(k)-dlrrf(l))

            else

              if (dlrit(j).gt.0.0d0) then

                phi((k-1)*rank+j,l) = (dlrit(j)-kfunf(1.0d0,dlrrf(k)))*&
                  kfunf_rel(dlrit(j),dlrrf(k))

              else

                phi((k-1)*rank+j,l) = (dlrit(j)+kfunf(0.0d0,dlrrf(k)))*&
                  kfunf_rel(dlrit(j),dlrrf(k))

              endif
            endif

          enddo
        enddo
      enddo

      phi = beta*phi

      end subroutine dlr_convtens


      subroutine dlr_convtens2(beta,rank,dlrrf,dlrit,it2cf,it2cfpiv,phi)

      ! Get tensor phi_{jkl} used to take a set of DLR coefficients to
      ! the matrix A of convolution by the corresponding DLR expansion.
      !
      ! A is applied to a vector of values of a function at the DLR
      ! imaginary time nodes, and returns the
      ! convolution at the DLR imaginary time nodes.
      !
      ! Given phi, the matrix A of convolution by a function with DLR
      ! coefficients rho_l is given by
      !
      ! A_jk = sum_l phi_jkl rho_l.
      !
      ! Input:
      !
      ! beta  - inverse temperature
      ! rank  - rank of DLR (# basis functions)
      ! dlrrf - selected real frequency nodes (omega points)
      ! dlrit - selected imaginary time nodes (tau points)
      !
      ! Output:
      !
      ! phi   - convolution tensor


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

      end subroutine dlr_convtens2



!      subroutine dlr_convtens3(beta,rank,dlrrf,dlrit,phi)
!
!      ! Get tensor phi_{jkl} used to take the values of a DLR expansion at
!      ! the DLR imaginary time nodes to the matrix A of convolution by
!      ! the corresponding DLR expansion.
!      !
!      ! A is applied to a vector of values of a function at the DLR
!      ! imaginary time nodes, and returns the
!      ! convolution at the DLR imaginary time nodes.
!      !
!      ! Given phi, the matrix A of convolution by a function with values
!      ! g_l at the DLR imaginary time nodes  is given by
!      !
!      ! A_jk = sum_l phi_jkl g_l.
!      !
!      ! We note that forming this tensor requires computations in
!      ! quadruple precision arithmetic to circumvent a numerical
!      ! instability, but dlrrf and dlrit do not need
!      ! to be computed to quadruple precision, and the tensor phi is
!      ! returned in double precision.
!      !
!      ! Input:
!      !
!      ! beta  - inverse temperature
!      ! rank  - rank of DLR (# basis functions)
!      ! dlrrf - selected real frequency nodes (omega points)
!      ! dlrit - selected imaginary time nodes (tau points)
!      !
!      ! Output:
!      !
!      ! phi   - convolution tensor
!
!
!      implicit none
!      integer rank
!      real *8 beta,dlrrf(rank),dlrit(rank)
!      real *8 phi(rank*rank,rank)
!      real *8, external :: kfun
!
!      integer j,k,l,info
!      integer, allocatable :: ipvt(:)
!      real *16, allocatable :: phitmp(:,:,:),phitmp2(:,:),it2cf(:,:)
!      real *16, allocatable :: qdlrit(:),qdlrrf(:)
!      real *16, external :: qkfunf,qkfunf_rel
!
!
!      allocate(qdlrit(rank),qdlrrf(rank))
!      allocate(phitmp(rank,rank,rank),phitmp2(rank,rank*rank))
!      allocate(it2cf(rank,rank),ipvt(rank))
!
!      qdlrit = dlrit
!      qdlrrf = dlrrf
!
!      do l=1,rank
!        do k=1,rank
!          do j=1,rank
!
!            if (k.ne.l) then
!
!              phitmp(j,k,l) = (qkfunf_rel(qdlrit(j),qdlrrf(l)) -&
!                qkfunf_rel(qdlrit(j),qdlrrf(k)))/(qdlrrf(k)-qdlrrf(l))
!
!            else
!
!              if (dlrit(j).gt.0.0d0) then
!
!                phitmp(j,k,l) = (qdlrit(j)-qkfunf(1.0q0,qdlrrf(k)))*&
!                  qkfunf_rel(qdlrit(j),qdlrrf(k))
!
!              else
!
!                phitmp(j,k,l) = (qdlrit(j)+qkfunf(0.0q0,qdlrrf(k)))*&
!                  qkfunf_rel(qdlrit(j),qdlrrf(k))
!
!              endif
!            endif
!
!          enddo
!        enddo
!      enddo
!
!
!
!      do k=1,rank
!        do j=1,rank
!          it2cf(j,k) = qkfunf_rel(qdlrit(j),qdlrrf(k))
!        enddo
!      enddo
!
!      call qgefa(it2cf,rank,rank,ipvt,info)
!
!
!
!
!      do l=1,rank
!        do k=1,rank
!          do j=1,rank
!            phitmp2(l,(k-1)*rank+j) = phitmp(j,k,l)
!          enddo
!        enddo
!      enddo
!
!      do k=1,rank*rank
!        call qgesl(it2cf,rank,rank,ipvt,phitmp2(:,k),1)
!      enddo
!            
!      do l=1,rank
!        do k=1,rank
!          do j=1,rank
!            phitmp(j,k,l) = phitmp2(l,(k-1)*rank+j)
!          enddo
!        enddo
!      enddo
!
!
!
!
!      do l=1,rank
!        do k=1,rank
!          do j=1,rank
!            phitmp2(k,(l-1)*rank+j) = phitmp(j,k,l)
!          enddo
!        enddo
!      enddo
!
!      do k=1,rank*rank
!        call qgesl(it2cf,rank,rank,ipvt,phitmp2(:,k),1)
!      enddo
!            
!      do l=1,rank
!        do k=1,rank
!          do j=1,rank
!            phitmp(j,k,l) = phitmp2(k,(l-1)*rank+j)
!          enddo
!        enddo
!      enddo
!
!
!
!      do l=1,rank
!        do k=1,rank
!          do j=1,rank
!            phi((k-1)*rank+j,l) = phitmp(j,k,l)
!          enddo
!        enddo
!      enddo
!
!
!      phi = beta*phi
!
!      end subroutine dlr_convtens3




!      subroutine dlr_convmat(rank,phi,it2cf,it2cfpiv,g,gmat)
!
!      ! Get matrix of convolution by a DLR expansion G in the DLR basis
!      ! -- that is, the matrix that this subroutine produces takes the
!      ! DLR coefficient representation of a function f to the DLR
!      ! coefficient representation of the convolution
!      !
!      ! int_0^1 G(t-t') f(t') dt'.
!      !
!      ! Input:
!      !
!      ! rank      - rank of DLR (# basis functions)
!      ! phi       - convolution tensor
!      ! it2cf  - imaginary time grid values -> DLR coefficients
!      !               transform matrix in lapack LU storage format
!      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
!      !               format
!      ! g         - DLR coefficients of a function G
!      !
!      ! Output:
!      !
!      ! gmat      - matrix of convolution by G in the DLR basis
!
!      implicit none
!      integer rank,it2cfpiv(rank)
!      real *8 phi(rank*rank,rank),it2cf(rank,rank),g(rank)
!      real *8 gmat(rank,rank)
!
!      integer i,j,info
!
!      call dgemv('N',rank*rank,rank,1.0d0,phi,rank*rank,g,1,0.0d0,&
!        gmat,1)
!
!      call dgetrs('N',rank,rank,it2cf,rank,it2cfpiv,gmat,rank,info)
!
!      end subroutine dlr_convmat


      subroutine dlr_convmat(rank,phi,it2cf,it2cfpiv,g,gmat)

      ! Get matrix of convolution by a DLR expansion G
      ! -- that is, the matrix that this subroutine produces takes the
      ! DLR imaginary time grid values of a function f to the DLR
      ! imaginary time grid values of the convolution
      !
      ! int_0^1 G(t-t') f(t') dt'.
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! phi       - convolution tensor
      ! it2cf  - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! g         - DLR imaginary time values of a Green's function G
      !
      ! Output:
      !
      ! gmat      - matrix of convolution by G

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 phi(rank*rank,rank),it2cf(rank,rank),g(rank)
      real *8 gmat(rank,rank)

      integer i,j,info
      real *8, allocatable :: gc(:)

      ! Get DLR coefficients of G

      allocate(gc(rank))

      call dlr_expnd(rank,it2cf,it2cfpiv,g,gc)

      ! Get convolution matrix taking coefficients -> values

      call dgemv('N',rank*rank,rank,1.0d0,phi,rank*rank,gc,1,0.0d0,&
        gmat,1)

      ! Precompose with matrix taking values -> coefficients

      gmat = transpose(gmat)

      call dgetrs('T',rank,rank,it2cf,rank,it2cfpiv,gmat,rank,info)

      gmat = transpose(gmat)

      end subroutine dlr_convmat



      subroutine dlr_convmat2(rank,phi,it2cf,it2cfpiv,g,gmat)

      ! Get matrix of convolution by a DLR expansion G
      ! -- that is, the matrix that this subroutine produces takes the
      ! DLR imaginary time grid values of a function f to the DLR
      ! imaginary time grid values of the convolution
      !
      ! int_0^1 G(t-t') f(t') dt'.
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! phi       - convolution tensor
      ! it2cf  - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! g         - DLR imaginary time values of a Green's function G
      !
      ! Output:
      !
      ! gmat      - matrix of convolution by G

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


      end subroutine dlr_convmat2



      subroutine dlr_convmat3(rank,phi,g,gmat)

      ! Get matrix of convolution by a DLR expansion G
      ! -- that is, the matrix that this subroutine produces takes the
      ! DLR imaginary time grid values of a function f to the DLR
      ! imaginary time grid values of the convolution
      !
      ! int_0^1 G(t-t') f(t') dt'.
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! phi       - convolution tensor
      ! it2cf  - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage
      !               format
      ! g         - DLR imaginary time values of a Green's function G
      !
      ! Output:
      !
      ! gmat      - matrix of convolution by G

      implicit none
      integer rank
      real *8 phi(rank*rank,rank),g(rank),gmat(rank,rank)

      integer info


      call dgemv('N',rank*rank,rank,1.0d0,phi,rank*rank,g,1,0.0d0,&
        gmat,1)


      end subroutine dlr_convmat3

