      !
      !
      ! This file contains subroutines to work with the discrete Lehmann
      ! representation in imaginary time
      !
      !





      !> Get DLR frequency nodes and imaginary time nodes
      !!
      !! @param[in]     lambda  dimensionless cutoff parameter
      !! @param[in]     eps     DLR error tolerance
      !! @param[in,out] rank    On input, maximum possible number of DLR
      !!                          basis functions, defining input size
      !!                          of various arrays; on output, number
      !!                          of DLR basis functions.
      !! @param[out]    dlrrf   DLR frequency nodes
      !! @param[out]    dlrit   DLR imaginary time nodes

      subroutine dlr_buildit(lambda,eps,rank,dlrrf,dlrit)

      implicit none
      integer rank
      real *8 lambda,eps,dlrrf(rank),dlrit(rank)

      integer p,npt,npo,nt,no
      integer, allocatable :: oidx(:)
      real *8 kerr(2)
      real *8, allocatable :: kmat(:,:),t(:),om(:)


      ! Set parameters for the fine grid based on lambda

      call gridparams(lambda,p,npt,npo,nt,no)


      ! Get fine composite Chebyshev discretization of K(tau,omega)

      allocate(kmat(nt,no),t(nt),om(no))

      call kfine_cc(lambda,p,npt,npo,t,om,kmat,kerr)


      ! Select real frequency points for DLR basis

      rank = 500 ! Upper bound on possible rank

      allocate(oidx(rank))

      call dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx)


      ! Get DLR imaginary time grid

      call dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit)

      end subroutine dlr_buildit





      !> Get DLR imaginary time nodes
      !!
      !! @param[in]  lambda  dimensionless cutoff parameter
      !! @param[in]  nt      number of imaginary time fine grid points
      !! @param[in]  no      number of Matsubara frequency fine grid
      !!                       points
      !! @param[in]  t       imaginary time fine grid points
      !! @param[in]  kmat    kernel K(tau,omega), sampled at fine grid
      !!                       points
      !! @param[in]  rank    number of DLR basis functions
      !! @param[in]  oidx    column indices of kmat corresponding to
      !!                          DLR frequency nodes
      !! @param[out] dlrit   DLR imaginary time nodes

      subroutine dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit)

      implicit none
      integer nt,no,rank,oidx(rank)
      real *8 lambda,t(nt),kmat(nt,no),dlrit(rank)

      integer j,k
      integer, allocatable :: list(:),tidx(:)
      real *8, allocatable :: tmp(:,:),work(:)

      ! --- Select imaginary time nodes by pivoted QR on rows of 
      ! kmat ---

      ! Matrix of selected columns

      allocate(tmp(rank,nt),list(nt),work(nt),tidx(rank))

      do j=1,nt
        do k=1,rank
          tmp(k,j) = kmat(j,oidx(k))
        enddo
      enddo

      ! Pivoted QR

      call iddr_qrpiv(rank,nt,tmp,rank,list,work)

      ! Rearrange indices to get selected imaginary time node indices

      call ind_rearrange(nt,rank,list)

      ! Extract selected imaginary times. To maintain high precision for
      ! extremely large lambda and small eps calculations, if t was
      ! chosen which is close to 1, take the calculated value t*=1-t,
      ! which is known to full relative precision, and store -t*. Then t
      ! can either be recovered as 1+(-t*), resulting in a loss of
      ! relative precision, or we can use the high relative precision
      ! value directly if we have access to a high accuracy close-to-1
      ! evaluator.

      tidx = list(1:rank)

      do j=1,rank
        if (tidx(j).le.nt/2) then
          dlrit(j) = t(tidx(j))
        else
          dlrit(j) = -t(nt-tidx(j)+1)
        endif
      enddo
      
      end subroutine dlr_it





      !> Build transform matrix from DLR coefficients to values of DLR
      !! expansion on imaginary time grid
      !!
      !! To obtain the values of a DLR expansion on the imaginary time
      !! grid, apply the matrix cf2it to the vector of DLR coefficients
      !!
      !! @param[in]  rank   number of DLR basis functions
      !! @param[in]  dlrrf  DLR frequency nodes
      !! @param[in]  dlrit  DLR imaginary time nodes
      !! @param[out] cf2it  DLR coefficients -> imaginary time grid
      !!                      values transform matrix

      subroutine dlr_cf2it(rank,dlrrf,dlrit,cf2it)

      implicit none
      integer rank
      real *8 dlrrf(rank),dlrit(rank),cf2it(rank,rank)

      integer j,k
      real *8, external :: kfunf_rel

      ! Get the matrix K(tau_j,omega_k)

      do k=1,rank
        do j=1,rank
          cf2it(j,k) = kfunf_rel(dlrit(j),dlrrf(k))
        enddo
      enddo

      end subroutine dlr_cf2it





      !> Build transform matrix from values of a Green's function on
      !! imaginary time grid to its DLR coefficients; matrix is stored
      !! in LU factored form
      !!
      !! To obtain the DLR coefficients of a Green's function from its
      !! values on the imaginary time grid, use the dlr_expnd subroutine
      !! with arrays it2cf and it2cfpiv generated by this subroutine.
      !!
      !! @param[in]  rank      number of DLR basis functions
      !! @param[in]  dlrrf     DLR frequency nodes
      !! @param[in]  dlrit     DLR imaginary time nodes
      !! @param[out] it2cf     imaginary time grid values ->
      !!                         DLR coefficients transform matrix,
      !!                         stored in LAPACK LU factored format; LU
      !!                         factors
      !! @param[out] it2cfpiv  imaginary time grid values ->
      !!                         DLR coefficients transform matrix,
      !!                         stored in LAPACK LU factored format; LU
      !!                         pivots

      subroutine dlr_it2cf(rank,dlrrf,dlrit,it2cf,it2cfpiv)

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 dlrrf(rank),dlrit(rank),it2cf(rank,rank)

      integer j,k,info
      real *8, external :: kfunf_rel

      ! Get the matrix K(tau_j,omega_k)

      do k=1,rank
        do j=1,rank
          it2cf(j,k) = kfunf_rel(dlrit(j),dlrrf(k))
        enddo
      enddo

      ! LU factorize

      call dgetrf(rank,rank,it2cf,rank,it2cfpiv,info)

      end subroutine dlr_it2cf





      !> Build transform matrix from values of a Green's function G on
      !! imaginary time grid to values of reflection G(beta-tau) on
      !! imaginary time grid
      !!
      !! To obtain the values of a reflected Green's function on the
      !! imaginary time grid, apply the matrix it2itr to the vector of
      !! values of the Green's function on the imaginary time grid
      !!
      !! @param[in]  rank      number of DLR basis functions
      !! @param[in]  dlrrf     DLR frequency nodes
      !! @param[in]  dlrit     DLR imaginary time nodes
      !! @param[in]  it2cf     imaginary time grid values ->
      !!                         DLR coefficients transform matrix, stored in
      !!                         LAPACK LU factored format; LU factors
      !! @param[in]  it2cfpiv  imaginary time grid values ->
      !!                         DLR coefficients transform matrix, stored in
      !!                         LAPACK LU factored format; LU pivots
      !! @param[out] it2itr    imaginary time grid values -> reflected
      !!                         imaginary time grid values transform
      !!                         matrix

      subroutine dlr_it2itr(rank,dlrrf,dlrit,it2cf,it2cfpiv,it2itr)

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 dlrrf(rank),dlrit(rank),it2cf(rank,rank),it2itr(rank,rank)

      integer i,j,info
      real *8, external :: kfunf_rel

      ! Get matrix taking DLR coefficients to values of DLR expansion at
      ! imaginary time nodes reflected about tau = beta/2.

      do j=1,rank
        do i=1,rank
          it2itr(i,j) = kfunf_rel(-dlrit(i),dlrrf(j))
        enddo
      enddo

      ! Precompose with matrix taking DLR imaginary time grid values ->
      ! DLR coefficients

      it2itr = transpose(it2itr)

      call dgetrs('T',rank,rank,it2cf,rank,it2cfpiv,it2itr,rank,info)

      it2itr = transpose(it2itr)

      end subroutine dlr_it2itr





      !> Get DLR coefficients of a Green's function from its values on the
      !! imaginary time grid
      !!
      !! @param[in]  rank      number of DLR basis functions
      !! @param[in]  it2cf     imaginary time grid values ->
      !!                         DLR coefficients transform matrix, stored in
      !!                         LAPACK LU factored format; LU factors
      !! @param[in]  it2cfpiv  imaginary time grid values ->
      !!                         DLR coefficients transform matrix, stored in
      !!                         LAPACK LU factored format; LU pivots
      !! @param[in]  g         values of Green's function at imaginary
      !!                         time grid points
      !! @param[out] gc        DLR coefficients of Green's function

      subroutine dlr_expnd(rank,it2cf,it2cfpiv,g,gc)
      
      implicit none
      integer rank,it2cfpiv(rank)
      real *8 it2cf(rank,rank),g(rank),gc(rank)

      integer info

      ! Backsolve with imaginary time grid values -> DLR coefficients
      ! transform matrix stored in LU form

      gc = g

      call dgetrs('N',rank,1,it2cf,rank,it2cfpiv,gc,rank,info)

      end subroutine dlr_expnd





      !> Evaluate a DLR expansion at an imaginary time point
      !!
      !! @param[in]  rank   number of DLR basis functions
      !! @param[in]  dlrrf  DLR frequency nodes
      !! @param[in]  gc     DLR coefficients of expansion
      !! @param[in]  t      imaginary time point in relative format
      !! @param[out] gt     value of DLR expansion at t

      subroutine dlr_eval(rank,dlrrf,gc,t,gt)

      implicit none
      integer rank
      real *8 dlrrf(rank),gc(rank),t,gt

      integer i
      real *8 kval
      real *8, external :: kfunf

      gt = 0.0d0
      do i=1,rank

        ! For 0.5<t<1, corresponding to negative t', use symmetry of K
        ! to evaluate basis functions

        if (t.ge.0.0d0) then
          kval = kfunf(t,dlrrf(i))
        else
          kval = kfunf(-t,-dlrrf(i))
        endif

        gt = gt + gc(i)*kval

      enddo

      end subroutine dlr_eval



