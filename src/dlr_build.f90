      !
      !
      ! This file contains subroutines used to obtain the DLR basis
      ! functions
      !
      !





      !> Set parameters for composite Chebyshev fine grid
      !!
      !! @param[in]   lambda  dimensionless cutoff parameter
      !! @param[out]  p       Chebyshev degree in each subinterval
      !! @param[out]  npt     number of subintervals on [0,1/2] in
      !!                        imaginary time (total number of
      !!                        subintervals on [0,1] is 2*npt)
      !! @param[out]  npo     number of subintervals on [0,lambda] in
      !!                        real frequency (total number of
      !!                        subintervals on [-lambda,lambda is
      !!                        2*npo)
      !! @param[out]  nt      number of imaginary time fine grid points
      !!                        (=2*p*npt)
      !! @param[out]  no      number of Matsubara frequency fine grid
      !!                        points (=2*p*npo)

      subroutine ccfine_init(lambda,p,npt,npo,nt,no)

      implicit none
      integer p,npt,npo,nt,no
      real *8 lambda

      ! Chebyshev degree of panels

      p = 24

      
      ! Number of panels in imaginary time and real frequency
      !
      ! This choice has been made empirically to give fine
      ! discretization of Lehmann kernel to double machine precision
      ! accuracy
      
      npt = max(ceiling(log(lambda)/log(2.0d0))-2,1)
      npo = max(ceiling(log(lambda)/log(2.0d0)),1)


      ! Total number of fine grid points

      nt = 2*p*npt
      no = 2*p*npo

      end subroutine ccfine_init





      !> Get composite Chebyshev fine grids in tau and omega
      !!
      !! Note: for accurate discretization of Lehmann kernel, use
      !! parameters produced by subroutine gridparams.
      !!
      !! @param[in]   lambda  dimensionless cutoff parameter
      !! @param[in]   p       Chebyshev degree in each subinterval
      !! @param[in]   npt     number of subintervals on [0,1/2] in
      !!                        imaginary time (total number of
      !!                        subintervals on [0,1] is 2*npt)
      !! @param[in]   npo     number of subintervals on [0,lambda] in
      !!                        real frequency (total number of
      !!                        subintervals on [-lambda,lambda] is
      !!                        2*npo)
      !! @param[out]  t       imaginary time fine grid points in
      !!                        relative format
      !! @param[out]  om      real frequency fine grid points

      subroutine ccfine(lambda,p,npt,npo,t,om)

      implicit none
      integer p,npt,npo
      real *8 lambda,t(2*npt*p),om(2*npo*p)

      integer nt,no,i
      real *8 one,a,b
      real *8, allocatable :: xc(:),wc(:),pbpt(:),pbpo(:)


      one = 1.0d0


      ! Total number of fine grid points

      nt = 2*npt*p
      no = 2*npo*p


      ! Get Chebyshev nodes

      allocate(xc(p),wc(p))
      
      call barychebinit(p,xc,wc)


      ! Imaginary time panel break points

      allocate(pbpt(2*npt+1))

      pbpt(1) = 0*one
      do i=1,npt
        pbpt(i+1) = one/2**(npt-i+1)
      enddo


      ! Imaginary time grid points on (0,1/2)

      do i=1,npt
        a = pbpt(i)
        b = pbpt(i+1)
        t((i-1)*p+1:i*p) = a + (b-a)*(xc+one)/2
      enddo


      ! Imaginary time grid points on (1/2,1) in relative format

      t(nt/2+1:nt) = -t(nt/2:1:-1)



      ! Real frequency panel break points

      allocate(pbpo(2*npo+1))

      pbpo(npo+1) = 0*one
      do i=1,npo
        pbpo(npo+i+1) = lambda/2**(npo-i)
      enddo

      pbpo(1:npo) = -pbpo(2*npo+1:npo+2:-1)


      ! Real frequency grid points

      do i=1,2*npo
        a = pbpo(i)
        b = pbpo(i+1)
        om((i-1)*p+1:i*p) = a + (b-a)*(xc+one)/2
      enddo

      end subroutine ccfine





      !> Get fine discretization of Lehmann kernel
      !!
      !! Note: use parameters produced by subroutine ccfine.
      !!
      !! @param[in]   lambda  dimensionless cutoff parameter
      !! @param[in]   p       Chebyshev degree in each subinterval
      !! @param[in]   npt     number of subintervals on [0,1/2] in
      !!                        imaginary time (total number of
      !!                        subintervals on [0,1] is 2*npt)
      !! @param[in]   npo     number of subintervals on [0,lambda] in
      !!                        real frequency (total number of
      !!                        subintervals on [-lambda,lambda] is
      !!                        2*npo)
      !! @param[in]   t       imaginary time fine grid points in
      !!                        relative format
      !! @param[in]   om      real frequency fine grid points
      !! @param[out]  kmat    matrix of Lehmann kernel K(tau,omega)
      !!                        evaluated on fine grid
      !! @param[out]  err     Error of composite Chebyshev interpolant
      !!                        of K(tau,omega). err(1) is ~= max
      !!                        relative L^inf error in tau over all
      !!                        omega in fine grid. err(2) is ~= max
      !!                        L^inf error in omega over all tau in
      !!                        fine grid.

      subroutine dlr_kfine(lambda,p,npt,npo,t,om,kmat,err)

      implicit none
      integer p,npt,npo
      real *8 lambda,t(2*p*npt),om(2*p*npo),kmat(2*p*npt,2*p*npo),err(2)

      integer nt,no,i,j,k
      real *8 one,a,b,start,finish,xx,ktrue,ktest,errtmp
      real *8, allocatable :: xc(:),wc(:),pbpt(:),pbpo(:)
      real *8, allocatable :: xc2(:),wc2(:)
      real *8, external :: kfunf

      one = 1.0d0

      nt = 2*npt*p
      no = 2*npo*p

      ! Chebyshev nodes and interpolation weights

      allocate(xc(p),wc(p))
      
      call barychebinit(p,xc,wc)


      ! Evaluate Lehmann kernel K(tau,omega) on fine grid for tau in
      ! (0,1/2)

      do j=1,no
        do i=1,nt/2

          kmat(i,j) = kfunf(t(i),om(j))

        enddo
      enddo


      ! Get entries for tau in (1/2,1)

      kmat(nt/2+1:nt,1:no) = kmat(nt/2:1:-1,no:1:-1)



      ! Imaginary time panel break points

      allocate(pbpt(2*npt+1))

      pbpt(1) = 0*one
      do i=1,npt
        pbpt(i+1) = one/2**(npt-i+1)
      enddo


      ! Real frequency panel break points

      allocate(pbpo(2*npo+1))

      pbpo(npo+1) = 0*one
      do i=1,npo
        pbpo(npo+i+1) = lambda/2**(npo-i)
      enddo

      pbpo(1:npo) = -pbpo(2*npo+1:npo+2:-1)


      ! Check accuracy of Cheb interpolant on each panel in tau
      ! for fixed omega, and each panel in omega for fixed tau, by
      ! comparing with K(tau,omega) on Cheb grid of 2*p nodes

      allocate(xc2(2*p),wc2(2*p))

      call barychebinit(2*p,xc2,wc2)

      err(:) = 0.0d0

      do j=1,no

        errtmp = 0.0d0

        do i=1,npt
          
          a = pbpt(i)
          b = pbpt(i+1)

          do k=1,2*p

            xx = a+(b-a)*(xc2(k)+one)/2
            
            ktrue = kfunf(xx,om(j))

            call barycheb(p,xc,wc,xc2(k),kmat((i-1)*p+1:i*p,j),ktest)

            errtmp = max(errtmp,abs(ktrue-ktest))

          enddo
        enddo

        err(1) = max(err(1),errtmp/maxval(kmat(:,j)))

      enddo


      do j=1,nt/2

        errtmp = 0.0d0

        do i=1,2*npo
          
          a = pbpo(i)
          b = pbpo(i+1)

          do k=1,2*p

            xx = a+(b-a)*(xc2(k)+one)/2
            
            ktrue = kfunf(t(j),xx)

            call barycheb(p,xc,wc,xc2(k),kmat(j,(i-1)*p+1:i*p),ktest)

            errtmp = max(errtmp,abs(ktrue-ktest))

          enddo
        enddo

        err(2) = max(err(2),errtmp/maxval(kmat(j,:)))

      enddo

      end subroutine dlr_kfine




      !> Get DLR frequency nodes
      !!
      !! @param[in]     lambda  dimensionless cutoff parameter
      !! @param[in]     eps     DLR error tolerance
      !! @param[in]     nt      number of imaginary time fine grid points
      !! @param[in]     no      number of Matsubara frequency fine grid
      !!                          points
      !! @param[in]     om      Matsubara frequency fine grid points
      !! @param[in]     kmat    kernel K(tau,omega), sampled at fine grid
      !!                          points
      !! @param[in,out] r       On input, maximum possible number of DLR
      !!                          basis functions, defining input size
      !!                          of various arrays; on output, number
      !!                          of DLR basis functions.
      !! @param[out]    dlrrf   DLR frequency nodes
      !! @param[out]    oidx    column indices of kmat corresponding to
      !!                          DLR frequency nodes

      subroutine dlr_rf(lambda,eps,nt,no,om,kmat,r,dlrrf,oidx)

      implicit none
      integer nt,no,r,oidx(r)
      real *8 lambda,eps,om(no),kmat(nt,no),dlrrf(r)

      integer, allocatable :: list(:)
      real *8, allocatable :: tmp(:,:),work(:)

      ! Select real frequency nodes by pivoted QR on columns of kmat

      allocate(tmp(nt,no),list(no),work(max(nt,no)))

      tmp = kmat

      call iddp_qrpiv(eps,nt,no,tmp,r,list,work)


      ! Rearrange indices to get selected frequency point indices

      call ind_rearrange(no,r,list)


      ! Extract selected frequencies

      oidx(1:r) = list(1:r)
      dlrrf(1:r) = om(oidx(1:r))
          
      end subroutine dlr_rf



