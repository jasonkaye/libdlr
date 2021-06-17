      !
      !
      ! This file contains the core subroutines for working with the
      ! discrete Lehmann representation
      !
      !


      subroutine gridparams(lambda,p,npt,npo,nt,no)

      ! Set parameters for composite Chebyshev fine grid
      !
      ! Input:
      !
      ! lambda  - cutoff parameter
      !
      ! Output:
      !
      ! p       - Chebyshev degree in each subinterval
      ! npt     - # subintervals on [0,1/2] in tau space (# subintervals
      !             on [0,1] is 2*npt)
      ! npo     - # subintervals on [0,lambda] in omega space (#
      !             subintervals on [-lambda,lambda] is 2*npo)
      ! nt      - # fine grid points in tau = 2*npt*p
      ! no      - # fine grid points in omega = 2*npo*p

      implicit none
      integer p,npt,npo,nt,no
      real *8 lambda

      p = 24 ! Chebyshev degree of panels
      
      npt = max(ceiling(log(lambda)/log(2.0d0))-2,1)
      npo = max(ceiling(log(lambda)/log(2.0d0)),1)

      nt = 2*p*npt
      no = 2*p*npo

      end subroutine gridparams


      subroutine kfine_cc(lambda,p,npt,npo,t,om,kmat,err)

      ! Discretization of kernel K(tau,omega) on composite Chebyshev
      ! fine grids in tau and omega
      !
      ! Input:
      !
      ! lambda  - cutoff parameter
      ! p       - Chebyshev degree in each subinterval
      ! npt     - # subintervals on [0,1/2] in tau space (# subintervals
      !             on [0,1] is 2*npt)
      ! npo     - # subintervals on [0,lambda] in omega space (#
      !             subintervals on [-lambda,lambda] is 2*npo)
      !
      ! Output:
      !
      ! t       - tau fine grid points on (0,1/2) (half of full grid)
      ! om      - omega fine grid points
      ! kmat    - K(tau,omega) on fine grid 
      ! err     - Error of composite Chebyshev interpolant of
      !             K(tau,omega). err(1) is ~= max relative L^inf error
      !             in tau over all omega in fine grid. err(2) is ~= max
      !             L^inf error in omega over all tau in fine grid.

      implicit none
      integer p,npt,npo
      real *8 lambda,t(npt*p),om(2*npo*p)
      real *8 kmat(2*npt*p,2*npo*p),err(2)
      real *8, external :: kfunf

      integer nt,no,i,j,k
      real *8 one,a,b,start,finish,xx,ktrue,ktest,errtmp
      real *8, allocatable :: xc(:),wc(:),pbpt(:),pbpo(:)
      real *8, allocatable :: ktcoef(:,:),komcoef(:,:)
      real *8, allocatable :: xc2(:),wc2(:)

      one = 1.0d0

      nt = 2*npt*p
      no = 2*npo*p

      ! --- Chebyshev nodes and interpolation weights ---

      allocate(xc(p),wc(p))
      
      call barychebinit(p,xc,wc)

      ! -- Tau space discretization ---

      ! Panel break points

      allocate(pbpt(2*npt+1))

      pbpt(1) = 0*one
      do i=1,npt
        pbpt(i+1) = one/2**(npt-i+1)
      enddo

      !pbpt(npt+2:2*npt+1) = 1-pbpt(npt:1:-1)

      ! Grid points

      do i=1,npt
        a = pbpt(i)
        b = pbpt(i+1)
        t((i-1)*p+1:i*p) = a + (b-a)*(xc+one)/2
      enddo

      ! --- Omega space discretization ---

      ! Panel break points

      allocate(pbpo(2*npo+1))

      pbpo(npo+1) = 0*one
      do i=1,npo
        pbpo(npo+i+1) = lambda/2**(npo-i)
      enddo

      pbpo(1:npo) = -pbpo(2*npo+1:npo+2:-1)

      ! Grid points

      do i=1,2*npo
        a = pbpo(i)
        b = pbpo(i+1)
        om((i-1)*p+1:i*p) = a + (b-a)*(xc+one)/2
      enddo

      ! --- Sample K(tau,omega) on grid ---

      do j=1,no
        do i=1,nt/2

          kmat(i,j) = kfunf(t(i),om(j))

        enddo
      enddo

      ! Copy second half of matrix from first half to improve accuracy
      ! for extremely large npt: computing exp((1-t)*omega) loses digits
      ! if t is very close to 1 and (1-t)*omega ~ 1, but exp(-omega*t)
      ! is fine for small t and t*omega ~ 1.

      kmat(nt/2+1:nt,1:no) = kmat(nt/2:1:-1,no:1:-1)


      ! --- Check accuracy of Cheb interpolant on each panel in tau
      ! for fixed omega, and each panel in omega for fixed tau, by
      ! comparing with K(tau,omega) on Cheb grid of 2*p nodes ---

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

            call barycheb(p,xc2(k),kmat((i-1)*p+1:i*p,j),wc,xc,ktest)

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

            call barycheb(p,xc2(k),kmat(j,(i-1)*p+1:i*p),wc,xc,ktest)

            errtmp = max(errtmp,abs(ktrue-ktest))

          enddo
        enddo

        err(2) = max(err(2),errtmp/maxval(kmat(j,:)))

      enddo


      end subroutine kfine_cc


      
      subroutine dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx)

      ! Select real frequency nodes defining DLR basis
      !
      ! Input:
      !
      ! lambda  - cutoff parameter
      ! eps     - DLR error tolerance
      ! nt      - # fine grid points in tau
      ! no      - # fine grid points in omega
      ! om      - omega fine grid points
      ! kmat    - K(tau,omega) on fine grid
      ! rank    - max possible rank of DLR, defining input size of some
      !             arrays
      !
      ! Output :
      !
      ! rank    - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      ! oidx    - column indices of kmat corresponding to selected real
      !             frequency nodes

      implicit none
      integer nt,no,rank,oidx(rank)
      real *8 lambda,eps,om(no),kmat(nt,no),dlrrf(rank)

      integer, allocatable :: list(:)
      real *8, allocatable :: tmp(:,:),work(:)

      ! --- Select real frequency nodes by pivoted QR on columns of 
      ! kmat ---

      allocate(tmp(nt,no),list(no),work(max(nt,no)))

      tmp = kmat

      ! Pivoted QR 
      
      call iddp_qrpiv(eps,nt,no,tmp,rank,list,work)

      ! Rearrange indices to get selected frequency point indices

      call ind_rearrange(no,rank,list)

      ! Extract selected frequencies

      oidx(1:rank) = list(1:rank)
      dlrrf(1:rank) = om(oidx(1:rank))
          
      end subroutine dlr_rf


      subroutine dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit,tidx)

      ! Select imaginary time DLR nodes
      !
      ! Input:
      !
      ! lambda  - cutoff parameter
      ! nt      - # fine grid points in tau
      ! no      - # fine grid points in omega
      ! t       - tau fine grid points
      ! kmat    - K(tau,omega) on fine grid
      ! rank    - rank of DLR (# basis functions)
      ! oidx    - column indices of kmat corresponding to selected real
      !             frequency nodes
      !
      ! Output :
      !
      ! dlrit   - selected imaginary time nodes (tau points)
      ! tidx    - row indices of kmat corresponding to selected
      !             imaginary time nodes

      implicit none
      integer nt,no,rank,oidx(rank),tidx(rank)
      real *8 lambda,t(nt),kmat(nt,no),dlrit(rank)

      integer j,k
      integer, allocatable :: list(:)
      real *8, allocatable :: tmp(:,:),work(:)

      ! --- Select imaginary time nodes by pivoted QR on rows of 
      ! kmat ---

      ! Matrix of selected columns

      allocate(tmp(rank,nt),list(nt),work(nt))

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


      subroutine dlr_cf2it(nt,no,kmat,rank,oidx,tidx,cf2it)

      ! Build transform matrix from DLR coefficients to samples on
      ! imaginary time grid. To obtain the samples of a DLR expansion on
      ! the imaginary time grid, apply the matrix cf2it to the vector of
      ! DLR coefficients.
      !
      ! Input:
      !
      ! nt    - # fine grid points in tau
      ! no    - # fine grid points in omega
      ! kmat  - K(tau,omega) on fine grid
      ! rank  - rank of DLR (# basis functions)
      ! oidx  - column indices of kmat corresponding to selected
      !           real frequency nodes
      ! tidx  - row indices of kmat corresponding to selected imaginary
      !           time nodes
      !
      ! Output :
      !
      ! cf2it - DLR coefficients -> imaginary time grid values transform
      !           matrix


      implicit none
      integer nt,no,rank,tidx(rank),oidx(rank)
      real *8 kmat(nt,no),cf2it(rank,rank)

      integer j,k

      ! Extract select rows and columns of fine grid K matrix

      do k=1,rank
        do j=1,rank
          cf2it(j,k) = kmat(tidx(j),oidx(k))
        enddo
      enddo

      end subroutine dlr_cf2it


      subroutine dlr_it2cf(nt,no,kmat,rank,oidx,tidx,it2cf,it2cfpiv)

      ! Build transform matrix from samples on imaginary time grid to
      ! DLR coefficients, stored in LU form. To obtain the coefficients
      ! of a DLR expansion from samples on the imaginary time grid, use
      ! the outputs of this subroutine in conjunction with the dlr_expnd
      ! subroutine.
      !
      ! Input:
      !
      ! nt        - # fine grid points in tau
      ! no        - # fine grid points in omega
      ! kmat      - K(tau,omega) on fine grid
      ! rank      - rank of DLR (# basis functions)
      ! oidx      - column indices of kmat corresponding to selected
      !               real frequency nodes
      ! tidx      - row indices of kmat corresponding to selected
      !               imaginary time nodes
      !
      ! Output :
      !
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage format


      implicit none
      integer nt,no,rank,tidx(rank),oidx(rank),it2cfpiv(rank)
      real *8 kmat(nt,no),it2cf(rank,rank)

      integer j,k,info

      ! Extract select rows and columns of fine grid K matrix

      do k=1,rank
        do j=1,rank
          it2cf(j,k) = kmat(tidx(j),oidx(k))
        enddo
      enddo

      ! LU factorize

      call dgetrf(rank,rank,it2cf,rank,it2cfpiv,info)

      end subroutine dlr_it2cf


      subroutine dlr_expnd(rank,it2cf,it2cfpiv,g,gc)
      
      ! Get coefficients of DLR from samples on imaginary time DLR grid
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! it2cf     - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for it2cf in lapack LU storage format
      ! g         - Samples of a function G at imaginary time grid
      !               points
      !
      ! Output :
      !
      ! gc        - DLR coefficients of G
      
      implicit none
      integer rank,it2cfpiv(rank)
      real *8 it2cf(rank,rank),g(rank),gc(rank)

      integer info

      ! Backsolve with imaginary time grid values -> DLR coefficients
      ! transform matrix stored in LU form

      gc = g

      call dgetrs('N',rank,1,it2cf,rank,it2cfpiv,gc,rank,info)

      end subroutine dlr_expnd



      subroutine dlr_mfexpnd(rank,mf2cf,mf2cfpiv,g,gc)

      ! Get coefficients of DLR from samples on Matsubara frequency DLR
      ! grid
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! mf2cf     - Matsubara frequency grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! mf2cfpiv  - pivot matrix for mf2cf in lapack LU storage format
      ! g         - Samples of a function G at Matsubara frequency grid
      !               points
      !
      ! Output :
      !
      ! gc        - DLR coefficients of G

      implicit none
      integer rank,mf2cfpiv(rank)
      real *8 gc(rank)
      complex *16 mf2cf(rank,rank),g(rank)

      integer info
      complex *16, allocatable :: tmp(:)

      ! Backsolve with DLR transform matrix in factored form

      allocate(tmp(rank))

      tmp = g

      call zgetrs('N',rank,1,mf2cf,rank,mf2cfpiv,tmp,rank,info)

      gc = real(tmp)

      end subroutine dlr_mfexpnd


      subroutine dlr_eval(rank,dlrrf,g,t,val)

      ! Evaluate DLR expansion at a point t
      !
      ! Input:
      !
      ! rank    - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      ! g       - DLR coefficients of a function G
      ! t       - evaluation points in t' format
      !
      ! Output:
      !
      ! val     - value of DLR of G at t
      !
      ! Note: to evaluate at a point 0.5<t<1, input the value t' = t-1.
      ! If t' has been computed to high relative precision, this
      ! subroutine will avoid loss of digits for t very close to 1 by
      ! evaluating the kernel K using its symmetries.

      implicit none
      integer rank
      real *8 dlrrf(rank),g(rank),t,val

      integer i
      real *8 kval
      real *8, external :: kfunf

      val = 0.0d0
      do i=1,rank

        ! For 0.5<t<1, corresponding to negative t', use symmetry of K
        ! to evaluate basis functions

        if (t.ge.0.0d0) then
          kval = kfunf(t,dlrrf(i))
        else
          kval = kfunf(-t,-dlrrf(i))
        endif

        val = val + g(i)*kval

      enddo

      end subroutine dlr_eval



      subroutine dlr_mf_eval(rank,dlrrf,g,n,val)

      ! Evaluate DLR expansion at a point t
      !
      ! Input:
      !
      ! rank    - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      ! g       - DLR coefficients of a function G
      ! n       - evaluation point in Matsubara frequency
      !
      ! Output:
      !
      ! val     - value of DLR of G at n

      implicit none
      integer rank,n
      real *8 dlrrf(rank),g(rank)
      complex *16 val

      integer i
      complex *16 kval
      complex *16, external :: kfunf_mf

      val = 0.0d0
      do i=1,rank

        kval = kfunf_mf(n,dlrrf(i))

        val = val + g(i)*kval

      enddo

      end subroutine dlr_mf_eval


      subroutine dlr_mf(nmax,rank,dlrrf,dlrmf)

      ! Select Matsubara frequency DLR nodes
      !      
      ! Input:
      !
      ! nmax    - Matsubara frequency cutoff
      ! rank    - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      !
      ! Output :
      !
      ! dlrmf   - selected Matsubara frequency nodes

      implicit none
      integer nmax,rank,dlrmf(rank)
      real *8 dlrrf(rank)

      integer i,k,info
      integer, allocatable :: ns(:),list(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: poles(:,:)
      complex *16, external :: kfunf_mf

      ! Get matrix of Fourier transforms of DLR basis functions

      allocate(poles(rank,2*nmax+1),ns(2*nmax+1))

      ns = (/(i, i=-nmax,nmax)/)

      do i=1,2*nmax+1
        do k=1,rank
          
          poles(k,i) = kfunf_mf(ns(i),dlrrf(k))
          
        enddo
      enddo

      ! --- Select Matsubara frequency nodes by pivoted QR on rows of
      ! Fourier transformed K matrix ---

      allocate(list(2*nmax+1),work(2*nmax+1))

      ! Pivoted QR

      call idzr_qrpiv(rank,2*nmax+1,poles,rank,list,work)

      ! Rearrange indices to get selected frequency point indices

      call ind_rearrange(2*nmax+1,rank,list)

      ! Extract selected frequencies

      dlrmf = ns(list(1:rank))

      end subroutine dlr_mf


      subroutine dlr_mf2cf(nmax,rank,dlrrf,dlrmf,mf2cf,mf2cfpiv)

      ! Build transform matrix from samples on Matsubara frequency grid
      ! to DLR coefficients in LU form
      !
      ! Input:
      !
      ! nmax    - Matsubara frequency cutoff
      ! rank      - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      ! dlrmf   - selected Matsubara frequency nodes
      !
      ! Output :
      !
      ! mf2cf   - Matsubara frequency grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! mf2cfpiv  - pivot matrix for mf2cf in lapack LU storage format

      implicit none
      integer nmax,rank,dlrmf(rank),mf2cfpiv(rank)
      real *8 dlrrf(rank)
      complex *16 mf2cf(rank,rank)

      integer j,k,info
      complex *16, external :: kfunf_mf

      ! Extract selected rows and columns of Fourier transformed K
      ! matrix

      do k=1,rank
        do j=1,rank
          mf2cf(j,k) = kfunf_mf(dlrmf(j),dlrrf(k))
        enddo
      enddo

      ! LU factorize

      call zgetrf(rank,rank,mf2cf,rank,mf2cfpiv,info)

      end subroutine dlr_mf2cf


      subroutine dlr_cf2mf(rank,dlrrf,dlrmf,cf2mf)

      ! Build transform matrix from DLR coefficients to samples on
      ! Matsubara frequency grid. To obtain the samples of a DLR
      ! expansion on the Matsubara frequency grid, apply the matrix
      ! cf2mf to the vector of DLR coefficients.
      !
      ! Input:
      !
      ! rank  - rank of DLR (# basis functions)
      ! dlrrf - selected real frequency nodes (omega points)
      ! dlrmf - selected Matsubara frequency nodes
      !
      ! Output :
      !
      ! cf2mf - DLR coefficients -> Matsubara frequency grid values
      !           transform matrix


      implicit none
      integer rank,dlrmf(rank)
      real *8 dlrrf(rank)
      complex *16 cf2mf(rank,rank)

      complex *16, external :: kfunf_mf

      integer i,j

      ! Evaluated Matsubara frequency kernel at selected real
      ! frequencies and Matsubara frequencies 

      do j=1,rank
        do i=1,rank
          cf2mf(i,j) = kfunf_mf(dlrmf(i),dlrrf(j))
        enddo
      enddo

      end subroutine dlr_cf2mf




      subroutine dlr_convtens(rank,dlrrf,dlrit,phi)

      ! Get tensor phi_{jkl} used to take a set of DLR coefficients to
      ! the matrix of convolution by the corresponding DLR expansion.
      ! The matrix A of convolution on [0,beta] by a function
      ! represented by a DLR expansion with coefficients rho_l is given
      ! by
      !
      ! A_jk = beta * sum_l phi_jkl rho_l.
      !
      ! This is for the fermionic case only.
      !
      ! Input:
      !
      ! rank  - rank of DLR (# basis functions)
      ! dlrrf - selected real frequency nodes (omega points)
      ! dlrit - selected imaginary time nodes (tau points)
      !
      ! Output:
      !
      ! phi   - convolution tensor


      implicit none
      integer rank
      real *8 dlrrf(rank),dlrit(rank)
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

      end subroutine dlr_convtens


      subroutine dlr_convmat(rank,phi,it2cf,it2cfpiv,g,gmat)

      ! Get matrix of convolution by a DLR expansion G in the DLR basis
      ! -- that is, the matrix that this subroutine produces takes the
      ! DLR coefficient representation of a function f to the DLR
      ! coefficient representation of the convolution
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
      ! g         - DLR coefficients of a function G
      !
      ! Output:
      !
      ! gmat      - matrix of convolution by G in the DLR basis

      implicit none
      integer rank,it2cfpiv(rank)
      real *8 phi(rank*rank,rank),it2cf(rank,rank),g(rank)
      real *8 gmat(rank,rank)

      integer i,j,info

      call dgemv('N',rank*rank,rank,1.0d0,phi,rank*rank,g,1,0.0d0,&
        gmat,1)

      call dgetrs('N',rank,rank,it2cf,rank,it2cfpiv,gmat,rank,info)

      end subroutine dlr_convmat





      subroutine eqpts_rel(n,t)

      ! Get equispaced points on [0,1] in relative format
      !
      ! Relative format means that points 0.5<t<1 are computed and stored as the negative
      ! distance from 1; that is, t* = t-1 in exact arithmetic. This is
      ! to used to maintain full relative precision for calculations
      ! with large lambda and small eps.
      !
      ! Input:
      !
      ! n - Number of points on [0,1]
      !
      ! Output :
      !
      ! t - n equispaced points on [0,1], including endpoints, in
      !       relative format

      implicit none
      integer n
      real *8 t(n)

      integer i

      do i=1,n-1

        if (i.le.n/2) then
          t(i) = (i-1)*1.0d0/(n-1)
        else
          t(i) = -(n-i)*1.0d0/(n-1)
        endif

      enddo

      t(n) = 1.0d0

      end subroutine eqpts_rel



      subroutine rel2abs(n,t,tabs)

      ! Convert points on [0,1] from relative format to absolute format
      !
      ! Relative format means that points 0.5<t<1 are computed and stored as the negative
      ! distance from 1; that is, t* = t-1 in exact arithmetic. This is
      ! to used to maintain full relative precision for calculations
      ! with large lambda and small eps. Absolute format means that all
      ! points are stored as normal.
      !
      ! Note: converting a point from relative to absolute format will,
      ! in general, constitute a loss of relative accuracy in the
      ! location of the point if the point is close to t = 1. For
      ! example, in three-digit arithmetic, the point t = 0.999111 could
      ! be stored as t* = -0.889e-3 in the relative format, but only as
      ! t = 0.999 in the absolute format.
      !
      ! Input:
      !
      ! n     - Number of points
      ! t     - Array of points on [0,1] stored in relative format
      !
      ! Output:
      !
      ! trel  - Array of points t in absolute format

      implicit none
      integer n
      real *8 t(n),tabs(n)

      integer i

      do i=1,n

        if (t(i).lt.0.0d0) then
          tabs(i) = t(i)+1.0d0
        else
          tabs(i) = t(i)
        endif

      enddo

      end subroutine rel2abs


      subroutine abs2rel(n,tabs,t)

      ! Convert a point on [0,1] from absolute format to relative format
      !
      ! Relative format means that points 0.5<t<1 are computed and stored as the negative
      ! distance from 1; that is, t* = t-1 in exact arithmetic. This is
      ! to used to maintain full relative precision for calculations
      ! with large lambda and small eps. Absolute format means that all
      ! points are stored as normal.
      !
      ! If the user wishes to specify points -- for example points at
      ! which to sample or evaluate a DLR -- in absolute format, those
      ! points must first be converted to relative format using this
      ! subroutine before using them as inputs into any other
      ! subroutines. Of course, in order to maintain full relative
      ! precision in all calculations, the user must specify points in
      ! relative format from the beginning, but in most cases at most a
      ! mild loss of accuracy will result from using the absolute
      ! format.
      !
      ! Input:
      !
      ! n     - Number of points
      ! tabs  - Array of points on [0,1] stored in absolute format
      !
      ! Output:
      !
      ! t     - Array of points t in relative format

      implicit none
      integer n
      real *8 t(n),tabs(n)

      integer i

      do i=1,n

        if (t(i).gt.0.5d0.and.t(i).lt.1.0d0) then
          t(i) = tabs(i)-1.0d0
        else
          t(i) = tabs(i)
        endif

      enddo

      end subroutine abs2rel

