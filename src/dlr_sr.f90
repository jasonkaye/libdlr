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
      
      npt = ceiling(log(lambda)/log(2.0d0))-2
      npo = ceiling(log(lambda)/log(2.0d0))

      nt = 2*p*npt
      no = 2*p*npo

      end subroutine gridparams


      subroutine kfine_cc(fb,lambda,p,npt,npo,t,om,kmat,err)

      ! Discretization of kernel K(tau,omega) on composite Chebyshev
      ! fine grids in tau and omega
      !
      ! Note: this is a wrapper for the main subroutine, kfine_cc1
      !
      ! Input:
      !
      ! fb      - Fermionic (f) or bosonic (b) kernel
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
      real *8, external :: kfunf,kfunb
      character :: fb

      if (fb.eq.'f') then
        call kfine_cc1(kfunf,lambda,p,npt,npo,t,om,kmat,err)
      elseif (fb.eq.'b') then
        call kfine_cc1(kfunb,lambda,p,npt,npo,t,om,kmat,err)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine kfine_cc


      subroutine kfine_cc1(kfun,lambda,p,npt,npo,t,om,kmat,err)

      ! Main subroutine for kfine_cc
      !
      ! See subroutine kfine_cc for description of arguments, except for
      ! kfun, which indicates the kernel evaluator to be used, and
      ! depends on whether fb = 'f' or fb = 'b'.

      implicit none
      integer p,npt,npo
      real *8 lambda,t(npt*p),om(2*npo*p)
      real *8 kmat(2*npt*p,2*npo*p),err(2)
      real *8, external :: kfun

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

          kmat(i,j) = kfun(t(i),om(j))

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
            
            ktrue = kfun(xx,om(j))

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
            
            ktrue = kfun(t(j),xx)

            call barycheb(p,xc2(k),kmat(j,(i-1)*p+1:i*p),wc,xc,ktest)

            errtmp = max(errtmp,abs(ktrue-ktest))

          enddo
        enddo

        err(2) = max(err(2),errtmp/maxval(kmat(j,:)))

      enddo

      end subroutine kfine_cc1


      
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


      subroutine dlr_it2cf(nt,no,kmat,rank,oidx,tidx,dlrit2cf,it2cfpiv)

      ! Build transform matrix from samples on imaginary time grid to
      ! DLR coefficients in LU form
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
      ! dlrit2cf  - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for dlrit2cf in lapack LU storage format


      implicit none
      integer nt,no,rank,tidx(rank),oidx(rank),it2cfpiv(rank)
      real *8 kmat(nt,no),dlrit2cf(rank,rank)

      integer j,k,info

      ! Extract select rows and columns of fine grid K matrix

      do k=1,rank
        do j=1,rank
          dlrit2cf(j,k) = kmat(tidx(j),oidx(k))
        enddo
      enddo

      ! LU factorize

      call dgetrf(rank,rank,dlrit2cf,rank,it2cfpiv,info)

      end subroutine dlr_it2cf


      subroutine dlr_expnd(rank,dlrit2cf,it2cfpiv,g)
      
      ! Get coefficients of DLR from samples on imaginary time DLR grid
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! dlrit2cf  - imaginary time grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! it2cfpiv  - pivot matrix for dlrit2cf in lapack LU storage format
      ! g         - Samples of a function G at imaginary time grid
      !               points
      !
      ! Output :
      !
      ! g         - DLR coefficients of G
      
      implicit none
      integer rank,it2cfpiv(rank)
      real *8 dlrit2cf(rank,rank),g(rank)

      integer info

      ! Backsolve with imaginary time grid values -> DLR coefficients
      ! transform matrix stored in LU form

      call dgetrs('N',rank,1,dlrit2cf,rank,it2cfpiv,g,rank,info)

      end subroutine dlr_expnd



      subroutine dlr_mfexpnd(rank,dlrmf2cf,mf2cfpiv,g)

      ! Get coefficients of DLR from samples on Matsubara frequency DLR
      ! grid
      !
      ! Input:
      !
      ! rank      - rank of DLR (# basis functions)
      ! dlrmf2cf  - Matsubara frequency grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! mf2cfpiv  - pivot matrix for dlrmf2cf in lapack LU storage format
      ! g         - Samples of a function G at Matsubara frequency grid
      !               points
      !
      ! Output :
      !
      ! g         - DLR coefficients of G

      implicit none
      integer rank,mf2cfpiv(rank)
      complex *16 dlrmf2cf(rank,rank),g(rank)

      integer info

      ! Backsolve with DLR transform matrix in factored form

      call zgetrs('N',rank,1,dlrmf2cf,rank,mf2cfpiv,g,rank,info)

      end subroutine dlr_mfexpnd


      subroutine dlr_eval(fb,rank,dlrrf,g,t,val)

      ! Evaluate DLR expansion at a point t
      !
      ! Input:
      !
      ! fb      - Fermionic (f) or bosonic (b) kernel
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
      !
      ! This is a wrapper for the main subroutine, dlr_eval1

      implicit none
      integer rank
      real *8 dlrrf(rank),g(rank),t,val
      character :: fb

      real *8, external :: kfunf,kfunb

      if (fb.eq.'f') then
        call dlr_eval1(rank,kfunf,dlrrf,g,t,val)
      elseif (fb.eq.'b') then
        call dlr_eval1(rank,kfunb,dlrrf,g,t,val)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_eval


      subroutine dlr_eval1(rank,kfun,dlrrf,g,t,val)

      ! Main subroutine for dlr_eval
      !
      ! See subroutine dlr_eval for description of arguments, except for
      ! kfun, which indicates the kernel evaluator to be used, and
      ! depends on whether fb = 'f' or fb = 'b'.

      implicit none
      integer rank
      real *8 dlrrf(rank),g(rank),t,val
      real *8, external :: kfun

      integer i
      real *8 kval

      val = 0.0d0
      do i=1,rank

        ! For 0.5<t<1, corresponding to negative t', use symmetry of K
        ! to evaluate basis functions

        if (t.ge.0.0d0) then
          kval = kfun(t,dlrrf(i))
        else
          kval = kfun(-t,-dlrrf(i))
        endif

        val = val + g(i)*kval

      enddo

      end subroutine dlr_eval1



      subroutine dlr_mf_eval(fb,rank,dlrrf,g,n,val)

      ! Evaluate DLR expansion at a point t
      !
      ! Input:
      !
      ! fb      - Fermionic (f) or bosonic (b) kernel
      ! rank    - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      ! g       - DLR coefficients of a function G
      ! n       - evaluation point in Matsubara frequency
      !
      ! Output:
      !
      ! val     - value of DLR of G at n
      !
      ! This is a wrapper for the main subroutine, dlr_mf_eval1

      implicit none
      integer rank,n
      real *8 dlrrf(rank),g(rank)
      complex *16 val
      character :: fb

      complex *16, external :: kfunf_mf

      if (fb.eq.'f') then
        call dlr_mf_eval1(rank,kfunf_mf,dlrrf,g,n,val)
      elseif (fb.eq.'b') then
        stop 'choose fb = b not supported yet'
        !call dlr_mf_eval1(rank,kfunb,dlrrf,g,n,val)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_mf_eval


      subroutine dlr_mf_eval1(rank,kfun,dlrrf,g,n,val)

      ! Main subroutine for dlr_mf_eval
      !
      ! See subroutine dlr_mf_eval for description of arguments, except for
      ! kfun, which indicates the kernel evaluator to be used, and
      ! depends on whether fb = 'f' or fb = 'b'.

      implicit none
      integer rank,n
      real *8 dlrrf(rank),g(rank)
      complex *16 val
      complex *16, external :: kfun

      integer i
      complex *16 kval

      val = 0.0d0
      do i=1,rank

        kval = kfun(n,dlrrf(i))

        val = val + g(i)*kval

      enddo

      end subroutine dlr_mf_eval1



      subroutine dlr_mf(fb,nmax,rank,dlrrf,dlrmf)

      ! Select Matsubara frequency DLR nodes
      !      
      ! Input:
      !
      ! fb      - Fermionic (f) or bosonic (b) kernel
      ! nmax    - Matsubara frequency cutoff
      ! rank    - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      !
      ! Output :
      !
      ! dlrmf   - selected Matsubara frequency nodes
      !
      !
      ! This is a wrapper for the main subroutine, dlr_mf1

      implicit none
      integer nmax,rank,dlrmf(rank)
      real *8 dlrrf(rank)
      character :: fb

      complex *16, external :: kfunf_mf

      if (fb.eq.'f') then
        call dlr_mf1(kfunf_mf,nmax,rank,dlrrf,dlrmf)
      !elseif (fb.eq.'b') then
        !call dlr_mf1(kfunb_mf,nmax,rank,dlrrf,dlrmf)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_mf


      subroutine dlr_mf1(kfun,nmax,rank,dlrrf,dlrmf)

      ! Main subroutine for dlr_mf
      !
      ! See subroutine dlr_mf for description of arguments, except for
      ! kfun, which indicates the kernel evaluator to be used, and
      ! depends on whether fb = 'f' or fb = 'b'.

      implicit none
      integer rank,nmax,dlrmf(rank)
      real *8 dlrrf(rank)
      complex *16, external :: kfun

      integer i,k,info
      integer, allocatable :: ns(:),list(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: poles(:,:)

      ! Get matrix of Fourier transforms of DLR basis functions

      allocate(poles(rank,2*nmax+1),ns(2*nmax+1))

      ns = (/(i, i=-nmax,nmax)/)

      do i=1,2*nmax+1
        do k=1,rank
          
          poles(k,i) = kfun(ns(i),dlrrf(k))
          
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

      end subroutine dlr_mf1






      subroutine dlr_mf2cf(fb,nmax,rank,dlrrf,dlrmf,dlrmf2cf,mf2cfpiv)

      ! Build transform matrix from samples on Matsubara frequency grid
      ! to DLR coefficients in LU form
      !
      ! Input:
      !
      ! fb      - Fermionic (f) or bosonic (b) kernel
      ! nmax    - Matsubara frequency cutoff
      ! rank      - rank of DLR (# basis functions)
      ! dlrrf   - selected real frequency nodes (omega points)
      ! dlrmf   - selected Matsubara frequency nodes
      !
      ! Output :
      !
      ! dlrmf2cf  - Matsubara frequency grid values -> DLR coefficients
      !               transform matrix in lapack LU storage format
      ! mf2cfpiv  - pivot matrix for dlrmf2cf in lapack LU storage format
      !
      !
      ! This is a wrapper for the main subroutine, dlr_mf2cf1

      implicit none
      integer nmax,rank,dlrmf(rank),mf2cfpiv(rank)
      real *8 dlrrf(rank)
      complex *16 dlrmf2cf(rank,rank)
      character :: fb

      complex *16, external :: kfunf_mf

      if (fb.eq.'f') then
        call dlr_mf2cf1(kfunf_mf,nmax,rank,dlrrf,dlrmf,dlrmf2cf,&
          mf2cfpiv)
      !elseif (fb.eq.'b') then
        !call dlr_mf2cf1(kfunb_mf,nmax,rank,dlrrf,dlrmf,dlrmf2cf,&
        ! mf2cfpiv)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_mf2cf


      subroutine dlr_mf2cf1(kfun,nmax,rank,dlrrf,dlrmf,dlrmf2cf,&
          mf2cfpiv)

      ! Main subroutine for dlr_mf2cf
      !
      ! See subroutine dlr_mf2cf for description of arguments, except for
      ! kfun, which indicates the kernel evaluator to be used, and
      ! depends on whether fb = 'f' or fb = 'b'.

      implicit none
      integer nmax,rank,dlrmf(rank),mf2cfpiv(rank)
      real *8 dlrrf(rank)
      complex *16 dlrmf2cf(rank,rank)
      complex *16, external :: kfun

      integer j,k,info

      ! Extract selected rows and columns of Fourier transformed K
      ! matrix

      do k=1,rank
        do j=1,rank
          dlrmf2cf(j,k) = kfun(dlrmf(j),dlrrf(k))
        enddo
      enddo

      ! LU factorize

      call zgetrf(rank,rank,dlrmf2cf,rank,mf2cfpiv,info)

      end subroutine dlr_mf2cf1
