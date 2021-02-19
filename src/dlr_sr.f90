      ! ----- Subroutines for intermediate representation and discrete
      ! Lehmann representation -----


      ! ----- Kernel evaluator -----

      real *8 function kfunf(t,om)

      implicit none
      real *8 t,om

      if (om.ge.0.0d0) then
        kfunf = exp(-t*om)/(1.0d0+exp(-om))
      else
        kfunf = exp((1.0d0-t)*om)/(1.0d0+exp(om))
      endif

      end function kfunf


      real *8 function kfunb(t,om)

      implicit none
      real *8 t,om

      if (om.ge.0.0d0) then
        call evalexpfun(om,kfunb)
        kfunb = exp(-t*om)/kfunb
      else
        call evalexpfun(-om,kfunb)
        kfunb = exp((1.0d0-t)*om)/kfunb
      endif

      end function kfunb


      complex *16 function kfunf_mf(n,om)

      implicit none
      integer n
      real *8 om

      real *8 pi
      complex *16 eye

      pi = 4*atan(1.0d0)
      eye = (0.0d0,1.0d0)

      !kfunf_mf = tanh(om/2)/(2*pi*eye*n+om)
      kfunf_mf = 1.0d0/((2*n+1)*pi*eye+om)

      end function kfunf_mf


      real *8 function kfunf2(t,om)

      ! Stable evaluation of fermionic K(tau,omega) with tau<0 treated as 1+tau

      implicit none
      real *8 t,om

      real *8, external :: kfunf

        if (t.ge.0.0d0) then
          kfunf2 = kfunf(t,om)
        else
          kfunf2 = kfunf(-t,-om)
        endif

      end function kfunf2


      subroutine gridparams(lambda,p,npt,npo,nt,no)

      ! ----- Set fine grid parameters based on Lambda -----

      implicit none
      integer p,npt,npo,nt,no
      real *8 lambda

      p = 24 ! Chebyshev degree of panels
      
      npt = ceiling(log(lambda)/log(2.0d0))-2 ! Half the # panels in fine tau grid
      npo = ceiling(log(lambda)/log(2.0d0)) ! Half the # panels in fine omega grid

      nt = 2*p*npt ! Total # fine grid points in tau
      no = 2*p*npo ! Total # fine grid points in omega

      end subroutine gridparams


      subroutine kfine_cc(lambda,eps,fb,npt,npo,p,t,om,kmat,err)

      ! Get fine discretization of kernel K(tau,omega) on
      ! dyadically-refined composite Chebyshev grids on tau and omega
      !
      ! This is a wrapper for the main subroutine, kfine_cc1

      implicit none
      integer npt,npo,p
      real *8 lambda,eps,t(npt*p),om(2*npo*p)
      real *8 kmat(2*npt*p,2*npo*p),err(2)
      real *8, external :: kfunf,kfunb
      character :: fb

      if (fb.eq.'f') then
        call kfine_cc1(lambda,eps,kfunf,npt,npo,p,t,om,kmat,err)
      elseif (fb.eq.'b') then
        call kfine_cc1(lambda,eps,kfunb,npt,npo,p,t,om,kmat,err)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine kfine_cc


      subroutine kfine_cc1(lambda,eps,kfun,npt,npo,p,t,om,kmat,err)

      ! Get fine discretization of kernel K(tau,omega) on
      ! dyadically-refined composite Chebyshev grids on tau and omega
      !
      ! Main subroutine

      implicit none
      integer npt,npo,p
      real *8 lambda,eps,t(npt*p),om(2*npo*p)
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

      ! --- Get Chebyshev nodes and transform matrix ---

      allocate(xc(p),wc(p))
      
      call barychebinit(p,xc,wc)

      ! ----- Get tau space discretization -----

      ! --- Panel break points ---

      allocate(pbpt(2*npt+1))

      pbpt(1) = 0*one
      do i=1,npt
        pbpt(i+1) = one/2**(npt-i+1)
      enddo

      !pbpt(npt+2:2*npt+1) = 1-pbpt(npt:1:-1)

      ! --- Grid points ---

      do i=1,npt
        a = pbpt(i)
        b = pbpt(i+1)
        t((i-1)*p+1:i*p) = a + (b-a)*(xc+one)/2
      enddo

      ! ----- Get omega space discretization -----

      ! --- Panel break points ---

      allocate(pbpo(2*npo+1))

      pbpo(npo+1) = 0*one
      do i=1,npo
        pbpo(npo+i+1) = lambda/2**(npo-i)
      enddo

      pbpo(1:npo) = -pbpo(2*npo+1:npo+2:-1)

      ! --- Grid points ---

      do i=1,2*npo
        a = pbpo(i)
        b = pbpo(i+1)
        om((i-1)*p+1:i*p) = a + (b-a)*(xc+one)/2
      enddo

      ! ----- Sample K(tau,omega) on grid -----

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


      ! ----- Check accuracy of Cheb interpolant on each panel in tau
      ! for fixed omega, and each panel in omega for fixed tau, by
      ! comparing with K(tau,omega) on Cheb grid of 2*p nodes -----

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




      subroutine kcgl(lambda,eps,kfun,nt,no,p,t,om,twgt,pbpt,kmat)

      ! ----- Get discretization of kernel K(tau,omega) on
      ! dyadically-refined CGL grids on tau and omega -----

      implicit none
      integer nt,no,p
      real *8 lambda,eps,t(2*nt*p),om(2*no*p),twgt(2*nt*p)
      real *8 pbpt(2*nt+1),kmat(2*nt*p,2*no*p)
      real *8, external :: kfun

      integer nnt,nno,i,j
      real *8 one,a,b,start,finish
      real *8, allocatable :: xgl(:),legf(:,:),legb(:,:),wgl(:)
      real *8, allocatable :: pbpo(:)
      real *8, allocatable :: err1(:,:),err2(:,:)
      real *8, allocatable :: ktcoef(:,:),komcoef(:,:)

      one = 1.0d0

      ! --- Derived parameters ---

      nnt = 2*nt*p
      nno = 2*no*p

      ! --- Get Gauss-Legendre nodes and transforms ---

      allocate(xgl(p),legf(p,p),legb(p,p),wgl(p))
      
      call legeexps(2,p,xgl,legf,legb,wgl)

      ! ----- Get tau space discretization -----

      ! --- Panel break points ---

      pbpt(1) = 0*one
      do i=1,nt
        pbpt(i+1) = one/2**(nt-i+1)
      enddo

      pbpt(nt+2:2*nt+1) = 1-pbpt(nt:1:-1)

      ! --- Grid points ---

      do i=1,2*nt
        a = pbpt(i)
        b = pbpt(i+1)
        t((i-1)*p+1:i*p) = a + (b-a)*(xgl+one)/2
        twgt((i-1)*p+1:i*p) = wgl*(b-a)/2
      enddo

      ! ----- Get omega space discretization -----

      ! --- Panel break points ---

      allocate(pbpo(2*no+1))

      pbpo(no+1) = 0*one
      do i=1,no
        !pbpo(no+i+1) = lambda/2**(no-i+1)
        pbpo(no+i+1) = lambda/2**(no-i)
      enddo

      pbpo(1:no) = -pbpo(2*no+1:no+2:-1)


      ! --- Grid points ---

      do i=1,2*no
        a = pbpo(i)
        b = pbpo(i+1)
        om((i-1)*p+1:i*p) = a + (b-a)*(xgl+one)/2
      enddo



      ! ----- Sample K(tau,omega) on grid -----

      do j=1,nno
        do i=1,nnt/2

          kmat(i,j) = kfun(t(i),om(j))

        enddo
      enddo

      ! Copy second half of matrix from first half to improve accuracy
      ! for extremely large nt: computing exp((1-t)*omega) loses digits
      ! if t is very close to 1 and (1-t)*omega ~ 1, but exp(-omega*t)
      ! is fine for small t and t*omega ~ 1.

      kmat(nnt/2+1:nnt,1:nno) = kmat(nnt/2:1:-1,nno:1:-1)

      ! --- Estimate error of discretization ---

      allocate(ktcoef(nnt,nno),komcoef(nno,nnt))

      call cgl_pt2coef(p,2*nt,nno,legf,kmat,ktcoef)

      call cgl_pt2coef(p,2*no,nnt,legf,transpose(kmat),komcoef)

      allocate(err1(2*nt,nno),err2(2*no,nnt))
      
      do j=1,nno
        do i=1,2*nt
          err1(i,j) = sum(abs(ktcoef((i-1)*p+p-4:i*p,j)))
        enddo
      enddo


      do j=1,nnt
        do i=1,2*no
          err2(i,j) = sum(abs(komcoef((i-1)*p+p-4:i*p,j)))
        enddo
      enddo

      write(6,*) '---------------------------------------------------'
      write(6,*) 'Using ',nnt,' fine discretization nodes in tau'
      write(6,*) 'Using ',nno,' fine discretization nodes in omega'
      write(6,*) ''
      write(6,*) 'Error of K(tau,omega) discretization in tau ~= ',&
        maxval(err1)
      write(6,*) 'Error of K(tau,omega) discretization in omega ~= ',&
        maxval(err2)
      write(6,*) '---------------------------------------------------'
      write(6,*) ''


      end subroutine kcgl

      
      
      subroutine dlr_rf(lambda,eps,nt,no,om,kmat,rank,opts,oidx)

      ! Get real frequency DLR points

      ! On input, rank is maximum possible rank and defines the sizes of
      ! tpts and opts. On output, rank is the actual rank.

      implicit none
      integer nt,no,rank,oidx(rank)
      real *8 lambda,eps,om(no),kmat(nt,no),opts(rank)

      integer i,j,k,info
      integer, allocatable :: list(:)
      real *8, allocatable :: tmp(:,:),work(:)

      ! --- Get selected frequency points by pivoted QR on columns of
      ! fine grid K matrix ---

      allocate(tmp(nt,no),list(no),work(max(nt,no)))

      tmp = kmat

      ! Pivoted QR 
      
      call iddp_qrpiv(eps,nt,no,tmp,rank,list,work)

      ! Rearrange indices to get selected frequency point indices

      call ind_rearrange(no,rank,list)

      ! Extract selected frequencies

      oidx = list(1:rank)
      opts = om(oidx)
          
      end subroutine dlr_rf


      subroutine dlr_it(lambda,nt,no,rank,t,kmat,oidx,dlrit,tidx)

      ! Get DLR imaginary time points

      implicit none
      integer nt,no,rank,oidx(rank),tidx(rank)
      real *8 lambda,t(nt),kmat(nt,no),dlrit(rank)

      integer i,j,k,info
      integer, allocatable :: list(:)
      real *8, allocatable :: tmp(:,:),work(:)

      ! --- Get selected imaginary time points by pivoted QR on rows of
      ! selected columns ---

      ! Get matrix of selected columns, transposed

      allocate(tmp(rank,nt),list(nt),work(nt))

      do j=1,nt
        do k=1,rank
          tmp(k,j) = kmat(j,oidx(k))
        enddo
      enddo

      ! Pivoted QR

      call iddr_qrpiv(rank,nt,tmp,rank,list,work)

      ! Rearrange indices to get selected imaginary time point indices

      call ind_rearrange(nt,rank,list)

      ! Extract selected imaginary times. To maintain high precision for
      ! extremely large Lambda and small eps calculations, if a t was
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


      subroutine dlr_it2cf(nt,no,rank,kmat,tidx,oidx,dlrt2c,ipiv)

      ! Get transform matrix from imaginary time grid to DLR
      ! coefficients in LU form

      implicit none
      integer nt,no,rank,ipiv(rank),tidx(rank),oidx(rank)
      real *8 kmat(nt,no),dlrt2c(rank,rank)

      integer j,k,info

      ! Extract select rows and columns of fine grid K matrix

      do k=1,rank
        do j=1,rank
          dlrt2c(j,k) = kmat(tidx(j),oidx(k))
        enddo
      enddo

      ! LU factorize

      call dgetrf(rank,rank,dlrt2c,rank,ipiv,info)

      end subroutine dlr_it2cf






      subroutine dlr_expnd(rank,dlrmat,ipiv,g)

      ! Get coefficients of DLR from samples on imaginary time DLR grid

      ! On input, g contains samples of function at coarse grid points.
      ! On output, it contains the coefficients of the DLR expansion.

      implicit none
      integer rank,ipiv(rank)
      real *8 dlrmat(rank,rank),g(rank)

      integer info

      ! Backsolve with DLR transform matrix in factored form

      call dgetrs('N',rank,1,dlrmat,rank,ipiv,g,rank,info)

      end subroutine dlr_expnd



      subroutine dlr_mfexpnd(rank,dlrmf2c,ipiv,g)

      ! Get coefficients of DLR from samples on Matsubara frequency DLR grid

      ! On input, g contains samples of function at coarse grid points.
      ! On output, it contains the coefficients of the DLR expansion.

      implicit none
      integer rank,ipiv(rank)
      complex *16 dlrmf2c(rank,rank),g(rank)

      integer info

      ! Backsolve with DLR transform matrix in factored form

      call zgetrs('N',rank,1,dlrmf2c,rank,ipiv,g,rank,info)

      end subroutine dlr_mfexpnd


      subroutine dlr_eval(rank,fb,opts,g,t,val)

      ! Evaluate DLR expansion at a point t
      !
      ! Note: to evaluate at a point 0.5<t<=1, input the value t* = t-1.
      ! If t* has been computed to high relative precision, this
      ! subroutine will avoid loss of digits for t very close to 1 by
      ! evaluating the kernel K using its symmetries.
      !
      ! This is a wrapper for the main subroutine, dlr_eval1

      implicit none
      integer rank
      real *8 opts(rank),g(rank),t,val
      character :: fb

      real *8, external :: kfunf,kfunb

      if (fb.eq.'f') then
        call dlr_eval1(rank,kfunf,opts,g,t,val)
      elseif (fb.eq.'b') then
        call dlr_eval1(rank,kfunb,opts,g,t,val)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_eval


      subroutine dlr_eval1(rank,kfun,opts,g,t,val)

      ! Evaluate DLR expansion at a point
      !
      ! Main subroutine

      implicit none
      integer rank
      real *8 opts(rank),g(rank),t,val
      real *8, external :: kfun

      integer i
      real *8 kval

      val = 0.0d0
      do i=1,rank

        if (t.ge.0.0d0) then
          kval = kfun(t,opts(i))
        else
          kval = kfun(-t,-opts(i))
        endif

        val = val + g(i)*kval

      enddo

      end subroutine dlr_eval1



      !subroutine zdlreval(rank,kfun,ompts,g,t,val)

      !! ----- Evaluate DLR expansion at a point -----

      !implicit none
      !integer rank
      !real *8 ompts(rank),t,kval
      !complex *16 g(rank),val
      !real *8, external :: kfun

      !integer i

      !!val = 0
      !!do i=1,rank
      !!  val = val + g(i)*kfun(t,ompts(i))
      !!enddo

      !val = 0
      !do i=1,rank

      !  if (t.ge.0.0d0) then
      !    kval = kfun(t,ompts(i))
      !  else
      !    kval = kfun(-t,-ompts(i))
      !  endif

      !  val = val + g(i)*kval

      !enddo

      !end subroutine zdlreval

      subroutine dlr_mf(rank,dlro,nmax,fb,mfpts)

      ! Get Matsubara frequency grid for DLR
      !
      ! This is a wrapper for the main subroutine, dlr_mf1

      implicit none
      integer rank,nmax,mfpts(rank)
      real *8 dlro(rank)
      character :: fb

      complex *16, external :: kfunf_mf

      if (fb.eq.'f') then
        call dlr_mf1(rank,dlro,nmax,kfunf_mf,mfpts)
      !elseif (fb.eq.'b') then
        !call dlr_mf1(rank,dlro,nmax,kfun,mfpts)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_mf


      subroutine dlr_mf1(rank,dlro,nmax,kfun,mfpts)

      ! ----- Get Matsubara frequency grid for DLR basis -----

      implicit none
      integer rank,nmax,mfpts(rank)
      real *8 dlro(rank)
      complex *16, external :: kfun

      integer i,k,n,info
      integer, allocatable :: ns(:),list(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: poles(:,:),tmp(:,:)

      ! --- Get Fourier transforms of DLR basis functions ---

      allocate(poles(rank,2*nmax+1),ns(2*nmax+1))

      ns = (/(i, i=-nmax,nmax)/)

      do i=1,2*nmax+1
        do k=1,rank
          
          poles(k,i) = kfun(ns(i),dlro(k))
          
        enddo
      enddo

      ! --- Pivoted QR to select Matsubara frequency points ---

      allocate(tmp(rank,2*nmax+1),list(2*nmax+1),work(2*nmax+1))

      tmp = poles

      call idzr_qrpiv(rank,2*nmax+1,tmp,rank,list,work)

      call ind_rearrange(2*nmax+1,rank,list)

      mfpts = ns(list(1:rank))

      end subroutine dlr_mf1


      subroutine dlr_mf2cf(nmax,rank,dlro,dlrmf,fb,mf2cf,ipiv)

      ! Get Matsubara frequency values -> DLR coefficients transform matrix in
      ! LU form
      !
      ! This is a wrapper for the main subroutine, dlr_mf2cf1

      implicit none
      integer nmax,rank,dlrmf(rank),ipiv(rank)
      real *8 dlro(rank)
      complex *16 mf2cf(rank,rank)
      character :: fb

      complex *16, external :: kfunf_mf

      if (fb.eq.'f') then
        call dlr_mf2cf1(nmax,rank,dlro,dlrmf,kfunf_mf,mf2cf,ipiv)
      !elseif (fb.eq.'b') then
        !call dlr_mf2cf1(nmax,rank,dlro,dlrmf,kfun,mf2cf,ipiv)
      else
        stop 'choose fb = f or b'
      endif

      end subroutine dlr_mf2cf


      subroutine dlr_mf2cf1(nmax,rank,dlro,dlrmf,kfun,mf2cf,ipiv)

      implicit none
      integer nmax,rank,dlrmf(rank),ipiv(rank)
      real *8 dlro(rank)
      complex *16 mf2cf(rank,rank)
      complex *16, external :: kfun

      integer j,k,info
      integer, allocatable :: ns(:)

      allocate(ns(2*nmax+1))

      ns = (/(j, j=-nmax,nmax)/)

      do k=1,rank
        do j=1,rank
          mf2cf(j,k) = kfun(dlrmf(j),dlro(k))
        enddo
      enddo

      call zgetrf(rank,rank,mf2cf,rank,ipiv,info)


      end subroutine dlr_mf2cf1




      subroutine dlr_convinit(rank,dlrrf,dlrit,phi)

      ! Get matrix of convolutions of DLR basis functions,
      ! evaluated on DLR grid
      !
      ! Fermionic case only
      !
      ! Entries will be given by
      !
      ! int_0^1 phi_j(t_i-t') phi_k(t') dt'
      !
      ! for t_i a DLR grid point, and phi_j, phi_k DLR basis functions. 

      implicit none
      integer rank
      real *8 dlrrf(rank),dlrit(rank)
      real *8 phi(rank,rank,rank)
      real *8, external :: kfun

      integer j,k,l,ier,maxrec,numint
      real *8 one,rint1,rint2
      real *8, external :: kfunf,kfunf2

      one = 1.0d0

      do l=1,rank
        do k=1,rank
          do j=1,rank

            if (k.ne.l) then

              phi(j,k,l) = (kfunf2(dlrit(j),dlrrf(l)) -&
                kfunf2(dlrit(j),dlrrf(k)))/(dlrrf(k)-dlrrf(l))

            else

              if (dlrit(j).gt.0.0d0) then

                phi(j,k,l) = (dlrit(j)-kfunf(1.0d0,dlrrf(k)))*&
                  kfunf2(dlrit(j),dlrrf(k))

              else

                phi(j,k,l) = (dlrit(j)+kfunf(0.0d0,dlrrf(k)))*&
                  kfunf2(dlrit(j),dlrrf(k))

              endif
            endif

          enddo
        enddo
      enddo

      end subroutine dlr_convinit


      subroutine dlr_conv(rank,phi,it2cf,ipiv,g,gmat)

      ! Build matrix of convolution by a Green's function in its DLR
      ! coefficient representation

      implicit none
      integer rank,ipiv(rank)
      real *8 phi(rank,rank,rank),it2cf(rank,rank),g(rank)
      real *8 gmat(rank,rank)

      integer i,j,info

      do j=1,rank
        do i=1,rank
          gmat(i,j) = sum(g*phi(i,:,j))
        enddo
      enddo

      call dgetrs('N',rank,rank,it2cf,rank,ipiv,gmat,rank,info)

      end subroutine dlr_conv


      subroutine dlr_dyson(rank,mu,dlrit,it2cf,ipiv,cf2it,phi,&
          g,numit,w,fptol,useg0,sigeval)

      ! Solve the Dyson equation self-consistently using the DLR in
      ! imaginary time and weighted fixed point iteration

      ! g: on input, initial guess, unless useg0==1. On output, DLR
      ! coefficients of solution
      !
      ! numit: on input, max number of fixed point iterations. On
      ! output, actual number of fixed point iterations.
      !
      ! sigeval: subroutine with calling sequence sigeval(rank,g,sig),
      ! which takes in DLR coefficients of a Green's function G, and
      ! returns values of Sigma on the DLR imaginary time grid.

      implicit none
      integer rank,numit,useg0,ipiv(rank)
      real *8 dlrit(rank),mu,it2cf(rank,rank)
      real *8 cf2it(rank,rank),g(rank),w,fptol,phi(rank,rank,rank)

      integer i,j,info
      integer, allocatable :: ipiv1(:)
      real *8 one
      real *8, allocatable :: g0(:),g0mat(:,:),sig(:),sigmat(:,:)
      real *8, allocatable :: sysmat(:,:),gnew(:)
      real *8, allocatable :: tmp(:)
      real *8, external :: kfunf2

      one = 1.0d0

      ! --- Get DLR of G0, matrix of convolution by G0 ---

      ! Imaginary time grid representation of G0

      allocate(g0(rank))

      do i=1,rank

        g0(i) = -kfunf2(dlrit(i),mu)

      enddo

      ! DLR coefficient representation of G0

      call dlr_expnd(rank,it2cf,ipiv,g0)

      ! Matrix of convolution by G0

      allocate(g0mat(rank,rank))

      call dlr_conv(rank,phi,it2cf,ipiv,g0,g0mat)


      ! --- Fixed point iteration for G ---

      ! Either use input initial guess, or use G0 as initial guess

      if (useg0==1) then
        g = g0
      endif

      allocate(sig(rank),sigmat(rank,rank),gnew(rank))
      allocate(sysmat(rank,rank),ipiv1(rank))

      do i=1,numit

        ! Get Sigma in imaginary time grid representation from DLR
        ! coefficient representation of previous G

        call sigeval(rank,g,sig)

        ! DLR coefficient representation of Sigma

        call dlr_expnd(rank,it2cf,ipiv,sig)

        ! Matrix of convolution by Sigma

        call dlr_conv(rank,phi,it2cf,ipiv,sig,sigmat)

        ! Linear VIE system matrix

        sysmat = -matmul(g0mat,sigmat)

        do j=1,rank
          sysmat(j,j) = one + sysmat(j,j)
        enddo

        ! Solve linear VIE

        call dgetrf(rank,rank,sysmat,rank,ipiv1,info)

        gnew = g0

        call dgetrs('N',rank,1,sysmat,rank,ipiv1,gnew,rank,info)

        ! Check self-consistency

        if (maxval(abs(matmul(cf2it,gnew-g)))<fptol) then

          g = gnew
          numit = i

          return

        else

          g = w*gnew + (one-w)*g

        endif
      enddo

      write(6,*) 'Warning: fixed point iteration did not converge.'

      end subroutine dlr_dyson

      subroutine dlr_dyson_mf(beta,rank,mu,dlrit,it2cf,ipiv1,cf2it,&
          dlrmf,mf2cf,ipiv2,cf2mf,g,numit,w,fptol,useg0,sigeval)

      ! Solve the Dyson equation self-consistently using the DLR in
      ! Marsubara frequency and weighted fixed point iteration

      ! g: on input, initial guess, unless useg0==1. On output, DLR
      ! coefficients of solution
      !
      ! numit: on input, max number of fixed point iterations. On
      ! output, actual number of fixed point iterations.
      !
      ! sigeval: subroutine with calling sequence sigeval(rank,g,sig),
      ! which takes in DLR coefficients of a Green's function G, and
      ! returns values of Sigma on the DLR imaginary time grid.

      implicit none
      integer rank,numit,useg0,ipiv1(rank),ipiv2(rank),dlrmf(rank)
      real *8 beta,dlrit(rank),mu,it2cf(rank,rank)
      real *8 cf2it(rank,rank),g(rank),w,fptol,phi(rank,rank,rank)
      complex *16 mf2cf(rank,rank),cf2mf(rank,rank)

      integer i,j,info
      real *8 one
      real *8, allocatable :: sig(:)
      complex *16, allocatable :: g0(:),gmf(:),sigmf(:),gnew(:)
      real *8, external :: kfunf2
      complex *16, external :: kfunf_mf

      one = 1.0d0

      ! --- Get Matsubara frequency representation of G0 ---

      allocate(g0(rank))

      do i=1,rank

        g0(i) = -kfunf_mf(dlrmf(i),beta*mu)

      enddo


      ! --- Fixed point iteration for G ---

      ! Either use input initial guess, or use G0 as initial guess

      allocate(gmf(rank))

      if (useg0==1) then
        
        gmf = g0

        call dlr_mfexpnd(rank,mf2cf,ipiv2,gmf)

        g = real(gmf)

        !do i=1,rank

        !  g(i) = -kfunf2(dlrit(i),beta*mu)

        !enddo

        !call dlr_expnd(rank,it2cf,ipiv1,g)

      endif

      allocate(sig(rank),gnew(rank),sigmf(rank))

      do i=1,numit

        ! Get Sigma in imaginary time grid representation from DLR
        ! coefficient representation of previous G

        call sigeval(rank,g,sig)

        ! DLR coefficient representation of Sigma

        call dlr_expnd(rank,it2cf,ipiv1,sig)

        ! Matsubara frequency representation of Sigma

        sigmf = matmul(cf2mf,sig)

        ! Solve linear VIE

        gnew = g0/(one-beta**2*g0*sigmf)

        ! DLR coefficient representation of solution

        call dlr_mfexpnd(rank,mf2cf,ipiv2,gnew)

        ! Check self-consistency

        if (maxval(abs(matmul(cf2it,real(gnew)-g)))<fptol) then

          g = real(gnew)
          numit = i

          return

        else

          g = w*real(gnew) + (one-w)*g

        endif
      enddo

      write(6,*) 'Warning: fixed point iteration did not converge.'

      end subroutine dlr_dyson_mf


      subroutine leg_dyson(n,beta,mu,xgl,wgl,legf,g,numit,w,fptol,&
          useg0,sigeval)

      ! Solve the Dyson equation self-consistently using a Legendre
      ! representation in imaginary time and weighted fixed point
      ! iteration

      ! g: on input, initial guess, unless useg0==1. On output, solution
      ! on Legendre grid
      !
      ! numit: on input, max number of fixed point iterations. On
      ! output, actual number of fixed point iterations.
      !
      ! sigeval: subroutine with calling sequence sigeval(n,g,sig),
      ! which takes in Green's function G on Legendre grid, and
      ! returns values of Sigma on that grid.

      implicit none
      integer n,numit,useg0
      real *8 beta,mu,xgl(n),wgl(n),legf(n,n)
      real *8 g(n),w,fptol

      integer i,j,info
      integer, allocatable :: ipiv(:)
      real *8 one
      real *8, allocatable :: g0(:),t(:),g0mat(:,:),sig(:),sigmat(:,:)
      real *8, allocatable :: sysmat(:,:),gnew(:)

      one = 1.0d0

      ! --- Get G0, matrix of convolution by G0 ---

      ! Imaginary time grid representation of G0

      allocate(g0(n),t(n))

      t = (xgl+one)/2*beta

      g0 = -exp(-mu*t)/(one+exp(-mu*beta))

      ! Matrix of convolution by G0

      allocate(g0mat(n,n))

      call leg_conv(beta,n,t,xgl,wgl,legf,g0,g0mat)


      ! --- Fixed point iteration for G ---

      ! Either use input initial guess, or use G0 as initial guess

      if (useg0==1) then
        g = g0
      endif

      allocate(sig(n),sigmat(n,n),gnew(n))
      allocate(sysmat(n,n),ipiv(n))

      do i=1,numit

        ! Get Sigma in Legendre grid representation from Legendre 
        ! coefficient representation of previous G

        call sigeval(n,g,sig)

        ! Matrix of convolution by Sigma

        call leg_conv(beta,n,t,xgl,wgl,legf,sig,sigmat)

        ! Linear VIE system matrix

        sysmat = -matmul(g0mat,sigmat)

        do j=1,n
          sysmat(j,j) = one + sysmat(j,j)
        enddo

        ! Solve linear VIE

        call dgetrf(n,n,sysmat,n,ipiv,info)

        gnew = g0

        call dgetrs('N',n,1,sysmat,n,ipiv,gnew,n,info)

        ! Check self-consistency

        if (maxval(abs(gnew-g))<fptol) then

          g = gnew
          numit = i

          return

        else

          g = w*gnew + (one-w)*g

        endif
      enddo

      write(6,*) 'Warning: fixed point iteration did not converge.'

      end subroutine leg_dyson


      subroutine leg_conv(beta,p,tt,xgl,wgl,legf,g,convmat)
  
      ! Build matrix of convolution by kernel defined on [0,beta] by
      ! function g sampled at Legendre nodes

      implicit none
      integer p
      real *8 beta,tt(p),xgl(p),wgl(p),legf(p,p),g(p),convmat(p,p)

      integer i,j
      real *8 one,ttar
      real *8, allocatable :: c(:)
      real *8, allocatable :: tgl(:),p1(:,:),p2(:,:),sig1(:),sig2(:)

      one = 1.0d0

      allocate(c(p))

      c = matmul(legf,g)

      allocate(tgl(p),p1(p,p),p2(p,p),sig1(p),sig2(p))
        
      do j=1,p ! Loop through target points
      
        ttar = tt(j); ! Target point t
      
        ! Evaluate Sigma(t-t') and Pn(t') at GL points on [0,t]
      
        tgl = (xgl+1)/2*ttar; ! GL points on [0,t]
      
        ! Evaluate Legendre polynomials on [0,beta] at GL points on [0,t]

        do i=1,p
     
          call legepols(2*tgl(i)/beta-one,p-1,p1(:,i))

        enddo
      
        tgl = ttar-tgl ! Input to Sigma
      
        ! Evaluate Sigma

        do i=1,p
          
          call legeexev(2*tgl(i)/beta-one,sig1(i),c,p-1)
     
        enddo
      
        ! Evaluate Sigma(t-t') and Pn(t') at GL points on [t,beta]
      
        tgl = (xgl+one)/2*(beta-ttar)+ttar; ! GL points on [t,beta]
      
        ! Evaluate Legendre polynomials on [0,beta] at GL points on
        ! [t,beta]

        do i=1,p

          call legepols(2*tgl(i)/beta-one,p-1,p2(:,i))

        enddo
      
        tgl = ttar-tgl+beta; ! Input to Sigma

        do i=1,p

          call legeexev(2*tgl(i)/beta-one,sig2(i),c,p-1)
          sig2(i) = -sig2(i)
        
        enddo
      
        ! Weight and sum
      
        sig1 = sig1*wgl*ttar/2
        sig2 = sig2*wgl*(beta-ttar)/2
      
        convmat(j,1:p) = matmul(p1,sig1) + matmul(p2,sig2)
      
      enddo
      
      convmat = matmul(convmat,legf)
          
      end subroutine leg_conv




      subroutine cgl_pt2coef(p,npan,nfun,legf,val,coef)

      ! ----- Convert point value representation on composite
      ! Gauss-Legendre grid to coefficient representation -----

      integer p,npan
      real *8 val(p*npan,nfun),coef(p*npan,nfun),legf(p,p)

      do i=1,npan
        coef((i-1)*p+1:i*p,:) = matmul(legf,val((i-1)*p+1:i*p,:))
      enddo

      end subroutine cgl_pt2coef


      subroutine ind_rearrange(n,krank,ind)

      ! Rearrange output ind from iddp_qrpiv or iddr_qrpiv to give list
      ! of krank columns selected by a pivoted QR process
      !
      ! Output overwrites first krank entries of ind

      implicit none
      integer n,krank,ind(n)

      integer k,iswap
      integer, allocatable :: tmp(:)

      allocate(tmp(n))

      do k = 1,n
        tmp(k) = k
      enddo
      
      do k = 1,krank
      
        ! Swap rnorms(k) and rnorms(list(k)).
      
        iswap = tmp(k)
        tmp(k) = tmp(ind(k))
        tmp(ind(k)) = iswap
      
      enddo

      ind(1:krank) = tmp(1:krank)
     
      end subroutine ind_rearrange
     

      subroutine evalexpfun(x,val)

      ! ----- Evaluate the function (1-e^(-x))/x -----

      implicit none
      real *8 x,val

      integer i
      real *8 one,x1,c

      one = 1.0d0

      if (x>1.0d-1) then

        val = (one-exp(-x))/x

      else

        val = one
        c = one
        x1 = x
        do i=1,10
          c = c*(i+1)
          val = val + ((-1)**i)*x1/c
          x1 = x1*x
        enddo

      endif

      end subroutine evalexpfun


      subroutine barychebinit(n,x,w)

      ! Get Chebyshev nodes of first kind and corresponding barycentric
      ! Lagrange interpolation weights

      implicit none
      integer n
      real *8 x(n),w(n)

      integer j
      real *8 pi,c

      pi = 4*atan(1.0d0)

      do j=1,n
        c = (2.0d0*j-1)/(2*n)*pi
        x(n-j+1) = cos(c)
        w(n-j+1) = (1-2*mod(j-1,2))*sin(c)
        !w(j) = (-1)**(j-1)*sin(c)
      enddo

      end subroutine barychebinit


      subroutine barycheb(n,x,f,wc,xc,val)

      ! Barycentric Lagrange interpolation at Chebyshev nodes

      implicit none
      integer n
      real *8 x,f(n),wc(n),xc(n),val

      integer j
      real *8 dif,q,num,den

      do j=1,n
        if (x==xc(j)) then
          val = f(j)
          return
        endif
      enddo

      num = 0.0d0
      den = 0.0d0
      do j=1,n

        dif = x-xc(j)
        q = wc(j)/dif
        num = num + q*f(j)
        den = den + q

      enddo
      
      val = num/den

      end subroutine barycheb
