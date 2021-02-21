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


      
      subroutine irbasis(lambda,eps,nt,no,p,kmat,twgt,rank,irb)

      ! ----- Compute intermediate representation basis -----

      ! On input, rank is maximum possible rank and defines the sizes of
      ! tpts, ompts, and dlrmat. On output, rank is the actual rank.

      implicit none
      integer nt,no,p,rank,dosvd
      real *8 lambda,eps,kmat(2*p*nt,2*p*no),twgt(2*p*nt)
      real *8 irb(2*p*nt,rank)

      integer nnt,nno,i,j,iwhich(6),info,lwork
      integer, allocatable :: irows(:),icols(:)
      real *8 one,a,b,errout,start,finish,junk
      real *8 junk1,junk2,junk3,junk4
      real *8, allocatable :: xgl(:),legf(:,:),legb(:,:),wgl(:)
      real *8, allocatable :: pbpt(:),pbpo(:)
      real *8, allocatable :: ktcoef(:,:),komcoef(:,:)
      real *8, allocatable :: err1(:,:),err2(:,:)
      real *8, allocatable :: work(:),u(:,:),s(:)

      one = 1.0d0

      ! --- Derived parameters ---

      nnt = 2*nt*p
      nno = 2*no*p

      ! --- Reweight K(tau,omega) matrix in order to construct singular
      ! vectors which are orthonormal w.r.t. quadrature weights
      ! (approximate orthonormality in L2 rather than l2) ---

      !do i=1,2*nt*p
      !  kmat(i,:) = kmat(i,:)*sqrt(twgt(i))
      !enddo

      ! --- Compute SVD ---

      call cpu_time(start)

      lwork = 10*max(nnt,nno)
      allocate(s(min(nnt,nno)),u(nnt,min(nnt,nno)),work(lwork))

      call dgesvd('S','N',nnt,nno,kmat,nnt,s,u,nnt,junk,1,&
        work,lwork,info)

      ! --- Truncate SVD to get IR basis ---

      rank = min(nnt,nno)
      do i=1,min(nnt,nno)
        if (s(i).lt.eps*s(1)) then
          rank = i-1
          exit
        endif
      enddo

      irb(:,1:rank) = u(:,1:rank)

      call cpu_time(finish)
     
      write(6,*) '---------------------------------------------------'
      write(6,*) 'Built IR basis.'
      write(6,*) 'IR rank = ',rank
      write(6,*) 'Time for SVD = ',finish-start
      write(6,*) '---------------------------------------------------'
      write(6,*) ''

      ! --- Unweight singular vectors ---

      !do i=1,2*nt*p
      !  irb(i,:) = irb(i,:)/sqrt(twgt(i))
      !enddo

      end subroutine irbasis


      

      subroutine irskel(nt,no,p,t,twgt,rank,irb,tpts,irmat)

      ! ----- Compute intermediate representation basis -----

      ! On input, rank is maximum possible rank and defines the sizes of
      ! tpts, ompts, and dlrmat. On output, rank is the actual rank.

      !!! THIS ROUTINE USES ID_DIST !!!

      implicit none
      integer nt,no,p,rank
      real *8 t(2*nt*p),twgt(2*nt*p),irb(2*p*nt,rank),tpts(rank)
      real *8 irmat(rank,rank)

      integer nnt,nno,i,j
      integer, allocatable :: list(:)
      real *8 start,finish
      real *8, allocatable :: tmp(:,:),rnorms(:),pmat(:,:)

      ! --- Derived parameters ---

      nnt = 2*nt*p
      nno = 2*no*p
      
      ! ----- Compute interpolative decomposition and get IR grid -----

      allocate(tmp(rank,nnt),list(nnt),rnorms(nnt))
      allocate(pmat(rank,nnt))

      call cpu_time(start)

      tmp = transpose(irb)

      call iddr_id(rank,nnt,tmp,rank,list,rnorms)

      tpts = t(list(1:rank))

      call cpu_time(finish)


      write(6,*) '---------------------------------------------------'
      write(6,*) 'Built IR grid.'
      write(6,*) 'Time for ID = ',finish-start
      write(6,*) '---------------------------------------------------'
      write(6,*) ''

      ! --- Get coarse grid pts -> IR coefs matrix ---

      call idd_reconint(nnt,list,rank,tmp,pmat)

      irmat = transpose(matmul(pmat,irb))

      end subroutine irskel



      subroutine ireval(p,n,neval,pbp,irb_coef,c,x,y)

      ! ----- Evaluate IR expansion at a point -----

      implicit none
      integer p,n,neval
      real *8 pbp(2*n+1),irb_coef(2*n*p,neval),c(neval),x,y

      integer i,pan
      real *8 one,xx,a,b
      real *8, allocatable :: vals(:)

      one = 1.0d0

      ! --- Figure out what panel the point is in and get panel
      ! endpoints ---

      if (x<0.5d0) then
        
        xx = 2*x
        xx = log(xx)/log(2.0d0)
        pan = max(ceiling(n+xx),1)

      else

        xx = 1-x
        xx = 2*xx
        xx = log(xx)/log(2.0d0)
        pan = 2*n-max(ceiling(n+xx),1)+1

      endif

      a = pbp(pan)
      b = pbp(pan+1)

      ! --- Get panel endpoints and transform point such that panel
      ! becomes [-1,1] ---

      xx = -one + 2*(x-a)/(b-a)

      ! --- Evaluate IR basis functions ---

      !!! TODO: THIS CAN BE SPED UP BY PRECOMPUTING LEGENDRE POLYNOMIALS

      allocate(vals(neval))

      do i=1,neval
      
        call legeexev(xx,vals(i),irb_coef((pan-1)*p+1:pan*p,i),p-1)

      enddo

      ! --- Sum up expansion ---

      y = sum(c*vals)

      end subroutine ireval


      subroutine ireval1(p,n,pbp,irbc,x,y)

      ! ----- Evaluate single IR basis function at a point -----

      implicit none
      integer p,n,neval
      real *8 pbp(2*n+1),irbc(2*n*p),x,y

      integer i,pan
      real *8 one,xx,a,b
      real *8, allocatable :: vals(:)

      one = 1.0d0

      ! --- Figure out what panel the point is in and get panel
      ! endpoints ---

      if (x<0.5d0) then
        
        xx = 2*x
        xx = log(xx)/log(2.0d0)
        pan = max(ceiling(n+xx),1)

      else

        xx = 1-x
        xx = 2*xx
        xx = log(xx)/log(2.0d0)
        pan = 2*n-max(ceiling(n+xx),1)+1

      endif

      a = pbp(pan)
      b = pbp(pan+1)

      ! --- Get panel endpoints and transform point such that panel
      ! becomes [-1,1] ---

      xx = -one + 2*(x-a)/(b-a)

      ! --- Evaluate IR basis function ---

      !!! TODO: THIS CAN BE SPED UP BY PRECOMPUTING LEGENDRE POLYNOMIALS

      call legeexev(xx,y,irbc((pan-1)*p+1:pan*p),p-1)

      end subroutine ireval1



      
      subroutine irexpand(rank,irmat,g)

      ! ----- Get coefficients of IR expansion from samples on coarse
      ! grid -----

      ! On input, g contains samples of function at coarse grid points.
      ! On output, it contains the coefficients of the IR expansion.

      implicit none
      integer rank
      real *8 irmat(rank,rank),g(rank)

      g = matmul(irmat,g)

      end subroutine irexpand



      subroutine irconv(p,nt,pbp,rank,irbc,irt,phi)

      ! ----- Get matrix of convolutions of IR basis functions,
      ! evaluated on IR grid -----

      !
      ! Entries will be given by
      !
      ! int_0^1 phi_j(t_i-t') phi_k(t') dt'
      !
      ! for t_i an IR grid points, and phi_j, phi_k IR basis functions. 

      implicit none
      integer p,nt,rank
      real *8 pbp(2*nt+1),irbc(2*p*nt,rank),irt(rank)
      real *8 phi(rank,rank,rank)

      integer i,j,k,ier,maxrec,numint
      real *8 one,par1,par2,rint1,rint2

      one = 1.0d0

      do k=1,rank
        write(6,*) k
        do j=1,rank
          do i=1,rank

            call adapgaus(ier,0*one,irt(i),intgrd1,par1,par2,&
              16,1.0d-14,rint1,maxrec,numint)

            if (ier.ne.0) then

              write(6,*) 'ier = ',ier

            endif

            call adapgaus(ier,irt(i),one,intgrd2,par1,par2,&
              16,1.0d-14,rint2,maxrec,numint)

            phi(i,j,k) = rint1-rint2

          enddo
        enddo
      enddo

      contains

        real *8 function intgrd1(tp,par1,par2)

        implicit none
        real *8 tp,par1,par2

        integer ier,maxrec,numint
        real *8 val1,val2

        call ireval1(p,nt,pbp,irbc(:,j),irt(i)-tp,val1)
        
        call ireval1(p,nt,pbp,irbc(:,k),tp,val2)

        intgrd1 = val1*val2

        end function intgrd1


        real *8 function intgrd2(tp,par1,par2)

        implicit none
        real *8 tp,par1,par2

        integer ier,maxrec,numint
        real *8 val1,val2

        call ireval1(p,nt,pbp,irbc(:,j),one+irt(i)-tp,val1)
        
        call ireval1(p,nt,pbp,irbc(:,k),tp,val2)

        intgrd2 = val1*val2

        end function intgrd2

      end subroutine irconv





      subroutine irpt2coef(p,nt,rank,irb,irbc)

      ! ----- Convert point values representation of IR basis functions
      ! to composite Gauss-Legendre coefficients representation -----

      implicit none
      integer p,nt,rank
      real *8 irb(2*p*nt,rank),irbc(2*p*nt,rank)

      real *8, allocatable :: xgl(:),legf(:,:),legb(:,:),wgl(:)

      allocate(xgl(p),legf(p,p),legb(p,p),wgl(p))
      
      call legeexps(2,p,xgl,legf,legb,wgl)

      call cgl_pt2coef(p,2*nt,rank,legf,irb,irbc)

      end subroutine irpt2coef


