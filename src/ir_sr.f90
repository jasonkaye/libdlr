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


