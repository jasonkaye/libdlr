      !
      !
      ! This file contains subroutines for solving the Dyson equation
      ! using the discrete Lehmann representation 
      !
      !


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

      call dlr_convmat(rank,phi,it2cf,ipiv,g0,g0mat)


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

        call dlr_convmat(rank,phi,it2cf,ipiv,sig,sigmat)

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
          dlrmf,dlrmf2cf,ipiv2,cf2mf,g,numit,w,fptol,useg0,sigeval)

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
      complex *16 dlrmf2cf(rank,rank),cf2mf(rank,rank)

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

        call dlr_mfexpnd(rank,dlrmf2cf,ipiv2,gmf)

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

        call dlr_mfexpnd(rank,dlrmf2cf,ipiv2,gnew)

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


