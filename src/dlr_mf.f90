      ! -------------------------------------------------------------
      !
      ! This file contains subroutines to work with the DLR in Matsubara
      ! frequency
      !
      ! -------------------------------------------------------------
      !
      ! Copyright (C) 2021 The Simons Foundation
      ! 
      ! Author: Jason Kaye
      ! 
      ! -------------------------------------------------------------
      ! 
      ! libdlr is licensed under the Apache License, Version 2.0 (the
      ! "License"); you may not use this file except in compliance with
      ! the License.  You may obtain a copy of the License at
      ! 
      !     http://www.apache.org/licenses/LICENSE-2.0
      ! 
      ! Unless required by applicable law or agreed to in writing,
      ! software distributed under the License is distributed on an "AS
      ! IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
      ! express or implied.  See the License for the specific language
      ! governing permissions and limitations under the License.
      ! 
      ! -------------------------------------------------------------



      !> Get DLR frequency nodes and Matsubara frequency nodes
      !!
      !! Note: Matsubara frequencies are given by
      !! i*omega_n = i*pi*(2n+1) for fermionic Green's functions, and
      !! i*omega_n = i*pi*2n for bosonic Green's functions. The integers
      !! returned in the array dlrmf are the indices n of these
      !! frequencies. The input xi controls whether these are indices
      !! into fermionic or bosonic Matsubara frequencies. For example,
      !! if dlrmf(j) = 3 for some index j, and xi = -1, then the
      !! corresponding Matsubara frequency is i*omega_n = i*pi*(2*3+1); if
      !! xi = 1, then it is i*omega_n = i*pi*(2*3).
      !!
      !! @param[in]     lambda  dimensionless cutoff parameter
      !! @param[in]     eps     DLR error tolerance
      !! @param[in]     nmax    Matsubara frequency cutoff
      !! @param[in]     xi      xi=-1 for fermionic frequencies; xi=1
      !!                          for bosonic frequencies
      !! @param[in,out] r       On input, maximum possible number of DLR
      !!                          basis functions, defining input size
      !!                          of various arrays; on output, number
      !!                          of DLR basis functions.
      !! @param[out]    dlrrf   DLR frequency nodes
      !! @param[out]    dlrmf   DLR Matsubara frequency nodes

      subroutine dlr_mf_build(lambda,eps,nmax,xi,r,dlrrf,dlrmf)

      implicit none
      integer nmax,xi,r,dlrmf(r)
      real *8 lambda,eps,dlrrf(r)

      integer p,npt,npo,nt,no
      integer, allocatable :: oidx(:)
      real *8 kerr(2)
      real *8, allocatable :: kmat(:,:),t(:),om(:)


      ! Set parameters for fine grid

      call ccfine_init(lambda,p,npt,npo,nt,no)


      ! Get composite Chebyshev fine discretization of Lehmann kernel

      allocate(kmat(nt,no),t(nt),om(no))

      call ccfine(lambda,p,npt,npo,t,om)

      call dlr_kfine(lambda,p,npt,npo,t,om,kmat,kerr)


      ! Select DLR frequency nodes

      r = 500 ! Upper bound on possible DLR rank

      allocate(oidx(r))

      call dlr_rf(lambda,eps,nt,no,om,kmat,r,dlrrf,oidx)


      ! Get DLR Matsubara frequency nodes

      call dlr_mf(nmax,r,dlrrf,xi,dlrmf)

      end subroutine dlr_mf_build





      !> Get DLR Matsubara frequency nodes
      !! 
      !! Note: Matsubara frequencies are given by
      !! i*omega_n = i*pi*(2n+1) for fermionic Green's functions, and
      !! i*omega_n = i*pi*2n for bosonic Green's functions. The integers
      !! returned in the array dlrmf are the indices n of these
      !! frequencies. The input xi controls whether these are indices
      !! into fermionic or bosonic Matsubara frequencies. For example,
      !! if dlrmf(j) = 3 for some index j, and xi = -1, then the
      !! corresponding Matsubara frequency is i*omega_n = i*pi*(2*3+1); if
      !! xi = 1, then it is i*omega_n = i*pi*(2*3).
      !!
      !! @param[in]  nmax   Matsubara frequency cutoff
      !! @param[in]  r      number of DLR basis functions
      !! @param[in]  dlrrf  DLR frequency nodes
      !! @param[in]  xi     xi=-1 for fermionic frequencies; xi=1
      !!                      for bosonic frequencies
      !! @param[out] dlrmf  DLR Matsubara frequency nodes

      subroutine dlr_mf(nmax,r,dlrrf,xi,dlrmf)

      implicit none
      integer nmax,r,xi,dlrmf(r)
      real *8 dlrrf(r)

      integer i,k,info
      integer, allocatable :: ns(:),list(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: tmp(:,:)
      complex *16, external :: kfunmf

      ! Get matrix of Fourier transforms of DLR basis functions

      allocate(tmp(r,2*nmax+1),ns(2*nmax+1))

      ns = (/(i, i=-nmax,nmax)/)

      do i=1,2*nmax+1
        do k=1,r
          
          tmp(k,i) = kfunmf(2*ns(i)+(1-xi)/2,dlrrf(k))
          
        enddo
      enddo


      ! Pivoted QR to select Matsubara frequency nodes

      allocate(list(2*nmax+1),work(2*nmax+1))

      call idzr_qrpiv(r,2*nmax+1,tmp,r,list,work)


      ! Rearrange indices to get selected Matsubara frequency nodes
      ! indices

      call ind_rearrange(2*nmax+1,r,list)


      ! Extract selected Matsubara frequency nodes

      dlrmf = ns(list(1:r))

      end subroutine dlr_mf





      !> Build transform matrix from DLR coefficients to values of DLR
      !! expansion on Matsubara frequency grid
      !!
      !! Use the subroutine dlr_cf2mf to apply the transformation.
      !!
      !! @param[in]  r      number of DLR basis functions
      !! @param[in]  dlrrf  DLR freq nodes
      !! @param[in]  dlrmf  DLR Matsubara freq nodes
      !! @param[in]  xi      xi=-1 for fermionic frequencies; xi=1
      !!                       for bosonic frequencies
      !! @param[out] cf2mf  DLR coeffs -> Matsubara freq grid
      !!                      values transform matrix



      subroutine dlr_cf2mf_init(r,dlrrf,dlrmf,xi,cf2mf)

      implicit none
      integer r,dlrmf(r),xi
      real *8 dlrrf(r)
      complex *16 cf2mf(r,r)

      complex *16, external :: kfunmf

      integer i,j

      ! Get the matrix of Matsbuara frequency DLR basis functions at
      ! DLR Matsubara frequency nodes

      do j=1,r
        do i=1,r
          cf2mf(i,j) = kfunmf(2*dlrmf(i)+(1-xi)/2,dlrrf(j))
        enddo
      enddo

      end subroutine dlr_cf2mf_init





      !> Transform DLR coefficients to values of DLR expansion on
      !! Matsubara frequency grid
      !!
      !! @param[in]  r      number of DLR basis functions
      !! @param[in]  n      number of orbital indices
      !! @param[in]  cf2mf  DLR coefficients -> imaginary time grid
      !!                      values transform matrix
      !! @param[in]  beta   Inverse temperature
      !! @param[in]  gc     DLR coefficients of Green's function
      !! @param[out] gn     values of DLR expansion at Matsubara
      !!                      frequency grid points



      subroutine dlr_cf2mf(r,n,cf2mf,beta,gc,gn)

      implicit none
      integer r,n
      real *8 beta,gc(r,n,n)
      complex *16 cf2mf(r,r),gn(r,n,n)

      complex *16, allocatable :: tmp(:,:,:)

      ! Apply transformation matrix to coefficient vector

      tmp = beta*gc
      call zgemm('N','N',r,n*n,r,(1.0d0,0.0d0),cf2mf,r,tmp,r,(0.0d0,0.0d0),gn,r)

      end subroutine dlr_cf2mf





      !> Build transform matrix from values of a Green's function on
      !! Matsubara frequency grid to its DLR coefficients; matrix is
      !! stored in LU factored form
      !!
      !! Use the subroutine dlr_mf2cf to apply the transformation.
      !!
      !! @param[in]  nmax    Matsubara frequency cutoff
      !! @param[in]  r         number of DLR basis functions
      !! @param[in]  dlrrf   DLR frequency nodes
      !! @param[in]  dlrmf   DLR Matsubara frequency nodes
      !! @param[in]  xi      xi=-1 for fermionic frequencies; xi=1
      !!                       for bosonic frequencies
      !! @param[out] mf2cf   Matsubara frequency grid values ->
      !!                       DLR coefficients transform matrix,
      !!                       stored in LAPACK LU factored format; LU
      !!                       factors
      !! @param[out] mf2cfp  Matsubra frequency grid values ->
      !!                       DLR coefficients transform matrix,
      !!                       stored in LAPACK LU factored format; LU
      !!                       pivots

      subroutine dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,xi,mf2cf,mf2cfp)

      implicit none
      integer nmax,r,xi,dlrmf(r),mf2cfp(r)
      real *8 dlrrf(r)
      complex *16 mf2cf(r,r)

      integer i,j,info
      complex *16, external :: kfunmf

      ! Get the matrix of Matsbuara frequency DLR basis functions at
      ! DLR Matsubara frequency nodes

      call dlr_cf2mf_init(r,dlrrf,dlrmf,xi,mf2cf)


      ! LU factorize

      call zgetrf(r,r,mf2cf,r,mf2cfp,info)

      end subroutine dlr_mf2cf_init





      !> Transform values of DLR expansion on Matsubara frequency grid
      !! to DLR coefficients
      !!
      !! @param[in]  r       number of DLR basis functions
      !! @param[in]  n       number of orbital indices
      !! @param[in]  mf2cf   Matsubara frequency grid values ->
      !!                       DLR coefficients transform matrix,
      !!                       stored in LAPACK LU factored format; LU
      !!                       factors
      !! @param[in]  mf2cfp  Matsubra frequency grid values ->
      !!                       DLR coefficients transform matrix,
      !!                       stored in LAPACK LU factored format; LU
      !!                       pivots
      !! @param[in]  beta    Inverse temperature
      !! @param[in]  g       values of Green's function at Matsubara
      !!                       freq grid points
      !! @param[out] gc      DLR coefficients of Green's function

      subroutine dlr_mf2cf(r,n,mf2cf,mf2cfp,beta,g,gc)

      implicit none
      integer r,n,mf2cfp(r)
      real *8 gc(r,n,n),beta
      complex *16 mf2cf(r,r),g(r,n,n)

      integer info
      complex *16, allocatable :: tmp(:,:,:)

      ! Solve interpolation problem using DLR coefficients -> Matsubara
      ! frequency grid values matrix stored in LU form

      allocate(tmp(r,n,n))
      tmp = g/beta

      call zgetrs('N',r,n*n,mf2cf,r,mf2cfp,tmp,r,info)

      ! Remove imaginary part

      gc = real(tmp)

      end subroutine dlr_mf2cf





      !> Evaluate a DLR expansion at a Matsubara frequency point
      !! 
      !! Note: Matsubara frequencies are given by
      !! i*omega_n = i*pi*(2n+1) for fermionic Green's functions, and
      !! i*omega_n = i*pi*2n for bosonic Green's functions. The
      !! Matsubara frequency at which the DLR expansion is evaluated it
      !! determined by the input index nmf, and xi, which determines whether
      !! n is an index into a fermionic or bosonic Matsubara frequency.
      !! For example, if nmf=11 and xi = -1, then the DLR expansion will
      !! be evaluated at the Matsubara frequency 
      !! i*omega_n = i*pi*(2*11+1); if xi = 1, then it will be evaluated
      !! at i*omega_n = i*pi*(2*11).
      !!
      !! @param[in]  r      number of DLR basis functions
      !! @param[in]  n      number of orbital indices
      !! @param[in]  dlrrf  DLR frequency nodes
      !! @param[in]  xi     xi=-1 for fermionic frequencies; xi=1
      !!                      for bosonic frequencies
      !! @param[in]  gc     DLR coefficients of expansion
      !! @param[in]  nmf    Matsubara frequency integer
      !! @param[in]  beta   Inverse temperature
      !! @param[out] gn     value of DLR expansion at i*omega_n

      subroutine dlr_mf_eval(r,n,dlrrf,xi,gc,nmf,beta,gn)

      implicit none
      integer r,n,xi,nmf
      real *8 dlrrf(r),gc(r,n,n),beta
      complex *16 gn(n,n)

      integer i
      complex *16 kval
      complex *16, external :: kfunmf

      ! Evaluate Matsubara frequency DLR basis functions and sum against
      ! DLR coefficients

      gn = (0.0d0,0.0d0)
      do i=1,r

        kval = kfunmf(2*nmf+(1-xi)/2,dlrrf(i))

        gn = gn + gc(i,:,:)*kval

      enddo

      gn = beta*gn

      end subroutine dlr_mf_eval





      !> Get DLR coefficients from scattered data by least squares
      !! fitting
      !!
      !! @param[in]  r        number of DLR basis functions
      !! @param[in]  dlrrf    DLR frequency nodes
      !! @param[in]  xi       xi=-1 for fermionic frequencies; xi=1 for
      !!                        bosonic frequencies
      !! @param[in]  m        number of Matsubara frequency points at
      !!                        which Green's function G is sampled
      !! @param[in]  nsamp    Matsubara frequency integers at which G is
      !!                        sampled
      !! @param[in]  gsamp    values of G at sampling points
      !! @param[out] gc       DLR coefficients of Green's function

      subroutine dlr_mf_fit(r,dlrrf,xi,m,nsamp,gsamp,gc)

      implicit none
      integer r,xi,m,nsamp(m)
      real *8 dlrrf(r),gc(r)
      complex *16 gsamp(m)

      integer i,j,rank,lwork,info
      real *8 rcond
      integer, allocatable :: jpvt(:)
      real *8, allocatable :: rwork(:)
      complex *16, allocatable :: kls(:,:),work(:),tmp(:)
      complex *16, external :: kfunmf

      
      ! Get system matrix for least squares fitting; columns are DLR
      ! basis functions evaluated at Matsubara frequency sampling points

      allocate(kls(m,r))

      do j=1,r
        do i=1,m
          kls(i,j) = kfunmf(2*nsamp(i)+(1-xi)/2,dlrrf(j))
        enddo
      enddo


      ! Get size of work array for least squares fitting

      allocate(work(1),jpvt(r),tmp(m),rwork(2*r))

      call zgelsy(m,r,1,kls,m,tmp,m,jpvt,rcond,rank,work,-1,rwork,info)

      lwork = work(1)

      deallocate(work)


      ! Least squares fitting of data to determine DLR coefficients

      allocate(work(lwork))

      tmp = gsamp

      call zgelsy(m,r,1,kls,m,tmp,m,jpvt,rcond,rank,work,lwork,rwork,&
        info)

      gc = real(tmp(1:r))

      end subroutine dlr_mf_fit
