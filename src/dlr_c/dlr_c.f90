module libdlr_c
  use iso_c_binding, only: c_int, c_double, c_char
  implicit none
contains

  subroutine c_ccfine_init(lambda,p,npt,npo,nt,no) bind(C)
    integer(c_int), intent(out) :: p,npt,npo,nt,no
    real(c_double), intent(in) :: lambda
    call ccfine_init(lambda,p,npt,npo,nt,no)
  end subroutine c_ccfine_init

  subroutine c_ccfine(lambda,p,npt,npo,nt,no,t,om) bind(C)
    integer(c_int), intent(in) :: p,npt,npo,nt,no
    real(c_double), intent(in) :: lambda
    real(c_double), intent(out) :: t(npt*p), om(2*npo*p)
    call ccfine(lambda,p,npt,npo,nt,no,t,om)
  end subroutine c_ccfine
  
  subroutine c_dlr_kfine(lambda,p,npt,npo,t,om,kmat,err) bind(C)
    integer(c_int), intent(in) :: p, npt, npo
    real(c_double), intent(in) :: lambda, t(npt*p), om(2*npo*p)
    real(c_double), intent(out) :: kmat(2*npt*p,2*npo*p), err(2)
    call dlr_kfine(lambda,p,npt,npo,t,om,kmat,err)
  end subroutine c_dlr_kfine

  subroutine c_dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx) bind(C)
    integer(c_int), intent(in) :: nt,no
    integer(c_int), intent(inout) :: rank
    integer(c_int), intent(out) :: oidx(rank)
    real(c_double), intent(in) :: lambda,eps,om(no),kmat(nt,no)
    real(c_double), intent(out) :: dlrrf(rank)
    call dlr_rf(lambda,eps,nt,no,om,kmat,rank,dlrrf,oidx)
  end subroutine c_dlr_rf

  subroutine c_dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit,tidx) bind(C)
    integer(c_int), intent(in) :: nt,no,rank,oidx(rank)
    integer(c_int), intent(out) :: tidx(rank)
    real(c_double), intent(in) :: lambda,t(nt),kmat(nt,no)
    real(c_double), intent(out) :: dlrit(rank)
    call dlr_it(lambda,nt,no,t,kmat,rank,oidx,dlrit,tidx)
  end subroutine c_dlr_it

  subroutine c_dlr_cf2it_init(rank,dlrrf,dlrit,cf2it) bind(C)
    integer(c_int), intent(in) :: rank
    real(c_double), intent(in) :: dlrrf(rank),dlrit(rank)
    real(c_double), intent(out) :: cf2it(rank,rank)
    call dlr_cf2it_init(rank,dlrrf,dlrit,cf2it)
  end subroutine c_dlr_cf2it_init

  subroutine c_dlr_it2cf_init(rank,dlrrf,dlrit,dlrit2cf,it2cfpiv) bind(C)
    integer(c_int), intent(in) :: rank
    integer(c_int), intent(out) :: it2cfpiv(rank)
    real(c_double), intent(in) :: dlrrf(rank),dlrit(rank)
    real(c_double), intent(out) :: dlrit2cf(rank,rank)
    call dlr_it2cf_init(rank,dlrrf,dlrit,dlrit2cf,it2cfpiv)
  end subroutine c_dlr_it2cf_init

  subroutine c_dlr_mf(nmax,rank,dlrrf,xi,dlrmf) bind(C)
    integer(c_int), intent(in) :: nmax,rank,xi
    integer(c_int), intent(out) :: dlrmf(rank)
    real(c_double), intent(in) :: dlrrf(rank)
    call dlr_mf(nmax,rank,dlrrf,xi,dlrmf)
  end subroutine c_dlr_mf

  subroutine c_dlr_cf2mf_init(rank,dlrrf,dlrmf,xi,cf2mf) bind(C)
    integer(c_int), intent(in) :: rank,dlrmf(rank),xi
    real(c_double), intent(in) :: dlrrf(rank)
    !complex(c_double_complex), intent(out) :: cf2mf(rank,rank)
    complex *16, intent(out) :: cf2mf(rank,rank)
    call dlr_cf2mf_init(rank,dlrrf,dlrmf,xi,cf2mf)
  end subroutine c_dlr_cf2mf_init

  subroutine c_dlr_mf2cf_init(nmax,rank,dlrrf,dlrmf,xi,dlrmf2cf,mf2cfpiv) bind(C)
    integer(c_int), intent(in) :: nmax,rank,dlrmf(rank),xi
    integer(c_int), intent(out) :: mf2cfpiv(rank)
    real(c_double), intent(in) :: dlrrf(rank)
    !complex(c_double_complex), intent(out) :: dlrmf2cf(rank,rank)
    complex *16, intent(out) :: dlrmf2cf(rank,rank)
    call dlr_mf2cf_init(nmax,rank,dlrrf,dlrmf,xi,dlrmf2cf,mf2cfpiv)
  end subroutine c_dlr_mf2cf_init

end module libdlr_c
