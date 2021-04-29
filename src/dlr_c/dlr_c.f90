module libdlr_c
  use iso_c_binding, only: c_int, c_double, c_char
  implicit none
contains

  subroutine c_gridparams(lambda,p,npt,npo,nt,no) bind(C)
    integer(c_int), intent(out) :: p,npt,npo,nt,no
    real(c_double), intent(in) :: lambda
    call gridparams(lambda,p,npt,npo,nt,no)
  end subroutine c_gridparams

  subroutine c_kfine_cc(fb,lambda,p,npt,npo,t,om,kmat,err) bind(C)
    integer(c_int), intent(in) :: p, npt, npo
    real(c_double), intent(out) :: lambda, t(npt*p), om(2*npo*p)
    real(c_double), intent(out) :: kmat(2*npt*p,2*npo*p), err(2)
    character(c_char), intent(in) :: fb
    call kfine_cc(fb,lambda,p,npt,npo,t,om,kmat,err)
  end subroutine c_kfine_cc

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

  subroutine c_dlr_it2cf(nt,no,kmat,rank,oidx,tidx,dlrit2cf,it2cfpiv) bind(C)
    integer(c_int), intent(in) :: nt,no,rank,tidx(rank),oidx(rank)
    integer(c_int), intent(out) :: it2cfpiv(rank)
    real(c_double), intent(in) :: kmat(nt,no)
    real(c_double), intent(out) :: dlrit2cf(rank,rank)
    call dlr_it2cf(nt,no,kmat,rank,oidx,tidx,dlrit2cf,it2cfpiv)
  end subroutine c_dlr_it2cf

  subroutine c_dlr_mf(fb,nmax,rank,dlrrf,dlrmf) bind(C)
    integer(c_int), intent(in) :: nmax,rank
    integer(c_int), intent(out) :: dlrmf(rank)
    real(c_double), intent(in) :: dlrrf(rank)
    character(c_char), intent(in) :: fb
    call dlr_mf(fb,nmax,rank,dlrrf,dlrmf)
  end subroutine c_dlr_mf

  subroutine c_dlr_mf2cf(fb,nmax,rank,dlrrf,dlrmf,dlrmf2cf,mf2cfpiv) bind(C)
    integer(c_int), intent(in) :: nmax,rank,dlrmf(rank)
    integer(c_int), intent(out) :: mf2cfpiv(rank)
    real(c_double), intent(in) :: dlrrf(rank)
    !complex(c_double_complex), intent(out) :: dlrmf2cf(rank,rank)
    !complex(c_double), intent(out) :: dlrmf2cf(rank,rank)
    complex *16, intent(out) :: dlrmf2cf(rank,rank)
    character(c_char), intent(in) :: fb
    call dlr_mf2cf(fb,nmax,rank,dlrrf,dlrmf,dlrmf2cf,mf2cfpiv)
  end subroutine c_dlr_mf2cf

end module libdlr_c
