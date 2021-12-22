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

  subroutine c_dlr_rf(lambda,eps,nt,no,om,kmat,r,dlrrf,oidx) bind(C)
    integer(c_int), intent(in) :: nt,no
    integer(c_int), intent(inout) :: r
    integer(c_int), intent(out) :: oidx(r)
    real(c_double), intent(in) :: lambda,eps,om(no),kmat(nt,no)
    real(c_double), intent(out) :: dlrrf(r)
    call dlr_rf(lambda,eps,nt,no,om,kmat,r,dlrrf,oidx)
  end subroutine c_dlr_rf

  subroutine c_dlr_it_build(lambda,eps,r,dlrrf,dlrit) bind(C)
    integer(c_int), intent(in) :: r
    real(c_double), intent(in) :: lambda,eps
    real(c_double), intent(out) :: dlrrf(r),dlrit(r)
    call dlr_it_build(lambda,eps,r,dlrrf,dlrit)
  end subroutine c_dlr_it_build

  subroutine c_dlr_it(lambda,nt,no,t,kmat,r,oidx,dlrit,tidx) bind(C)
    integer(c_int), intent(in) :: nt,no,r,oidx(r)
    integer(c_int), intent(out) :: tidx(r)
    real(c_double), intent(in) :: lambda,t(nt),kmat(nt,no)
    real(c_double), intent(out) :: dlrit(r)
    call dlr_it(lambda,nt,no,t,kmat,r,oidx,dlrit,tidx)
  end subroutine c_dlr_it

  subroutine c_dlr_cf2it_init(r,dlrrf,dlrit,cf2it) bind(C)
    integer(c_int), intent(in) :: r
    real(c_double), intent(in) :: dlrrf(r),dlrit(r)
    real(c_double), intent(out) :: cf2it(r,r)
    call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)
  end subroutine c_dlr_cf2it_init

  subroutine c_dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp) bind(C)
    integer(c_int), intent(in) :: r
    integer(c_int), intent(out) :: it2cfp(r)
    real(c_double), intent(in) :: dlrrf(r),dlrit(r)
    real(c_double), intent(out) :: it2cf(r,r)
    call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)
  end subroutine c_dlr_it2cf_init

  subroutine c_dlr_it2cf(r,it2cf,it2cfp,g,gc) bind(C)
    integer(c_int), intent(in) :: r,it2cfp(r)
    real(c_double), intent(in) :: it2cf(r,r),g(r)
    real(c_double), intent(out) :: gc(r)
    call dlr_it2cf(r,it2cf,it2cfp,g,gc)
  end subroutine c_dlr_it2cf

  subroutine c_dlr_it_eval(r,dlrrf,gc,t,g) bind(C)
    integer(c_int), intent(in) :: r
    real(c_double), intent(in) :: dlrrf(r),gc(r),t
    real(c_double), intent(out) :: g
    call dlr_it_eval(r,dlrrf,gc,t,g)
  end subroutine c_dlr_it_eval

  subroutine c_dlr_mf(nmax,r,dlrrf,xi,dlrmf) bind(C)
    integer(c_int), intent(in) :: nmax,r,xi
    integer(c_int), intent(out) :: dlrmf(r)
    real(c_double), intent(in) :: dlrrf(r)
    call dlr_mf(nmax,r,dlrrf,xi,dlrmf)
  end subroutine c_dlr_mf

  subroutine c_dlr_cf2mf_init(r,dlrrf,dlrmf,xi,cf2mf) bind(C)
    integer(c_int), intent(in) :: r,dlrmf(r),xi
    real(c_double), intent(in) :: dlrrf(r)
    !complex(c_double_complex), intent(out) :: cf2mf(r,r)
    complex *16, intent(out) :: cf2mf(r,r)
    call dlr_cf2mf_init(r,dlrrf,dlrmf,xi,cf2mf)
  end subroutine c_dlr_cf2mf_init

  subroutine c_dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,xi,dlrmf2cf,mf2cfpiv) bind(C)
    integer(c_int), intent(in) :: nmax,r,dlrmf(r),xi
    integer(c_int), intent(out) :: mf2cfpiv(r)
    real(c_double), intent(in) :: dlrrf(r)
    !complex(c_double_complex), intent(out) :: dlrmf2cf(r,r)
    complex *16, intent(out) :: dlrmf2cf(r,r)
    call dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,xi,dlrmf2cf,mf2cfpiv)
  end subroutine c_dlr_mf2cf_init

  subroutine c_dlr_convtens(beta,xi,r,dlrrf,dlrit,it2cf,it2cfp,phi) bind(C)
    integer(c_int), intent(in) :: xi,r,it2cfp(r)
    real(c_double), intent(in) :: beta,dlrrf(r),dlrit(r),it2cf(r,r)
    real(c_double), intent(out) :: phi(r*r,r)
    call dlr_convtens(beta,xi,r,dlrrf,dlrit,it2cf,it2cfp,phi)
  end subroutine c_dlr_convtens

  subroutine c_dlr_convmat(r,it2cf,it2cfp,phi,g,gmat) bind(C)
    integer(c_int), intent(in) :: r,it2cfp(r)
    real(c_double), intent(in) :: it2cf(r,r),phi(r*r,r),g(r)
    real(c_double), intent(out) :: gmat(r,r)
    call dlr_convmat(r,it2cf,it2cfp,phi,g,gmat)
  end subroutine c_dlr_convmat

  subroutine c_eqpts_rel(n,t) bind(C)
    integer(c_int), intent(in) :: n
    real(c_double), intent(out) :: t(n)
    call eqpts_rel(n,t)
  end subroutine c_eqpts_rel


end module libdlr_c
