module libdlr_c
  use iso_c_binding, only: c_int, c_double, c_double_complex, c_char
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

  subroutine c_dlr_cf2it(r,n,cf2it,gc,g) bind(C)
    integer(c_int), intent(in) :: r,n
    real(c_double), intent(in) :: cf2it(r,r),gc(r,n,n)
    real(c_double), intent(out) :: g(r,n,n)
    call dlr_cf2it(r,n,cf2it,gc,g)
  end subroutine c_dlr_cf2it

  subroutine c_dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp) bind(C)
    integer(c_int), intent(in) :: r
    integer(c_int), intent(out) :: it2cfp(r)
    real(c_double), intent(in) :: dlrrf(r),dlrit(r)
    real(c_double), intent(out) :: it2cf(r,r)
    call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)
  end subroutine c_dlr_it2cf_init

  subroutine c_dlr_it2cf(r,n,it2cf,it2cfp,g,gc) bind(C)
    integer(c_int), intent(in) :: r,n,it2cfp(r)
    real(c_double), intent(in) :: it2cf(r,r),g(r,n,n)
    real(c_double), intent(out) :: gc(r,n,n)
    call dlr_it2cf(r,n,it2cf,it2cfp,g,gc)
  end subroutine c_dlr_it2cf

  subroutine c_dlr_it2itr_init(r,dlrrf,dlrit,it2cf,it2cfp,it2itr) bind(C)
    integer(c_int), intent(in) :: r,it2cfp(r)
    real(c_double), intent(in) :: dlrrf(r),dlrit(r),it2cf(r,r)
    real(c_double), intent(out) :: it2itr(r,r)
    call dlr_it2itr_init(r,dlrrf,dlrit,it2cf,it2cfp,it2itr)
  end subroutine c_dlr_it2itr_init

  subroutine c_dlr_it2itr(r,n,it2itr,g,gr) bind(C)
    integer(c_int), intent(in) :: r,n
    real(c_double), intent(in) :: it2itr(r,r),g(r,n,n)
    real(c_double), intent(out) :: gr(r,n,n)
    call dlr_it2itr(r,n,it2itr,g,gr)
  end subroutine c_dlr_it2itr

  subroutine c_dlr_it_eval(r,n,dlrrf,gc,t,gt) bind(C)
    integer(c_int), intent(in) :: r,n
    real(c_double), intent(in) :: dlrrf(r),gc(r,n,n),t
    real(c_double), intent(out) :: gt(n,n)
    call dlr_it_eval(r,n,dlrrf,gc,t,gt)
  end subroutine c_dlr_it_eval

  subroutine c_dlr_it_fit(r,n,dlrrf,m,tsamp,gsamp,gc) bind(C)
    integer(c_int), intent(in) :: r,n,m
    real(c_double), intent(in) :: dlrrf(r),tsamp(m),gsamp(m,n,n)
    real(c_double), intent(out) :: gc(r,n,n)
    call dlr_it_fit(r,n,dlrrf,m,tsamp,gsamp,gc)
  end subroutine c_dlr_it_fit

  subroutine c_dlr_mf_build(lambda,eps,nmax,xi,r,dlrrf,dlrmf) bind(C)
    integer(c_int), intent(in) :: nmax,xi,r
    real(c_double), intent(in) :: lambda,eps,dlrrf(r)
    integer(c_int), intent(out) :: dlrmf(r)
    call dlr_mf_build(lambda,eps,nmax,xi,r,dlrrf,dlrmf)
  end subroutine c_dlr_mf_build

  subroutine c_dlr_mf(nmax,r,dlrrf,xi,dlrmf) bind(C)
    integer(c_int), intent(in) :: nmax,r,xi
    integer(c_int), intent(out) :: dlrmf(r)
    real(c_double), intent(in) :: dlrrf(r)
    call dlr_mf(nmax,r,dlrrf,xi,dlrmf)
  end subroutine c_dlr_mf

  subroutine c_dlr_cf2mf_init(r,dlrrf,dlrmf,xi,cf2mf) bind(C)
    integer(c_int), intent(in) :: r,dlrmf(r),xi
    real(c_double), intent(in) :: dlrrf(r)
    complex(c_double_complex), intent(out) :: cf2mf(r,r)
    call dlr_cf2mf_init(r,dlrrf,dlrmf,xi,cf2mf)
  end subroutine c_dlr_cf2mf_init

  subroutine c_dlr_cf2mf(r,n,cf2mf,gc,gn) bind(C)
    integer(c_int), intent(in) :: r,n
    real(c_double), intent(in) :: gc(r,n,n)
    complex(c_double_complex), intent(in) :: cf2mf(r,r)
    complex(c_double_complex), intent(out) :: gn(r,n,n)
    call dlr_cf2mf(r,n,cf2mf,gc,gn)
  end subroutine c_dlr_cf2mf

  subroutine c_dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,xi,mf2cf,mf2cfp) bind(C)
    integer(c_int), intent(in) :: nmax,r,dlrmf(r),xi
    integer(c_int), intent(out) :: mf2cfp(r)
    real(c_double), intent(in) :: dlrrf(r)
    complex(c_double_complex), intent(out) :: mf2cf(r,r)
    call dlr_mf2cf_init(nmax,r,dlrrf,dlrmf,xi,mf2cf,mf2cfp)
  end subroutine c_dlr_mf2cf_init

  subroutine c_dlr_mf2cf(r,n,mf2cf,mf2cfp,g,gc) bind(C)
    integer(c_int), intent(in) :: r,n,mf2cfp(r)
    complex(c_double_complex), intent(in) :: mf2cf(r,r),g(r,n,n)
    real(c_double), intent(out) :: gc(r,n,n)
    call dlr_mf2cf(r,n,mf2cf,mf2cfp,g,gc)
  end subroutine c_dlr_mf2cf

  subroutine c_dlr_mf_eval(r,n,dlrrf,xi,gc,nmf,gn) bind(C)
    integer(c_int), intent(in) :: r,n,xi,nmf
    real(c_double), intent(in) :: dlrrf(r),gc(r,n,n)
    complex(c_double_complex), intent(out) :: gn(n,n)
    call dlr_mf_eval(r,n,dlrrf,xi,gc,nmf,gn)
  end subroutine c_dlr_mf_eval

  subroutine c_dlr_mf_fit(r,dlrrf,xi,m,nsamp,gsamp,gc) bind(C)
    integer(c_int), intent(in) :: r,xi,m,nsamp(m)
    real(c_double), intent(in) :: dlrrf(r)
    complex(c_double_complex), intent(in) :: gsamp(m)
    real(c_double), intent(out) :: gc(r)
    call dlr_mf_fit(r,dlrrf,xi,m,nsamp,gsamp,gc)
  end subroutine c_dlr_mf_fit

  subroutine c_dlr_convtens(beta,xi,r,dlrrf,dlrit,it2cf,it2cfp,phi) bind(C)
    integer(c_int), intent(in) :: xi,r,it2cfp(r)
    real(c_double), intent(in) :: beta,dlrrf(r),dlrit(r),it2cf(r,r)
    real(c_double), intent(out) :: phi(r*r,r)
    call dlr_convtens(beta,xi,r,dlrrf,dlrit,it2cf,it2cfp,phi)
  end subroutine c_dlr_convtens

  subroutine c_dlr_convmat(r,n,it2cf,it2cfp,phi,g,gmat) bind(C)
    integer(c_int), intent(in) :: r,n,it2cfp(r)
    real(c_double), intent(in) :: it2cf(r,r),phi(r*r,r),g(r,n,n)
    real(c_double), intent(out) :: gmat(r*n,r*n)
    call dlr_convmat(r,n,it2cf,it2cfp,phi,g,gmat)
  end subroutine c_dlr_convmat

  subroutine c_dlr_conv(r,n,gmat,f,h) bind(C)
    integer(c_int), intent(in) :: r,n
    real(c_double), intent(in) :: gmat(r*n,r*n),f(r,n,n)
    real(c_double), intent(out) :: h(r,n,n)
    call dlr_conv(r,n,gmat,f,h)
  end subroutine c_dlr_conv

  subroutine c_dlr_fstconv_init(beta,r,dlrrf,dlrit,cf2it,fstconv) bind(C)
    integer(c_int), intent(in) :: r
    real(c_double), intent(in) :: beta,dlrrf(r),dlrit(r),cf2it(r,r)
    real(c_double), intent(out) :: fstconv(r,2*r)
    call dlr_fstconv_init(beta,r,dlrrf,dlrit,cf2it,fstconv)
  end subroutine c_dlr_fstconv_init

  subroutine c_dlr_fstconv(r,n,cf2it,it2cf,it2cfp,fstconv,f,g,h) bind(C)
    integer(c_int), intent(in) :: r,n,it2cfp(r)
    real(c_double), intent(in) :: cf2it(r,r),it2cf(r,r),fstconv(r,2*r)
    real(c_double), intent(in) :: f(r,n,n),g(r,n,n)
    real(c_double), intent(out) :: h(r,n,n)
    call dlr_fstconv(r,n,cf2it,it2cf,it2cfp,fstconv,f,g,h)
  end subroutine c_dlr_fstconv

  subroutine c_dlr_ipmat(beta,r,dlrit,dlrrf,it2cf,it2cfp,ipmat) bind(C)
    integer(c_int), intent(in) :: r,it2cfp(r)
    real(c_double), intent(in) :: beta,dlrit(r),dlrrf(r),it2cf(r,r)
    real(c_double), intent(out) :: ipmat(r,r)
    call dlr_ipmat(beta,r,dlrit,dlrrf,it2cf,it2cfp,ipmat)
  end subroutine c_dlr_ipmat

  subroutine c_dyson_it(r,n,it2cf,it2cfp,phi,g0,g0mat,sig,g) bind(C)
    integer(c_int), intent(in) :: r,n,it2cfp(r)
    real(c_double), intent(in) :: it2cf(r,r)
    real(c_double), intent(in) :: phi(r*r,r),g0(r,n,n),g0mat(r*n,r*n)
    real(c_double), intent(in) :: sig(r,n,n)
    real(c_double), intent(out) :: g(r,n,n)
    call dyson_it(r,n,it2cf,it2cfp,phi,g0,g0mat,sig,g)
  end subroutine c_dyson_it

  subroutine c_dyson_mf(beta,r,n,g0,sigmf,gmf) bind(C)
    integer(c_int), intent(in) :: r,n
    real(c_double), intent(in) :: beta
    complex(c_double_complex), intent(in) :: g0(r,n,n),sigmf(r,n,n)
    complex(c_double_complex), intent(out) :: gmf(r,n,n)
    call dyson_mf(beta,r,n,g0,sigmf,gmf)
  end subroutine c_dyson_mf

  subroutine c_eqpts_rel(n,t) bind(C)
    integer(c_int), intent(in) :: n
    real(c_double), intent(out) :: t(n)
    call eqpts_rel(n,t)
  end subroutine c_eqpts_rel

  subroutine c_rel2abs(n,t,tabs) bind(C)
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: t(n)
    real(c_double), intent(out) :: tabs(n)
    call rel2abs(n,t,tabs)
  end subroutine c_rel2abs

  subroutine c_abs2rel(n,tabs,t) bind(C)
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: tabs(n)
    real(c_double), intent(out) :: t(n)
    call abs2rel(n,tabs,t)
  end subroutine c_abs2rel

  subroutine c_kfunf_rel(t,om,val) bind(C)
    real(c_double), intent(in) :: t,om
    real(c_double), intent(out) :: val

    real *8, external :: kfunf_rel

    val = kfunf_rel(t,om)

  end subroutine c_kfunf_rel

  subroutine c_kfunmf(n,om,val) bind(C)
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: om
    complex(c_double_complex), intent(out) :: val

    complex(c_double_complex), external :: kfunmf

    val = kfunmf(n,om)

  end subroutine c_kfunmf

end module libdlr_c
