/* Header file for libdlr C interface */

void c_ccfine_init(double *lambda, int *p, int *npt, int *npo, int *nt, int *no);

void c_ccfine(double *lambda, int *p, int *npt, int *npo, double *t, double *om);

void c_dlr_kfine(double *lambda, int *p, int *npt, int *npo, double *t, double *om, double *kmat, double *err);

void c_dlr_rf(double *lambda, double *eps, int *nt, int *no, double *om, double *kmat, int *r, double *dlrrf, int *oidx);

void c_dlr_it_build(double *lambda, double *eps, int *r, double* dlrrf, double* dlrit);

void c_dlr_it(double *lambda, int *nt, int *no, double *t, double *kmat, int *r, int *oidx, double* dlrit, int *tidx);

void c_dlr_cf2it_init(int *r, double *dlrrf, double *dlrit, double *cf2it);

void c_dlr_cf2it(int *r, int *n, double *cf2it, double *gc, double *g);

void c_dlr_it2cf_init(int *r, double *dlrrf, double *dlrit, double *it2cf, int *it2cfp);

void c_dlr_it2cf(int *r, int *n, double *it2cf, int *it2cfp, double *g, double *gc);

void c_dlr_it2itr_init(int *r, double *dlrrf, double *dlrit, double *it2cf, int *it2cfp, double *it2itr);

void c_dlr_it2itr(int *r, int *n, double *it2itr, double *g, double *gr);

void c_dlr_it_eval(int *r, int *n, double *dlrrf, double *gc, double *t, double *g);

void c_dlr_it_fit(int *r, int *n, double *dlrrf, int *m, double *tsamp, double *gsamp, double *gc);

void c_dlr_mf_build(double *lambda, double *eps, int *nmax, int *xi, int *r, double *dlrrf, int *dlrmf);

void c_dlr_mf(int *nmax, int *r, double *dlrrf, int *xi, int *dlrmf);

void c_dlr_cf2mf_init(int *r, double *dlrrf ,int *dlrmf, int *xi, double _Complex *cf2mf);

void c_dlr_cf2mf(int *r, int *n, double _Complex *cf2mf, double *gc, double _Complex *gn);

void c_dlr_mf2cf_init(int *nmax, int *r, double *dlrrf, int *dlrmf, int *xi, double _Complex *mf2cf, int *mf2cfp);

void c_dlr_mf2cf(int *r, int *n, double _Complex *mf2cf, int *mf2cfp, double _Complex *g, double *gc);

void c_dlr_mf_eval(int *r, int *n, double *dlrrf, int *xi, double *gc, int *nmf, double _Complex *gn);

void c_dlr_mf_fit(int *r, double *dlrrf, int *xi, int *m, int *nsamp, double _Complex *gsamp, double *gc);

void c_dlr_convtens(double *beta, int *xi, int *r, double *dlrrf, double *dlrit, double *it2cf, int *it2cfp, double *phi);

void c_dlr_convmat(int *r, int *n, double *it2cf, int *it2cfp, double *phi, double *g, double *gmat);

void c_dlr_conv(int *r,int *n,double *gmat,double *f,double *h);

void c_dlr_fstconv_init(double *beta, int *r, double *dlrrf, double *dlrit, double *cf2it, double *fstconv);

void c_dlr_fstconv(int *r, int *n, double *cf2it, double *it2cf, int *it2cfp, double *fstconv, double *f, double *g, double *h);

void c_dlr_ipmat(double *beta, int *r, double *dlrit, double *dlrrf, double *it2cf, int *it2cfp, double *ipmat);

void c_dyson_it(int *r, int *n, double *it2cf, int *it2cfp, double *phi, double *g0, double *g0mat, double *sig, double *g);

void c_dyson_mf(double *beta, int *r, int *n, double _Complex *g0, double _Complex *sigmf, double _Complex *gmf);

void c_eqpts_rel(int *n, double *t);

void c_rel2abs(int *n,double *t,double *tabs);

void c_abs2rel(int *n,double *tabs,double *t);

void c_kfunf_rel(double *t, double *om, double *val);

void c_kfunmf(int *n, double *om, double _Complex *val);
