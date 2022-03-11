/* Header file for libdlr C interface */

void c_ccfine_init(double *lambda, int *p, int *npt, int *npo, int *nt, int *no);

void c_ccfine(double *lambda, int *p, int *npt, int *npo, double *t, double *om);

void c_dlr_kfine(double *lambda, int *p, int *npt, int *npo, double *t, double *om, double *kmat, double *err);

void c_dlr_rf(double *lambda, double *eps, int *nt, int *no, double *om, double *kmat, int *rank, double *dlrrf, int *oidx);

void c_dlr_it_build(double *lambda, double *eps, int *r, double* dlrrf, double* dlrit);

void c_dlr_it(double *lambda, int *nt, int *no, double *t, double *kmat, int *rank, int *oidx, double* dlrit, int *tidx);

void c_dlr_cf2it_init(int *rank, double *dlrrf, double *dlrit, double *cf2it);

void c_dlr_it2cf_init(int *rank, double *dlrrf, double *dlrit, double *it2cf, int *it2cfp);

void c_dlr_it2cf(int *r, double *it2cf, int *it2cfp, double *g, double *gc);

void c_dlr_it2itr_init(int *r, double *dlrrf, double *dlrit, double *it2cf, int *it2cfp, double *it2itr);

void c_dlr_it2itr(int *r, double *it2itr, double *g, double *gr);

void c_dlr_it_eval(int *r, double *dlrrf, double *gc, double *t, double *g);

void c_dlr_mf(int *nmax, int *rank, double *dlrrf, int *xi, int *dlrmf);

void c_dlr_cf2mf_init(int *rank, double *dlrrf,int *dlrmf, int *xi, double _Complex *cf2mf);

void c_dlr_mf2cf_init(int *nmax, int *rank, double *dlrrf,int *dlrmf, int *xi, double _Complex *dlrmf2cf, int *mf2cfpiv);

void c_dlr_convtens(double *beta, int *xi, int *r, double *dlrrf, double *dlrit, double *it2cf, int *it2cfp, double *phi);

void c_dlr_convmat(int *r, int*n, double *it2cf, int *it2cfp, double *phi, double *g, double *gmat);

void c_dlr_ipmat(double *beta, int *r, double *dlrit, double *dlrrf, double *it2cf, int *it2cfp, double *ipmat);

void c_dlr_dysonit(int *r, int *n, double *it2cf, int *it2cfp, double *phi, double *g0, double *g0mat, double *sig, double *g);

void c_eqpts_rel(int *n, double *t);

void c_kfunf_rel(double *t, double *om, double *val);
