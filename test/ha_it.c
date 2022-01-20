      /* -------------------------------------------------------------
       *
       * Copyright (C) 2021 The Simons Foundation
       * 
       * Author: Jason Kaye, Nils Wentzell
       * 
       * -------------------------------------------------------------
       * 
       * libdlr is licensed under the Apache License, Version 2.0 (the
       * "License"); you may not use this file except in compliance with
       * the License.  You may obtain a copy of the License at
       * 
       *     http://www.apache.org/licenses/LICENSE-2.0
       * 
       * Unless required by applicable law or agreed to in writing,
       * software distributed under the License is distributed on an "AS
       * IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
       * express or implied.  See the License for the specific language
       * governing permissions and limitations under the License.
       * 
       * ------------------------------------------------------------- */
      
      #include <stdlib.h>
      #include <stdio.h>
      #include <math.h>

      #include "dlr_c/dlr_c.h"

      void ha_it_main(double, double, int, double);
      double gfun(double, double);
      double max(double a, double b){ return a > b ? a : b; }

      int main() {

      // Test discrete Lehmann representation using Green's function
      // generated from Lehmann representation with density which is a
      // sum of two delta functions. Recover DLR coefficients from
      // values of Green's function on imaginary time grid, and measure
      // the error of the resulting DLR expansion on a test grid.
      
      int ntst;
      double lambda,eps,beta;

      // Input parameters

      lambda = 1000.0;   // DLR high energy cutoff
      eps = 1.0e-14;     // DLR error tolerance

      ntst = 10000;      // # imaginary time test points
      beta = 1000.0;     // Inverse temperature


      // Main test subroutine

      ha_it_main(lambda,eps,ntst,beta);
      
      }

      void ha_it_main(double lambda, double eps, int ntst, double beta) {
        
      int i,j,r;
      int *it2cfp;
      double one,gtrue,gtest,errl2,errlinf,gmax,gl2;
      double *ttst,*it2cf,*dlrit,*dlrrf;
      double *g,*gc;

      one = 1.0;


      // Get DLR frequencies, imaginary time grid

      r = 500; // Upper bound on DLR rank

      dlrrf = malloc(r*sizeof(double));
      dlrit = malloc(r*sizeof(double));

      c_dlr_it_build(&lambda,&eps,&r,dlrrf,dlrit);


      // Get imaginary time values -> DLR coefficients matrix (LU form)

      it2cfp = malloc(r*sizeof(int));
      it2cf  = malloc(r*r*sizeof(double));

      c_dlr_it2cf_init(&r,dlrrf,dlrit,it2cf,it2cfp);


      // Sample G at imaginary time nodes

      g = malloc(r*sizeof(double));

      for(i = 0; i < r; ++i) {

        g[i] = gfun(beta,dlrit[i]);

      }


      // Get DLR coefficients of G

      gc = malloc(r*sizeof(double));

      c_dlr_it2cf(&r,it2cf,it2cfp,g,gc);


      // Get test points in relative format

      ttst = malloc(ntst*sizeof(double));

      c_eqpts_rel(&ntst,ttst);


      // Measure L^inf and L^2 errors

      errlinf = 0*one;
      errl2 = 0*one;
      gmax = 0*one;
      gl2 = 0*one;

      for(i = 0; i < ntst; ++i) {

        // Evaluate Green's function

        gtrue = gfun(beta,ttst[i]);

        // Evaluate DLR

        c_dlr_it_eval(&r,dlrrf,gc,ttst+i,&gtest);

        // Update L^inf and L^2 errors, norms

        errlinf = max(errlinf,fabs(gtrue-gtest));
        errl2 = errl2 + (gtrue-gtest)*(gtrue-gtest);

        gmax = max(gmax,fabs(gtrue));
        gl2 = gl2 + gtrue*gtrue;

      }

      errl2 = sqrt((ttst[1]-ttst[0])*errl2);
      gl2 = sqrt((ttst[1]-ttst[0])*gl2);


      printf("\n");
      printf("DLR rank = %i\n",r);
      printf("Abs L^inf err = %2.2e\n",errlinf);
      printf("Abs L^2 err   = %2.2e\n",errl2);
      printf("Rel L^inf err = %2.2e\n",errlinf/gmax);
      printf("Rel L^2 err   = %2.2e\n",errl2/gl2);
      printf("\n");

      free(it2cfp);
      free(ttst); free(it2cf); free(dlrit);  free(dlrrf);
      free(g); free(gc);

      // Return failed status if error is not sufficiently small

      if (errlinf > 1.0e-13) {
        exit(1);
      }

      }



      double gfun(double beta, double t) {

      // Evaluate Green's function with sum-of-delta-functions spectral
      // density 

      double a1,a2,a3,a4,a5,tmp,val=0.;

      a1 = -0.804*beta;
      a2 = -0.443*beta;
      a3 =  0.093*beta;
      a4 =  0.915*beta;
      a5 =  0.929*beta;

      c_kfunf_rel(&t, &a1, &tmp); val += tmp;
      c_kfunf_rel(&t, &a2, &tmp); val += tmp;
      c_kfunf_rel(&t, &a3, &tmp); val += tmp;
      c_kfunf_rel(&t, &a4, &tmp); val += tmp;
      c_kfunf_rel(&t, &a5, &tmp); val += tmp;
      
      return val;
      }
