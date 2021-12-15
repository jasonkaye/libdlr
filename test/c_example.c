#include <stdio.h>
#include "dlr_c/dlr_c.h"

int main() {
    int p, npt, npo, nt, no;
    double lambda = 100.;
    printf("DLR example usage from C.\n");
    printf("Settings for lambda = %2.2f\n", lambda);
    c_ccfine_init(&lambda, &p, &npt, &npo, &nt, &no);
    printf("p = %i, npt = %i, npo = %i, nt = %i, no = %i\n", p, npt, npo, nt, no);
    return 0;
}
