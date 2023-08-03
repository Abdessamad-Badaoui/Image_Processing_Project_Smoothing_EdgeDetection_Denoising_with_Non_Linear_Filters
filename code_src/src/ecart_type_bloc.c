#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

static int doubleComparator(const void *first, const void *second) {
    double firstDouble = * (const double *) first;
    double secondDouble = * (const double *) second;
    if (firstDouble < secondDouble) {
        return -1;
    }
    else if (firstDouble > secondDouble) {
        return 1;
    }
    else {
        return 0;
    }
}

double *ecart_type_bloc(double **entree, int t, int nl, int nc) {

    double *sortie = calloc(nc*nl, sizeof(double));
    assert(sortie != NULL);

    int x, y, p, q;

    for (x = 0; x < nl; x++) {
        for (y = 0; y < nc; y++) {
            double ecart_type_xy = 0;
            double esperance_xy = 0;
            for (p = -t; p <= t; p++) {
                int indice_x = x - p;
                if ((x - p) < 0) {
                    indice_x = abs(x - p) -1;
                }
                else if ((x - p) >= nl) {
                    indice_x = 2*nl - (x - p) -1;
                }
                for (int q = -t; q <= t; q++) {
                    int indice_y = y - q;
                    if ((y - q) < 0) {
                        indice_y = abs(y - q) - 1;
                    }
                    else if ((y - q) >= nc) {
                        indice_y = 2*nc - (y - q) -1;
                    }

                    esperance_xy += entree[indice_x][indice_y];
                }
            }

            esperance_xy /= (2*t + 1)*(2*t + 1);

            for (p = -t; p <= t; p++) {
                int indice_x = x - p;
                if ((x - p) < 0) {
                    indice_x = abs(x - p) -1;
                }
                else if ((x - p) >= nl) {
                    indice_x = 2*nl - (x - p) -1;
                }
                for (int q = -t; q <= t; q++) {
                    int indice_y = y - q;
                    if ((y - q) < 0) {
                        indice_y = abs(y - q) - 1;
                    }
                    else if ((y - q) >= nc) {
                        indice_y = 2*nc - (y - q) -1;
                    }

                    ecart_type_xy += (entree[indice_x][indice_y] - esperance_xy)*(entree[indice_x][indice_y] - esperance_xy);
                }
            }

            ecart_type_xy /= (2*t + 1)*(2*t + 1);
            ecart_type_xy = sqrt(ecart_type_xy);
            sortie[x*nc + y] = ecart_type_xy;
        }
    }

    qsort(sortie, nc*nl, sizeof(double), doubleComparator);

    return sortie;
}