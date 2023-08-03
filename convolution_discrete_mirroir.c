#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

double** convolution_discrete_mirroir(double **entree, double **filtre, int nl, int nc) {
    
    double **sortie = alloue_image_double(nl, nc);
    
    int x = 0, y = 0;

    for (x = 0; x < nl; x++) {
        for (y = 0; y < nc; y++) {
            double valeur_x_y = 0;
            for (int p = -1; p < 2; p++) {
                double somme_inter = 0;
                int indice_x = x - p;
                if ((x - p) < 0) {
                    indice_x = abs(x - p) -1;
                }
                else if ((x - p) >= nl) {
                    indice_x = 2*nl - (x - p) -1;
                }

                for (int q = -1; q < 2; q++) {
                    int indice_y = y - q;
                    if ((y - q) < 0) {
                        indice_y = abs(y - q) - 1;
                    }
                    else if ((y - q) >= nc) {
                        indice_y = 2*nc - (y - q) -1;
                    }

                    somme_inter += filtre[p+1][q+1]*entree[indice_x][indice_y];
                }
                valeur_x_y += somme_inter;
            }
            sortie[x][y] = valeur_x_y;
        }
    }

    return sortie;
}