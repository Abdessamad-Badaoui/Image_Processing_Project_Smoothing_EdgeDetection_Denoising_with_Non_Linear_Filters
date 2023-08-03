#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

double **filtre_adaptatif_recursif(double **entree, int nl, int nc, double k, int number_iterations) {

    double **sortie = alloue_image_double(nl, nc);
    double **image_intermediaire = alloue_image_double(nl, nc);

    for (int i = 0; i < nl; i++) {
        for (int j = 0; j < nc; j++) {
            image_intermediaire[i][j] = entree[i][j];
        }
    }

    double **alternate[2] = {image_intermediaire, sortie};

    int index = 0;

    int iteration = 0;

    while (iteration < number_iterations) {
        double **source = alternate[index];
        double **destination = alternate[(index + 1)%2];

        for (int x = 0; x < nl; x++) {
            for (int y = 0; y < nc; y++) {
                double valeur_x_y = 0;
                double normalisation = 0;
                for (int p = -1; p < 2; p++) {
                    double somme_inter = 0;
                    double normalisation_inter = 0;
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
                        
                        int indice_x_minus_1 = (indice_x - 1) < 0 ? abs(indice_x - 1) - 1 : indice_x - 1 ;
                        int indice_x_plus_1 = (indice_x + 1) >= nl ? 2*nl - (indice_x + 1) - 1 : indice_x + 1;
                        int indice_y_minus_1 = (indice_y - 1) < 0 ? abs(indice_y - 1) - 1 : indice_y - 1;
                        int indice_y_plus_1 = (indice_y + 1) >= nc ? 2*nc - (indice_y + 1) - 1 : indice_y + 1;

                        double w_t = exp(-(((source[indice_x_plus_1][indice_y] - source[indice_x_minus_1][indice_y])*(source[indice_x_plus_1][indice_y] - source[indice_x_minus_1][indice_y]) + (source[indice_x][indice_y_plus_1] - source[indice_x][indice_y_minus_1])*(source[indice_x][indice_y_plus_1] - source[indice_x][indice_y_minus_1]))/2/k/k));
                        
                        normalisation_inter += w_t;
                        somme_inter += w_t*source[indice_x][indice_y];
                    }
                    normalisation += normalisation_inter;
                    valeur_x_y += somme_inter;
                }
                destination[x][y] = valeur_x_y/normalisation;
            }
        }
        iteration++;
        index = (index + 1)%2;
    }
    sortie = alternate[index];
    image_intermediaire = alternate[(index + 1)%2];
    libere_image(image_intermediaire);
    return sortie;
}