#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

static int mirroir(int indice, int n) {
    int indice_mirroir = indice;
    if (indice < 0) {
        indice_mirroir = abs(indice) -1;
    }
    else if (indice >= n) {
        indice_mirroir = 2*n - indice -1;
    }
    return indice_mirroir;
}

double difference_patch(double **entree, int nl, int nc, int p_x, int p_y, int q_x, int q_y, int demi_taille_patch) {
    int r = demi_taille_patch;
    double result = 0.0;
    int q_x_i_mirroir;
    int q_y_j_mirroir;
    int p_x_i_mirroir;
    int p_y_j_mirroir;
    for (int i = -r; i <= r; i++) {
        for (int j = -r; j <= r; j++) {
            q_x_i_mirroir = mirroir(q_x + i, nl);
            q_y_j_mirroir = mirroir(q_y + j, nc);
            p_x_i_mirroir = mirroir(p_x + i, nl);
            p_y_j_mirroir = mirroir(p_y + j, nc);
            result += (entree[q_x_i_mirroir][q_y_j_mirroir] - entree[p_x_i_mirroir][p_y_j_mirroir]) * (entree[q_x_i_mirroir][q_y_j_mirroir] - entree[p_x_i_mirroir][p_y_j_mirroir]);
        }
    }
    result = result/(2*r + 1)/(2*r + 1);
    return result;
}



unsigned char **nl_means_patch(unsigned char **entree, int nl, int nc, int demi_taille_region, int demi_taille_patch, double h, double sigma) {
    double **entree_double = imuchar2double(entree, nl, nc);
    int t = demi_taille_region;

    double **sortie_double = alloue_image_double(nl, nc);

    int q_x;
    int q_y;
    double difference_patch_value;
    double *poids = (double *) malloc((2*t + 1) * (2*t + 1) * sizeof(double));
    double *im_q = (double *) malloc((2*t + 1) * (2*t + 1) * sizeof(double));
    for (int p_x = 0; p_x < nl; p_x++) {
        for (int p_y = 0; p_y < nc; p_y++) {
            double w_p_q = 0;
            for (int i = -t; i <= t; i++) {
                for (int j = -t; j <= t; j++) {
                    q_x = mirroir(p_x + i, nl);
                    q_y = mirroir(p_y + j, nc);
                    difference_patch_value = difference_patch(entree_double, nl, nc, p_x, p_y, q_x, q_y, demi_taille_patch);
                    w_p_q = exp(-fmax(difference_patch_value - 2*sigma*sigma, 0)/h/h);
                    poids[(i + t)*(2*t + 1) + j + t] = w_p_q;
                    im_q[(i + t)*(2*t + 1) + j + t] = entree_double[q_x][q_y];
                }
            }
            double ims_p = 0;
            double facteur_normalisation = 0;
            for (int i = 0; i < (2*t + 1) * (2*t + 1); i++) {
                ims_p += poids[i] * im_q[i];
                facteur_normalisation += poids[i];
            }
            ims_p /= facteur_normalisation;

            sortie_double[p_x][p_y] = ims_p;
        }
    }

    free(poids);
    free(im_q);
    unsigned char **sortie = imdouble2uchar(sortie_double, nl, nc);
    libere_image(entree_double);
    libere_image(sortie_double);
    
    return sortie;
}
