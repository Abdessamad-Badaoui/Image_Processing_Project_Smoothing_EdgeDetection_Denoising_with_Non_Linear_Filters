#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

int main(int ac, char **av) {
    
    int nb,nl,nc, oldnl,oldnc;
    unsigned char **im2=NULL,** im1=NULL, **im3=NULL, **im4 = NULL, **im5 = NULL;
    
    if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
        /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

    double **sortie;
    sortie = calloc(nl, sizeof(double *));
    for (int i = 0; i < nl; i++) {
        sortie[i] = calloc(nc, sizeof(double));
    }

    double sigma = 5;
    int taille_masque = (int) 4*sigma;

    laplacien_gaussien(sortie, im1, nl, nc, sigma, taille_masque);

    FILE *laplacien_file = fopen("laplacien_gaussien_formes2.txt", "w");

    for (int i = 0; i < nl; i++) {
        for (int j = 0; j < nl ; j++) {
            fprintf(laplacien_file, "%f\t", sortie[i][j]);
        }
        fprintf(laplacien_file, "\n");
    }
    fclose(laplacien_file);
    
    double signe_i_j;
    double signe_i_j_plus_1;
    double signe_i_plus_1_j;
    for (int i = 0; i < nl - 1; i++) {
        for (int j = 0; j < nc - 1 ; j++) {
            if (sortie[i][j] < 1e-18) {
                signe_i_j = 0;
            }
            else {
                signe_i_j = sortie[i][j]/fabs(sortie[i][j]);
            }
            if (sortie[i][j+1] < 1e-18) {
                signe_i_j_plus_1 = 0;
            }
            else {
                signe_i_j_plus_1 = sortie[i][j+1]/fabs(sortie[i][j+1]);
            }
            if (sortie[i+1][j] < 1e-18) {
                signe_i_plus_1_j = 0;
            }
            else {
                signe_i_plus_1_j = sortie[i+1][j]/fabs(sortie[i+1][j]);
            }

            if (signe_i_j != signe_i_j_plus_1) {
                im1[i][j+1] = (unsigned char) 255;
            }
            if (signe_i_j != signe_i_plus_1_j) {
                im1[i+1][j] = (unsigned char) 255;
            }        
        }
    }

    ecritureimagepgm(av[2], im1, nl, nc);    

    return 0;
}