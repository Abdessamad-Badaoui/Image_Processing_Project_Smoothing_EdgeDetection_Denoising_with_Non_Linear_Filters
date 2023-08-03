#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

int main(int ac, char **av) {

    int nb, nl, nc, oldnl, oldnc;
    unsigned char **im2=NULL,** im1=NULL, **im3=NULL;
    double **im4;

    if (ac < 4) {
        printf("Usage : %s entree t pourcentile \n",av[0]);
        exit(1);
    }
    
    float p = atof(av[3]);
    int t = atoi(av[2]);

    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1 = lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  {
        puts("Lecture image impossible");
        exit(1);
    }
    
    double **entree = imuchar2double(im1, nl, nc);
    libere_image(im1);

    double values_filtre[3][3] = {{0, -1.0/5, 0}, {-1.0/5, 1, -1.0/5}, {0, -1.0/5, 0}};

    double **filtre = (double **) malloc(3*sizeof(double*));
    assert(filtre != NULL);
    for (int i = 0; i < 3; i++) {
        filtre[i] = (double *) malloc(3*sizeof(double));
        assert(filtre[i] != NULL);
    }
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            filtre[i][j] = values_filtre[i][j];
        }
    }

    double **sortie = convolution_discrete_mirroir(entree, filtre, nl, nc);

    double *ecart_type_vector = ecart_type_bloc(sortie, t, nl, nc);

    p /= 100;
    if (p <= 0.0) {
        printf("%f\n", 1.13*ecart_type_vector[0]);
    } 
    else if ( p >= 1.0) {
        printf("%f\n", 1.13*ecart_type_vector[nl*nc - 1]);
    }
    else {
        unsigned int indice = (unsigned int) floor(p*nc*nl);
        printf("%f\n", 1.13*ecart_type_vector[indice]);
    }

    FILE *sorted_values_file = fopen("hist_values.txt", "w");

    for (int i = 0; i < nl*nc; i++) {
        fprintf(sorted_values_file, "%f\n", ecart_type_vector[i]);
    }

    fclose(sorted_values_file);

    for (int i = 0; i < 3; i++) {
        free(filtre[i]);
    }

    libere_image(entree);
    free(filtre);
    free(ecart_type_vector);

    libere_image(sortie);
    
    return EXIT_SUCCESS;
}