#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

int main(int ac, char **av) {
    int nl, nc;
    unsigned char ** im1=NULL;
    double sigma = 45;
    int r = 7;
    int t = 35;
    double h = 0.35*sigma;

    if (ac < 3) {
        printf("Usage : %s entree sortie \n",av[0]);
        exit(1);
    }

    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  {
        puts("Lecture image impossible");
        exit(1);
    }

    unsigned char **sortie = nl_means_patch(im1, nl, nc, t, r, h, sigma);

    ecritureimagepgm(av[2],sortie,nl,nc);

    libere_image(im1);
    libere_image(sortie);

    return EXIT_SUCCESS;
}