#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

int main(int ac, char **av) {
    
    int nl, nc;
    unsigned char ** im1=NULL;

    if (ac < 4) {
        printf("Usage : %s entree sortie demi_taille_filtre\n",av[0]);
        exit(1);
    }

    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  {
        puts("Lecture image impossible");
        exit(1);
    }

    unsigned char **sortie = filtre_median(im1, nl, nc, atoi(av[3]));

    ecritureimagepgm(av[2],sortie,nl,nc);

    libere_image(im1);
    libere_image(sortie);
    
    return EXIT_SUCCESS;
}