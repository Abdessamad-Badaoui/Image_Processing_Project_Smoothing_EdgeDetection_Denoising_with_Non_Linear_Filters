#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

void gradient_lisse_separable(double*** sortie,  unsigned char** entree, int nl, int nc, double sigma, int taille_masque) 
{
    double ** im1, **lissage, **im2;
    im1=imuchar2double(entree,nl,nc);
    lissage=alloue_image_double(nl,nc);

    double tableau[taille_masque + 1];
    double sum = 0;
    for (int i=0;i<taille_masque +1 ;i++)
    {
        tableau[i] = exp(-i*i/(2*sigma*sigma));
        sum += tableau[i];
    }

    int n = taille_masque;
    for (int x = 0; x < nl; x++) {
        for (int y = 0; y < nc; y++) 
        {
            lissage[x][y] = 0;
            for (int j=-n;j<n+1;j++)
            {
                lissage[x][y] += tableau[abs(j)]*im1[(((j+x)%nl) + nl)%nl][y];
            }
            lissage[x][y] /= (2*sum -1);
        }
    }

    for (int x = 0; x < nl; x++) {
        for (int y = 0; y < nc; y++) {
            sortie[x][y][1] = 0;
            for (int j = -n; j < n+1; j++) {
                sortie[x][y][1] += lissage[x][(((j+y)%nc) + nc)%nc] * (-j/(sigma*sigma)) * tableau[abs(j)];
            }
        }
    }

    for (int x = 0; x < nl; x++) {
        for (int y = 0; y < nc; y++) 
        {
            lissage[x][y] = 0;
            for (int j=-n;j<n+1;j++)
            {
                lissage[x][y] += tableau[abs(j)]*im1[x][(((j+y)%nc) + nc)%nc];
            }
            lissage[x][y] /= (2*sum -1);
        }
    }

    for (int x = 0; x < nl; x++) {
        for (int y = 0; y < nc; y++) {
            sortie[x][y][0] = 0;
            for (int j = -n; j < n+1; j++) {
                sortie[x][y][0] += lissage[(((j+x)%nl) + nl)%nl][y] * (-j/(sigma*sigma)) * tableau[abs(j)];
            }
        }
    }
}