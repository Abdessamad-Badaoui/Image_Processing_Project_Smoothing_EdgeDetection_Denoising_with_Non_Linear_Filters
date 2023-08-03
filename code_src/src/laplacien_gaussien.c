#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

void laplacien_gaussien( double** sortie,  unsigned char** entree, int nl, int nc, double sigma, int taille_masque) 
{
    double ** im1, **im2;
    im1=imuchar2double(entree,nl,nc);
    
    double tableau[taille_masque + 1];
    for (int i=0;i<taille_masque +1 ;i++)
    {
        tableau[i] = exp(-i*i/(2*sigma*sigma));
    }

    double **temp;
    
    temp = sortie;

    int n = taille_masque;
    for (int x=0; x<nl; x++) {
        for (int y=0;y<nc; y++) {
            temp[x][y] = 0;
            for (int i =-n;i<n+1;i++)
            {
                for (int j=-n;j<n+1;j++)
                {
                    temp[x][y] += (tableau[abs(i)]*tableau[abs(j)]*im1[(((i+x)%nl)+nl)%nl][(((j+y)%nc)+nc)%nc]*((i*i + j*j)/(2*sigma*sigma) - 1)/(M_PI*sigma*sigma*sigma*sigma) );
                }
            }
        }
    }
    
    libere_image(im1);
}