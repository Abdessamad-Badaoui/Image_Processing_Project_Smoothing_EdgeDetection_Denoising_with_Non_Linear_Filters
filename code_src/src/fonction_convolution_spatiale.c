#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"


unsigned char** convolution_masque2D( unsigned char** sortie,  unsigned char** entree, int nl, int nc, double sigma, int taille_masque) 
{
    double ** im1, **im2;
    im1=imuchar2double(entree,nl,nc);
    im2=alloue_image_double(nl,nc);
    
    double tableau[taille_masque + 1];
    double sum = 0;
    for (int i=0;i<taille_masque +1 ;i++)
    {
        tableau[i] = exp(-i*i/(2*sigma*sigma));
        sum += tableau[i];
    }

    int n = taille_masque;
    for (int x=0; x<nl; x++) {
        for (int y=0;y<nc; y++) 
        {
            im2[x][y] = 0;
            for (int i =-n;i<n+1;i++)
            {
                for (int j=-n;j<n+1;j++)
                {
                    im2[x][y] += ( tableau[abs(i)]*tableau[abs(j)]*im1[(((i+x)%nl)+nl)%nl][(((j+y)%nc)+nc)%nc] );
                }
            }

            im2[x][y] = im2[x][y] / ((2*sum -1)*(2*sum -1));
        }
    }

    sortie = imdouble2uchar(im2,nl,nc);
    libere_image(im1);
    libere_image(im2);
    return sortie;
}


unsigned char** convolution_masque_separable( unsigned char** sortie,  unsigned char** entree, int nl, int nc, double sigma, int taille_masque) 
{
    double ** im1, **g1, **im2;
    im1=imuchar2double(entree,nl,nc);
    g1=alloue_image_double(nl,nc);
    im2 = alloue_image_double(nl,nc);


    double tableau[taille_masque + 1];
    double sum = 0;
    for (int i=0;i<taille_masque +1 ;i++)
    {
        tableau[i] = exp(-i*i/(2*sigma*sigma));
        sum += tableau[i];
    }

    int n = taille_masque;
    for (int x=0; x<nl; x++) {
        for (int y=0;y<nc; y++) 
        {
            for (int j=-n;j<n+1;j++)
            {
                g1[x][y] += (tableau[abs(j)]*im1[x][(((j+y)%nc)+nc)%nc])/(2*sum - 1);
            }
        }
    }

    for (int x=0; x<nl; x++) {
        for (int y=0;y<nc; y++) 
        {
            for (int j=-n;j<n+1;j++)
            {
                im2[x][y] += (tableau[abs(j)]*g1[(((j+x)%nl)+nl)%nl][y])/(2*sum - 1);
            }
        }
    }

    sortie = imdouble2uchar(im2,nl,nc);
    libere_image(im1);
    libere_image(g1);
    libere_image(im2);
    return sortie;

}
