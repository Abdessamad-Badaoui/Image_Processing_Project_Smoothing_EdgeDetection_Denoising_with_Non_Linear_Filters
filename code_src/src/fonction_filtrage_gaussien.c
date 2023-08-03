#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

/* Filtrage gaussien par FFT */
unsigned char** filtrage_gauss( unsigned char** sortie,  unsigned char** entree, int nl, int nc, double sigma) 
{ 
    int nb, oldnl,oldnc;
    double ** im1,** im2,**tmp2, ** im3,** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;
    

    im1=imuchar2double(entree,nl,nc);
    oldnl=nl; oldnc=nc;
    im2=padimdforfft(im1,&nl,&nc);

    im3=alloue_image_double(nl,nc); 
    im4=alloue_image_double(nl,nc); 
    im5=alloue_image_double(nl,nc);

    fft(im2,im3,im4,im5,nl,nc);
    tmp2 = alloue_image_double(nl,nc);
    im6=alloue_image_double(nl,nc);
    im7=alloue_image_double(nl,nc);
    
    fftshift(im4,im5,im6,im7,nl,nc);

    double x,y;
    for (int u=0; u<nl; u++) {
        for (int v=0;v<nc; v++) 
        {
            double x  = (((double) u-(nl/2))/nl) * (((double) u-(nl/2))/nl);
            double y  = (((double) v-(nc/2))/nc) * (((double) v-(nc/2))/nc);
            tmp2[u][v] = exp(-2*M_PI*M_PI*sigma*sigma*(x + y)) * im6[u][v];
            im7[u][v] = exp(-2*M_PI*M_PI*sigma*sigma*(x + y)) * im7[u][v];
            
        }
    }
    im8=alloue_image_double(nl,nc);
    im9=alloue_image_double(nl,nc);

    fftshift(tmp2,im7,im8,im9,nl,nc);

    im10 = alloue_image_double(nl,nc);
    ifft(im8,im9,im10,im7,nl,nc);

    unsigned char **im11 = imdouble2uchar(im10,nl,nc);
    
    sortie = crop(im11,0,0,oldnl,oldnc);

    libere_image(im1);
    libere_image(im3);
    libere_image(im4);
    libere_image(im5);
    libere_image(im6);
    libere_image(im7);
    libere_image(im8);
    libere_image(im9);
    libere_image(im10);
    libere_image(im11);
    libere_image(tmp2);

    return sortie;
}


double PSNR(unsigned char ** image1, unsigned char ** image2, int nc, int nl) {
    double sum = 0;
    double sum_inter;
    for (int u=0; u<nl; u++) {
        sum_inter = 0;
        for (int v=0;v<nc; v++) {
            sum_inter += (image1[u][v] - image2[u][v])*(image1[u][v] - image2[u][v]);
        }
        sum += sum_inter;
    }

    double log_arg = ((double)255*255*nc*nl)/sum;
    return 10*log10(log_arg);
}