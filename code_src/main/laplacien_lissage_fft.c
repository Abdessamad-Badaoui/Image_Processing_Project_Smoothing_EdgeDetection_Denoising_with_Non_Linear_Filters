#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"


main(int ac, char **av) 
{
    int nl, nc;
    int nb, oldnl,oldnc;
    double **im1, ** im2,**tmp2, ** im3,** im4,** im5, ** im6, ** im7, **im8, **im9,**im10;
    unsigned char **im12, **im13;
    //unsigned char **sortie;
    double sigma = 5;

    if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
            
    im12=lectureimagepgm(av[1],&nl,&nc);
    im1=imuchar2double(im12, nl, nc);
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
            tmp2[u][v] = -4 * M_PI*M_PI *(x+y)*exp(-2*M_PI*M_PI*sigma*sigma*(x + y)) * im6[u][v];
            im7[u][v] = -4 * M_PI*M_PI *(x+y)*exp(-2*M_PI*M_PI*sigma*sigma*(x + y)) * im7[u][v];
            
        }
    }
    im8=alloue_image_double(nl,nc);
    im9=alloue_image_double(nl,nc);

    fftshift(tmp2,im7,im8,im9,nl,nc);

    im10 = alloue_image_double(nl,nc);
    ifft(im8,im9,im10,im7,nl,nc);

    //sortie = crop(imdouble2uchar(im10,nl,nc),0,0,oldnl,oldnc);

    for (int i=0; i<oldnl-1; i++){
        for (int j=0; j<oldnc-1; j++){
            //if (im10[i][j] >= seuil){
                if (im10[i][j] * im10[i][j+1] < 0 || im10[i][j] * im10[i+1][j] < 0){
                    im12[i][j] = (unsigned char) 255;
                }
            //}
            //sortie = crop(imdouble2uchar(im10,nl,nc),0,0,oldnl,oldnc);
        }
    }
    
    ecritureimagepgm(av[2], im12, oldnl, oldnc);

    libere_image(im1);
    libere_image(im3);
    libere_image(im4);
    libere_image(im5);
    libere_image(im6);
    libere_image(im7);
    libere_image(im8);
    libere_image(im9);
    libere_image(tmp2);
    libere_image(im10);
    libere_image(im12);
}
    