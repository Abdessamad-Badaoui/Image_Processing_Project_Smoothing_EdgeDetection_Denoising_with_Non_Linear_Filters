#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "pgm.h"



int main(int ac, char **av) {
    clock_t debut, fin ;
    unsigned char **im1 = NULL, **im2 = NULL;
    int nl,nc;
    int sigma[10] = {1,2,3,4,5,6,7,8,9,10};
    FILE *fichier1 = fopen("complexite_filtrage_gaussien.dat","w");
    FILE *fichier2 = fopen("complexite_filtrage_lineaire.dat","w");
    FILE *fichier3 = fopen("complexite_filtrage_separable.dat","w");
    im1 = lectureimagepgm("../imagestp/Archive/formes2g.pgm",&nl,&nc);
    //Filtrage gaussien par FFT
    for (int i=0;i<10;i++){
        debut=clock();
        im2 = filtrage_gauss(NULL, im1,nl,nc,sigma[i]);
        fin=clock();
        fprintf(fichier1,"%i %f\n",sigma[i],((double)fin-debut)/CLOCKS_PER_SEC);
    }

    //Filtrage linéaire gaussien (2D)
    for (int i=0;i<10;i++){
        debut=clock();
        im2 = convolution_masque2D(NULL, im1,nl,nc,sigma[i],4*sigma[i]-1);
        fin=clock();
        fprintf(fichier2,"%i %f\n",sigma[i],((double)fin-debut)/CLOCKS_PER_SEC);
    }

    //Filtrage linéaire séparable 
    for (int i=0;i<10;i++){
        debut=clock();
        im2 = convolution_masque_separable(NULL, im1,nl,nc,sigma[i],4*sigma[i]-1);
        fin=clock();
        fprintf(fichier3,"%i %f\n",sigma[i],((double)fin-debut)/CLOCKS_PER_SEC);
    }

    fclose(fichier1);
    fclose(fichier2);
    fclose(fichier3);

    system("gnuplot -p -e \" set terminal png; set output 'complexite.png'; plot 'complexite_filtrage_gaussien.dat' with linespoints, 'complexite_filtrage_lineaire.dat' with linespoints,'complexite_filtrage_separable.dat' with linespoints  \"");

    


}