#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"



main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL, **im3=NULL, **im4 = NULL, **im5 = NULL;
  
  if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

  FILE *data2D = fopen("data2D.dat","w");
  FILE *data_separable = fopen("data_separable.dat","w");

  int sigma = 5;
  double eqm2D, eqm_sep;

  for (int i = 1; i<6*sigma; i++)
  {
    im2 = filtrage_gauss(NULL,im1,nl,nc,sigma);
    im4 = convolution_masque2D(NULL,im1,nl,nc,sigma,i);
    im5 = convolution_masque_separable(NULL,im1,nl,nc,sigma,i);
    eqm2D = eqm(im2,im4,nl,nc);
    eqm_sep = eqm(im2,im5,nl,nc);
    fprintf(data2D, "%i %f\n",i,eqm2D);
    fprintf(data_separable, "%i %f\n",i,eqm_sep);
  }

  fclose(data2D);

  fclose(data_separable);

  system("gnuplot -p -e \" set terminal png; set output 'res2D.png'; plot 'data2D.dat' with linespoints title 'My Plot'\"");
  
  system("gnuplot -p -e \" set terminal png; set output 'res_separable.png'; plot 'data_separable.dat' with linespoints title 'My Plot'\"");

  libere_image(im1);
  libere_image(im2);
  libere_image(im4);
  libere_image(im5);
}