#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"



int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL, **im3=NULL, **im4 = NULL, **im5 = NULL;
  
  if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

  FILE *data2D = fopen("data2D.dat","w");
  FILE *data_separable = fopen("data_separable.dat","w");
  FILE *data_loi_W = fopen("data_loi_W.dat","w");

  int sigma;
  double eqm2D, eqm_sep;
  int n_min;
  double eqm_min = INFINITY;

  // Calcul de la taille du masque en fonction de sigma
  for (sigma=1; sigma<11; sigma++) {
    eqm_min = INFINITY;
    for (int i = 1; i<6*sigma; i++)
    {
      im2 = filtrage_gauss(NULL,im1,nl,nc,sigma);
      im5 = convolution_masque_separable(NULL,im1,nl,nc,sigma,i);
      eqm_sep = eqm(im2,im5,nl,nc);
      if (eqm_sep<eqm_min) {
        eqm_min = eqm_sep;
        n_min = i;
      }
    }
    libere_image(im2);
    libere_image(im5);
    fprintf(data_loi_W, "%i %i\n",sigma,n_min);
  }

  fclose(data_loi_W);

  system("gnuplot -p -e \" set terminal png; set output 'res_loi_W.png'; plot 'data_loi_W.dat' with linespoints title 'My Plot'\"");

  libere_image(im1);
  
  return EXIT_SUCCESS;
}