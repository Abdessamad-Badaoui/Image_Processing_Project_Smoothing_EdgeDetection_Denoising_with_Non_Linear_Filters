#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"



main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nb,nl,nc, oldnl,oldnc;
  unsigned char **im2=NULL,** im1=NULL, **im9=NULL, **im4;
  double** im5, ** im6, ** im7, **im8,**im10;
  
  if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
	/* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

  FILE *data = fopen("data.dat","w");
  im9 = lectureimagepgm("../imagestp/Archive/formes2g.pgm",&nl,&nc);
  double i = 0;
  double arg_max_PSNR = 0;
  double max_PSNR = 0; 
  while(i<=10)
  {
    im2 = filtrage_gauss(NULL, im1,nl,nc,i);
    double psnr = PSNR(im2,im9,nl,nc);
    fprintf(data, "%f %f\n",i,psnr);
    if (psnr > max_PSNR)
    {
        max_PSNR = psnr;
        arg_max_PSNR = i;
    }
    i = i + 0.1;
  }
  fclose(data);
  
  system("gnuplot -p -e \" set terminal png; set output 'res.png'; plot 'data.dat' with linespoints title 'My Plot'\"");

  im4 = filtrage_gauss(NULL,im1,nl,nc,arg_max_PSNR);

  ecritureimagepgm(av[2],im4,nl,nc);
  
  libere_image(im4);
  libere_image(im1);
  libere_image(im9);
  libere_image(im2);
}