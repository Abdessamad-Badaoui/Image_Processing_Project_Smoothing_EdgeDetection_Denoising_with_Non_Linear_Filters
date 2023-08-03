#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "pgm.h"

int main(int ac, char **av) {
    
    int nb,nl,nc,oldnl,oldnc;
    unsigned char **im2=NULL,** im1=NULL, **im3=NULL;
    double **im4;
    
    if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
        /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }

    double ***sortie;
    sortie = calloc(nl, sizeof(double **));
    for (int i = 0; i < nl; i++) {
        sortie[i] = calloc(nc, sizeof(double **));
        for (int j = 0; j < nc; j++) {
            sortie[i][j] = calloc(2, sizeof(double));
        }
    }

    double sigma = 5;
    int taille_masque = (int) 4*sigma;

    gradient_lisse_separable(sortie, im1, nl, nc, sigma, taille_masque);

    FILE* data_gradient = fopen("data_gradient_lisse_separable.dat", "w");
    
    for (int i = 0; i < nl; i++) {
        for (int j = 0; j < nc; j++) {
            fprintf(data_gradient, "%i %i %f %f\n", j, nl - i - 1, sortie[i][j][0], sortie[i][j][1]);
        }
    }
    fclose(data_gradient);
    
    double gx;
    double gy;
    double grad_P2;
    double grad_P1;
    int seuil_bas = 0;
    int seuil_haut = 10;

    im4 = alloue_image_double(nl, nc);

    for (int i = 1; i < nl - 1; i++) {
        for (int j = 1; j < nc - 1; j++) {
            double Gx = sortie[i][j][0];
            double Gy = sortie[i][j][1];
            if ((Gx != 0) && (Gy == 0)) {
                gx = Gx / fabs(Gx);
                int a_i = i + (int) gx/fabs(gx);
                int a_j = j;
                int c_i = i - (int) gx/fabs(gx);
                int c_j = j;
                grad_P2 = sqrt(sortie[a_i][a_j][0]*sortie[a_i][a_j][0] + sortie[a_i][a_j][1]* sortie[a_i][a_j][1]);
                grad_P1 = sqrt(sortie[c_i][c_j][0]*sortie[c_i][c_j][0] + sortie[c_i][c_j][1]* sortie[c_i][c_j][1]);
            }
            else if ((Gx == 0) && (Gy != 0)) {
                gy = Gy / fabs(Gy);
                int a_i = i;
                int a_j = j + (int) gy/fabs(gy);
                int c_i = i;
                int c_j = j - (int) gy/fabs(gy);
                grad_P2 = sqrt(sortie[a_i][a_j][0]*sortie[a_i][a_j][0] + sortie[a_i][a_j][1]* sortie[a_i][a_j][1]);
                grad_P1 = sqrt(sortie[c_i][c_j][0]*sortie[c_i][c_j][0] + sortie[c_i][c_j][1]* sortie[c_i][c_j][1]);
            }
            else if ((Gx != 0) && (Gy != 0)) {
                gx = Gx / sqrt(Gx*Gx + Gy*Gy);
                gy = Gy / sqrt(Gx*Gx + Gy*Gy);                
                int a_i = i + (int) gx/fabs(gx);
                int a_j = j + (int) gy/fabs(gy);
                int b_i = i + (int) gx/fmax(fabs(gx), fabs(gy));
                int b_j = j + (int) gy/fmax(fabs(gx), fabs(gy));
                int c_i = i - (int) gx/fabs(gx);
                int c_j = j - (int) gy/fabs(gy);
                int d_i = i - (int) gx/fmax(fabs(gx), fabs(gy));
                int d_j = j - (int) gy/fmax(fabs(gx), fabs(gy));
                grad_P2 = gx/gy*sqrt(sortie[a_i][a_j][0]*sortie[a_i][a_j][0] + sortie[a_i][a_j][1]* sortie[a_i][a_j][1]) + (gy - gx)/gy*sqrt(sortie[b_i][b_j][0]*sortie[b_i][b_j][0] + sortie[b_i][b_j][1]* sortie[b_i][b_j][1]);
                grad_P1 = gx/gy*sqrt(sortie[c_i][c_j][0]*sortie[c_i][c_j][0] + sortie[c_i][c_j][1]* sortie[c_i][c_j][1]) + (gy - gx)/gy*sqrt(sortie[d_i][d_j][0]*sortie[d_i][d_j][0] + sortie[d_i][d_j][1]* sortie[d_i][d_j][1]);
            }
            else {
                continue;
            }

            double grad_P = sqrt(Gx*Gx + Gy*Gy);
            if ((grad_P > grad_P1) && (grad_P > grad_P2)) {
                if (grad_P > seuil_haut) {
                    im4[i][j] = 1;
                    im1[i][j] = (unsigned char) 255;
                }
            }
        }
    }

    for (int i = 1; i < nl - 1; i++) {
        for (int j = 1; j < nc - 1; j++) {
            double Gx = sortie[i][j][0];
            double Gy = sortie[i][j][1];
            if ((Gx != 0) && (Gy == 0)) {
                gx = Gx / fabs(Gx);
                int a_i = i + (int) gx/fabs(gx);
                int a_j = j;
                int c_i = i - (int) gx/fabs(gx);
                int c_j = j;
                grad_P2 = sqrt(sortie[a_i][a_j][0]*sortie[a_i][a_j][0] + sortie[a_i][a_j][1]* sortie[a_i][a_j][1]);
                grad_P1 = sqrt(sortie[c_i][c_j][0]*sortie[c_i][c_j][0] + sortie[c_i][c_j][1]* sortie[c_i][c_j][1]);
            }
            else if ((Gx == 0) && (Gy != 0)) {
                gy = Gy / fabs(Gy);
                int a_i = i;
                int a_j = j + (int) gy/fabs(gy);
                int c_i = i;
                int c_j = j - (int) gy/fabs(gy);
                grad_P2 = sqrt(sortie[a_i][a_j][0]*sortie[a_i][a_j][0] + sortie[a_i][a_j][1]* sortie[a_i][a_j][1]);
                grad_P1 = sqrt(sortie[c_i][c_j][0]*sortie[c_i][c_j][0] + sortie[c_i][c_j][1]* sortie[c_i][c_j][1]);
            }
            else if ((Gx != 0) && (Gy != 0)) {
                gx = Gx / sqrt(Gx*Gx + Gy*Gy);
                gy = Gy / sqrt(Gx*Gx + Gy*Gy);                
                int a_i = i + (int) gx/fabs(gx);
                int a_j = j + (int) gy/fabs(gy);
                int b_i = i + (int) gx/fmax(fabs(gx), fabs(gy));
                int b_j = j + (int) gy/fmax(fabs(gx), fabs(gy));
                int c_i = i - (int) gx/fabs(gx);
                int c_j = j - (int) gy/fabs(gy);
                int d_i = i - (int) gx/fmax(fabs(gx), fabs(gy));
                int d_j = j - (int) gy/fmax(fabs(gx), fabs(gy));
                grad_P2 = gx/gy*sqrt(sortie[a_i][a_j][0]*sortie[a_i][a_j][0] + sortie[a_i][a_j][1]* sortie[a_i][a_j][1]) + (gy - gx)/gy*sqrt(sortie[b_i][b_j][0]*sortie[b_i][b_j][0] + sortie[b_i][b_j][1]* sortie[b_i][b_j][1]);
                grad_P1 = gx/gy*sqrt(sortie[c_i][c_j][0]*sortie[c_i][c_j][0] + sortie[c_i][c_j][1]* sortie[c_i][c_j][1]) + (gy - gx)/gy*sqrt(sortie[d_i][d_j][0]*sortie[d_i][d_j][0] + sortie[d_i][d_j][1]* sortie[d_i][d_j][1]);
            }
            else {
                continue;
            }

            double grad_P = sqrt(Gx*Gx + Gy*Gy);
            if ((grad_P > grad_P1) && (grad_P > grad_P2)) {
                if ((grad_P <= seuil_haut) && (grad_P > seuil_bas)) {
                    if ((im4[i+1][j] == 1) || (im4[i][j+1] == 1) || (im4[i-1][j] == 1) || (im4[i][j-1] == 1)) {
                        im1[i][j] = (unsigned char) 255;
                    }
                }
            }
        }
    }

    ecritureimagepgm(av[2], im1, nl, nc);
    
    libere_image(im4);

    return 0;
}