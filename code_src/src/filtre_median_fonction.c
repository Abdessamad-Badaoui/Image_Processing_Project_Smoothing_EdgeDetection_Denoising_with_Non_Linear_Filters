#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

unsigned char ** filtre_median(unsigned char **entree, int nl, int nc, int demi_taille_filtre) {
    int n = demi_taille_filtre;

    unsigned int *histogram = calloc(256, sizeof(int));

    unsigned char **sortie = alloue_image(nl, nc);

    int seuil = 2*n*n + 2*n;

    for (int i = n; i < nl-n; i++) {
        int j = n;
        unsigned int index_hist = 0;
        for (int u = -n; u <= n; u++) {
            for (int v = -n; v <= n; v++) {
                index_hist = (unsigned int) entree[i+u][j+v];
                histogram[index_hist]++;
            }
        }
        
        int somme = 0;
        for (int l = 0; l < 256; l++) {
            somme += histogram[l];
            if (somme >= seuil) {
                sortie[i][j] = (unsigned char) l;
                break;
            }
        }

        for (j = n+1; j < nc-n; j++) {
            for (int u = -n; u <= n; u++) {
                index_hist = (unsigned int) entree[i+u][j-n-1];
                histogram[index_hist]--;
                index_hist = (unsigned int) entree[i+u][j+n];
                histogram[index_hist]++;
            }

            int somme = 0;
            for (int l = 0; l < 256; l++) {
                somme += histogram[l];
                if (somme >= seuil) {
                    sortie[i][j] = (unsigned char) l;
                    break;
                }
            }
        }

        for (int l = 0; l < 256; l++) {
            histogram[l] = 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nc - n; j++) {
            sortie[i][j] = entree[i][j];
        }
    }
    for (int j = nc - n; j < nc; j++) {
        for (int i = 0; i < nl - n; i++) {
            sortie[i][j] = entree[i][j];
        }
    }
    for (int i = nl - n; i < nl; i++) {
        for (int j = n; j < nc; j++) {
            sortie[i][j] = entree[i][j];
        }
    }
    for (int j = 0; j < n; j++) {
        for (int i = n; i < nl; i++) {
            sortie[i][j] = entree[i][j];
        }
    }

    free(histogram);
    
    return sortie;
}