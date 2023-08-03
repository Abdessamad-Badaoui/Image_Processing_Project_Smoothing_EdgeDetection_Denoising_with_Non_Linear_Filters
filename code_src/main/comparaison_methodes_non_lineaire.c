#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "pgm.h"

static unsigned char** extract(unsigned char **im, int demie_taille_filtre, int nl, int nc){
    unsigned char** result = (unsigned char**) malloc(sizeof(unsigned char*) * (nl-2*demie_taille_filtre));
    for (int i = 0; i < nl-2*demie_taille_filtre; i++) {
        result[i] = (unsigned char*) malloc(sizeof(unsigned char) * (nc-2*demie_taille_filtre));
        for (int j = 0; j < nc-2*demie_taille_filtre; j++) {
            result[i][j] = im[demie_taille_filtre + i][demie_taille_filtre + j];
        }
    }
    return result;
}

int main(int ac, char **av) {
    unsigned char **image_reference;
    int nl, nc;

    image_reference = lectureimagepgm("../imagestp/Archive/formes2g.pgm", &nl, &nc);

    FILE *methodes_filtrage_non_lineaires = fopen("methodes_non_lineares.txt", "w");
    fputs("filtre_median\n", methodes_filtrage_non_lineaires);
    fputs("filtre_adaptatif_recursif\n", methodes_filtrage_non_lineaires);
    fputs("moyenne_non_locale\n", methodes_filtrage_non_lineaires);
    fclose(methodes_filtrage_non_lineaires);

    FILE *bruit_gaussien_file = fopen("resultat_psnr_gauss_30.txt", "w");
    FILE *bruit_poisson_file = fopen("resultat_psnr_poiss_3.txt", "w");
    FILE *bruit_poivre_sel_file = fopen("resultat_psnr_poivre_sel_15.txt", "w");
    FILE *bruit_speckle_file = fopen("resultat_psnr_speckle_35.txt", "w");

    FILE *bruit_gauss_exection_time = fopen("resultat_perf_gauss_30.txt", "w");
    FILE *bruit_poiss_exection_time = fopen("resultat_perf_poiss_3.txt", "w");
    FILE *bruit_poivre_sel_exection_time = fopen("resultat_perf_poivre_sel_15.txt", "w");
    FILE *bruit_speckle_exection_time = fopen("resultat_perf_speckle_35.txt", "w");

    FILE **fichier_psnr = (FILE **) malloc(4 * sizeof(FILE *));
    fichier_psnr[0] = bruit_gaussien_file;
    fichier_psnr[1] = bruit_poisson_file;
    fichier_psnr[2] = bruit_poivre_sel_file;
    fichier_psnr[3] = bruit_speckle_file;

    FILE **fichier_perf = (FILE **) malloc(4 * sizeof(FILE *));
    fichier_perf[0] = bruit_gauss_exection_time;
    fichier_perf[1] = bruit_poiss_exection_time;
    fichier_perf[2] = bruit_poivre_sel_exection_time;
    fichier_perf[3] = bruit_speckle_exection_time;

    char *bruit_gaussien_file_name = "../imagestp/Archive/formes2bb30.pgm";
    char *bruit_poisson_file_name = "../imagestp/Archive/formes2p3.pgm";
    char *bruit_poivre_sel_file_name = "../imagestp/Archive/formes2sp15.pgm";
    char *bruit_speckle_file_name = "../imagestp/Archive/formes2s35.pgm";

    char **fichiers_nom = (char **) malloc(4 * sizeof(char *));
    fichiers_nom[0] = bruit_gaussien_file_name;
    fichiers_nom[1] = bruit_poisson_file_name;
    fichiers_nom[2] = bruit_poivre_sel_file_name;
    fichiers_nom[3] = bruit_speckle_file_name;


    int demi_taille_filtre_median = 5;
    clock_t debut, fin;
    double duree;
    double psnr_value;
    unsigned char **image_reference_cropped = extract(image_reference, demi_taille_filtre_median, nl, nc);
    
    for (int i = 0; i < 4; i++) {
        char *bruit_file_name = fichiers_nom[i];
        FILE *psnr_file = fichier_psnr[i];
        FILE *per_file = fichier_perf[i];

        unsigned char **image_res_filtre_median;
        unsigned char **image_res_filtre_adaptatif_recursif;
        unsigned char **image_res_non_local_means;
        unsigned char **entree = lectureimagepgm(bruit_file_name, &nl, &nc);
        debut = clock();
        image_res_filtre_median = filtre_median(entree, nl, nc, demi_taille_filtre_median);
        fin = clock();
        duree = ((double) (fin - debut)) / CLOCKS_PER_SEC;

        unsigned char **image_res_filtre_median_cropped = extract(image_res_filtre_median, demi_taille_filtre_median, nl, nc);
        libere_image(image_res_filtre_median);

        psnr_value = PSNR(image_reference_cropped, image_res_filtre_median_cropped, nc - 2*demi_taille_filtre_median, nl - 2*demi_taille_filtre_median);

        fprintf(psnr_file, "%f\n", psnr_value);
        fprintf(per_file, "%f\n", duree);

        debut = clock();
        image_res_filtre_adaptatif_recursif = filtre_adaptatif_recursif(entree, nl, nc, 30, 150);
        fin = clock();
        duree = ((double) (fin - debut)) / CLOCKS_PER_SEC;

        psnr_value = PSNR(image_reference, image_res_filtre_adaptatif_recursif, nc, nl);

        fprintf(psnr_file, "%f\n", psnr_value);
        fprintf(per_file, "%f\n", duree);

        debut = clock();
        image_res_non_local_means = nl_means_patch(entree, nl, nc, 21, 5, 0.4*30, 30);
        fin = clock();
        duree = ((double) (fin - debut)) / CLOCKS_PER_SEC;

        psnr_value = PSNR(image_reference, image_res_non_local_means, nc, nl);

        fprintf(psnr_file, "%f\n", psnr_value);
        fprintf(per_file, "%f\n", duree);       
    }
    
    for (int i = 0; i < 4; i++) {
        fclose(fichier_psnr[i]);
        fclose(fichier_perf[i]);
    }

    return EXIT_SUCCESS;
}