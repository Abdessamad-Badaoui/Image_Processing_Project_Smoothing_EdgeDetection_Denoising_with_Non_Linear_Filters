#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "pgm.h"

int main(int ac, char **av) {
    int nl, nc;
    unsigned char ** im1=NULL;

    if (ac < 3) {
        printf("Usage : %s entree sortie \n",av[0]);
        exit(1);
    }

    /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
    im1=lectureimagepgm(av[1],&nl,&nc);
    if (im1==NULL)  {
        puts("Lecture image impossible");
        exit(1);
    }

    double **entree = imuchar2double(im1, nl, nc);
    libere_image(im1);

    int k;
    int number_iterations;
    
    unsigned char **reference = lectureimagepgm("../imagestp/Archive/formes2g.pgm", &nl, &nc);
    double **reference_double = imuchar2double(reference, nl, nc);
    libere_image(reference);

    char *data_file_name = malloc(256);
    assert(data_file_name != NULL);
    
    int longeur_nom = 0;
    for (int i = strlen(av[1])- 4 - 1; i >= 0; i--) {
        if (av[1][i] == '/') {
            break;
        }
        longeur_nom++;
    }

    char *sourcefile_name_no_ext = malloc(longeur_nom + 1);
    assert(sourcefile_name_no_ext != NULL);

    int l;
    for (l = strlen(av[1]) - 4 - longeur_nom; l < strlen(av[1]) - 4; l++) {
        sourcefile_name_no_ext[l - (strlen(av[1]) - 4 - longeur_nom)] = av[1][l];
    }
    sourcefile_name_no_ext[l - (strlen(av[1]) - 4 - longeur_nom)] = '\0';

    sprintf(data_file_name, "data_%s.dat", sourcefile_name_no_ext);

    FILE *data_psnr_result = fopen(data_file_name, "w");

    number_iterations = 100;
    for (k = 6; k < 10; k+=2) {
        double **sortie = filtre_adaptatif_recursif(entree, nl, nc, k, number_iterations);
        double psnr_score = psnr_double(sortie, reference_double, nl, nc);
        libere_image(sortie);
        fprintf(data_psnr_result, "%i %f\n", k, psnr_score);
    }
    fclose(data_psnr_result);
    char *plot_command = malloc(256);
    assert(plot_command != NULL);
    sprintf(plot_command, "gnuplot -p -e \" set terminal png; set output 'res_psnr_k_%s.png'; plot '%s' with linespoints title 'Filtre adaptatif : PSNR en fonction de k %s'\"", sourcefile_name_no_ext, data_file_name, sourcefile_name_no_ext);

    system(plot_command);

    FILE *data_performace = fopen("k_fct_d_ecart_type.dat", "w");
    FILE *bruit_gaussien_files = fopen("bruit_gaussien_files.log", "r");
    FILE *valeurs_bruit_gaussien = fopen("valeurs_bruit_gaussien.log", "r");

    char *buffer_file = NULL;
    char *buffer_value = NULL;
    size_t buffer_file_size = 0;
    size_t *buffer_value_size = 0;
    int value;

    while (getline(&buffer_file, &buffer_file_size, bruit_gaussien_files) != -1) {
        buffer_file[strcspn(buffer_file, "\n")] = '\0';
        printf("%s\n", buffer_file);
        if (getline(&buffer_value, &buffer_value_size, valeurs_bruit_gaussien) != -1) {
            value = atoi(buffer_value);
            char **bruit_gaussien_image = lectureimagepgm(buffer_file, &nl, &nc);
            if (bruit_gaussien_image==NULL)  {
                puts("Lecture image impossible");
                exit(1);
            }
            double **bruit_gaussien_image_double = imuchar2double(bruit_gaussien_image, nl, nc);
            libere_image(bruit_gaussien_image);
            int best_k = 0;
            double best_psnr_score = 0;
            for (k = 2; k < 100; k+=2) {
                double **sortie = filtre_adaptatif_recursif(bruit_gaussien_image_double, nl, nc, k, number_iterations);
                double psnr_score = psnr_double(sortie, reference_double, nl, nc);
                libere_image(sortie);
                if (psnr_score > best_psnr_score) {
                    best_psnr_score = psnr_score;
                    best_k = k;
                }
            }
            fprintf(data_performace, "%i %i\n", best_k, value);
        }
    }

    system("gnuplot -p -e \" set terminal png; set output 'meilleur_k_en_fct_ecart_type.png'; plot 'k_fct_d_ecart_type.dat' with linespoints title 'Meilleur k en fonction de l'ecart-type'\"");

    fclose(data_performace);
    fclose(bruit_gaussien_files);
    fclose(valeurs_bruit_gaussien);


    free(plot_command);
    free(sourcefile_name_no_ext);
    free(data_file_name);

    double **sortie = filtre_adaptatif_recursif(entree, nl, nc, 200, 100);

    unsigned char **sortie_uchar = imdouble2uchar(sortie, nl, nc);
    ecritureimagepgm(av[2],sortie_uchar,nl,nc);

    libere_image(sortie);
    libere_image(entree);
    libere_image(sortie_uchar);
    libere_image(reference_double);
    
    return EXIT_SUCCESS;
}