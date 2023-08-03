# Variables
CC = gcc
CFLAGS = -c -g -std=gnu99
LDFLAGS = -lm

# Dossiers
SRC_DIR = src
OBJ_DIR = obj
MAIN_DIR = main
INC_DIR = include

# Fichiers source
SRCS = $(wildcard $(SRC_DIR)/*.c)

# Fichiers objet
OBJS = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))

# Cibles
EXEC1 = essai1
EXEC2 = essai2
EXEC3 = filtrage_gaussien
EXEC4 = convolution_spatiale
EXEC5 = operateur_premier_ordre
EXEC6 = operateur_deuxieme_ordre_laplacien
EXEC7 = laplacien_lissage_fft
EXEC8 = loi_W_taille_masque
EXEC9 = complexite
EXEC10 = test_conv_discrete_mirroir
EXEC11 = main_filtre_median
EXEC12 = test_filtre_adaptatif_recursif
EXEC13 = test_moyenne_non_locale
EXEC14 = comparaison_methodes_non_lineaire

# Create main directory
$(shell mkdir -p $(MAIN_DIR))

all : $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11) $(EXEC12) $(EXEC13) $(EXEC14)

# Règles de compilation
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -o $@
	
$(OBJ_DIR)/%.o: $(MAIN_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INC_DIR) $< -o $@

$(EXEC1): $(OBJ_DIR)/essai1.o $(OBJ_DIR)/pgm.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC2): $(OBJ_DIR)/essai2.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC3): $(OBJ_DIR)/filtrage_gaussien.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/fonction_filtrage_gaussien.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC4): $(OBJ_DIR)/convolution_spatiale.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/divers.o $(OBJ_DIR)/fonction_filtrage_gaussien.o $(OBJ_DIR)/fonction_convolution_spatiale.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC5): $(OBJ_DIR)/operateur_premier_ordre.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/gradient_lisse_separable.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC6): $(OBJ_DIR)/operateur_deuxieme_ordre.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/laplacien_gaussien.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC7): $(OBJ_DIR)/laplacien_lissage_fft.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o
	$(CC) $^ -o $@  $(LDFLAGS)
	
$(EXEC8): $(OBJ_DIR)/loi_W_taille_masque.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/divers.o $(OBJ_DIR)/fonction_filtrage_gaussien.o $(OBJ_DIR)/fonction_convolution_spatiale.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC9): $(OBJ_DIR)/complexite.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/divers.o $(OBJ_DIR)/fonction_filtrage_gaussien.o $(OBJ_DIR)/fonction_convolution_spatiale.o
	$(CC) $^ -o $@ $(LDFLAGS)	

$(EXEC10): $(OBJ_DIR)/test_conv_discrete_mirroir.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/convolution_discrete_mirroir.o $(OBJ_DIR)/ecart_type_bloc.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC11): $(OBJ_DIR)/main_filtre_median.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/filtre_median_fonction.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC12): $(OBJ_DIR)/test_filtre_adaptatif_recursif.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/filtre_adaptatif_recursif.o $(OBJ_DIR)/divers.o
	$(CC) $^ -o $@ $(LDFLAGS)
	
$(EXEC13): $(OBJ_DIR)/test_moyenne_non_locale.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/non_local_means.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(EXEC14): $(OBJ_DIR)/comparaison_methodes_non_lineaire.o $(OBJ_DIR)/pgm.o $(OBJ_DIR)/fonction_filtrage_gaussien.o $(OBJ_DIR)/fft.o $(OBJ_DIR)/filtre_median_fonction.o $(OBJ_DIR)/filtre_adaptatif_recursif.o $(OBJ_DIR)/non_local_means.o
	$(CC) $^ -o $@ $(LDFLAGS)

# Nettoyage
clean:
	rm -f $(OBJ_DIR)/*.o $(EXEC1) $(EXEC2) $(EXEC3) $(EXEC4) $(EXEC5) $(EXEC6) $(EXEC7) $(EXEC8) $(EXEC9) $(EXEC10) $(EXEC11) $(EXEC12) $(EXEC13) $(EXEC14) *.pgm *.png *.dat *.txt

