# Compilateur utilisé
CC=mpicc

# Options en mode optimisé - La variable DEBUG est définie comme fausse
OPTIM_FLAG = -g -Wall -lm
# Options en mode debug - La variable est DEBUG est définie comme vraie
DEBUG_FLAG = -g -DDEBUG -lm

# On choisit comment on compile
CFLAGS = $(OPTIM_FLAG)

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = main.c mod_gradient.c mod_operations.c mod_scheme.c mod_function.c /home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_display.c /home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/parameter_and_display/mod_parameter.c /home/segal/Documents/MatMeca/S8/CHP/Project-CHP-2A/fonction/charge.c

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CFLAGS) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
