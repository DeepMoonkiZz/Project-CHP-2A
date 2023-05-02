CONSIGNES D'UTILISATION DU PROGRAMME

1 - COMPILATION :

Utiliser le Makefile :
  Compiler : make
  Effacer : make clean

--------------------------------------------------------------------------------

2 - INITIALISATION À L'AIDE DU FICHIER parameter.dat :

Entrer les valeurs des constantes et choisir la fonction f. Pour choisir
la fonction f, remplacer la dernière valeur par un entier compris entre 1 et 3.
- Le cas 1 correspond à la solution stationnaire
- Le cas 2 correspond à la deuxième solution stationnaire
- Le cas 3 correspond à la solution instationnaire

--------------------------------------------------------------------------------

3 - CALCUL ET RÉCUPÉRATION DES DONNÉES

a - Cas des solutions stationnaires :

Dans le module mod_display.c, remplacer la ligne 17 comme ceci :

- Pour le cas 1 :

  sprintf(filename, "Solutions/Stationnaire_1/sol_%d.dat", p);

- Pour le cas 2 :

  sprintf(filename, "Solutions/Stationnaire_2/sol_%d.dat", p);

Ensuite, dans main.c, vérifier que la ligne 75 n'est pas en commentaire et que
les lignes 71 et 72 le sont.
Faire ./run pour lancer le programme et récupérer les données.

b - Cas de la solution instationnaire

Dans le module mod_display.c, remplacer la ligne 17 comme ceci :

  sprintf(filename, "Solutions/Instationnaire/sol_%d.dat", p);

Ensuite, dans main.c, vérifier que la ligne 75 est en commentaire et que
les lignes 71 et 72 ne le sont pas.
Faire ./run pour lancer le programme et récupérer les données.

--------------------------------------------------------------------------------

4 - AFFICHAGE DES DONNÉES

a - Cas des solutions stationnaire

Dans le terminal, aller dans le fichier correspondant à la solution
stationnaire à visualiser. Ouvrir gnuplot et taper les commandes suivantes :

> set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)
> set terminal png
> set output "sol_0.png"
> plot "sol_0.dat" u 1:2:3 with image
> set title "Solution stationnaire : cas 1"
> show title

b - Cas de la solution instationnaire

Dans le terminal, aller dans le fichier correspondant à la solution
stationnaire à visualiser. Ouvrir gnuplot et taper les commandes suivantes :

set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)
set terminal png

dt=0.1

do for [i = 1:99] {
	t=i*dt
	set output "sol_".i.".png"
    plot "./sol_".i.".dat" u 1:2:3 with image
    #set title "./sol_".i.".dat"
    set title "t = ".sprintf("%f", t)." s"
    show title
    pause 0.2
}

NB : dt est à changer ainsi que le nombre d'itération en fonction des données
de parameter.dat. Imax = Tmax/dt.

--------------------------------------------------------------------------------

CONTENU DU DOSSIER CODE :
main.c / Makefile

mod_gradient.c / mod_gradient.h : Code du gradient conjugué

mod_display.c / mod_display.h : écriture des données du vecteur u dans un .dat

mod_scheme.c / mod_scheme.h : construction du second membre b dans Ax = b

mod_function.c / mod_function.h : définition des fonctions f, g et h

mod_operations.c / mod_operations.h : définition des différentes opérations
entre vecteurs et matrices, dont la construction de la matrice A non stockée.

Dossiers de données : Instationnaire, Stationnaire_1, Stationnaire_2
