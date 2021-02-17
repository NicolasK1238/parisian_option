/////////////////////////////////////
Projet Option Parisienne - FR
/////////////////////////////////////

Nicolas Klaeylé - Timothée Fabre

README.txt
Ce fichier décrit le fonctionnement du programme, ainsi que la règle de compilation
à utiliser.
----------------------------------------------------------------------------------
Fonctionnement du code: 

 - Règle de Compilation
Le code doit être compilé une seule fois avec la commande 
"g++ -o main main.cpp"
Il est nécessaire d'avoir un compilateur g++ supportant le c++14.

Cette commande construit l'exécutable main qui est appelé par la suite.

- Exécution
Tous les paramètres d'input sont spécifiés dans le fichier input.dat
Après avoir compilé le programme, il n'est plus nécessaire de le recompiler.
Pour chaque nouvelle simulation, il suffit de mettre à jour les paramètres 
d'input dans le fichier input.dat, puis de sauvegarder le fichier, et enfin 
d'exécuter le programme avec la commande "./main"

Cette architecture permet de n'avoir qu'à compiler une seule fois le programme, 
pour ensuite l'exécuter avec des paramètres différents chaque fois que l'on 
veut pricer une option.
Le programme calcule également des surfaces de prix.
Ce processus est assez long, on peut donc l'éteindre depuis le fichier input.dat
Nous allons par la suite détailler le fonctionnement du fichier input.dat, puis 
expliquer le role des différents fichiers du programme

----------------------------------------------------------------------------------
Organisation du fichier input.dat

Nous reproduisons ici le fichier input.dat, ainsi que les différentes contraintes 
associée à l'écriture de la valeur de chacune des variables

        # //Caractéristiques de l'option -----------------------------------
        #
        0.001 //  Taux d'intéret (%) r
        0.0 //  Taux de dividende (%) delta
        0.2 //  Volatilité (%) sigma
        1000 //  Valeur initiale (en cash) x
        3 //  Maturité (en année) T
        1000 //  Strike (en cash) K
        1 //  type_pc (type_pc = 1 pour Call, type_pc = -1 pour Put)
        1 //  type_ud (type_ud = 1 pour Up, type_ud = -1 pour Down)
        1 //  type_io (type_io = 1 pour In, type_io = -1 pour Out)
# //////////////////////////////////////////////////////////////////////////
IMPORTANT : Pour les options à doubles barrières, la variable 
                type_ud n'a plus d'importance (pas besoin de la mettre à jour)
        #
        # //Caractéristiques des barrières ---------------------------------
        #
        1 //  nombre barrière (n = 1 simple barrière, n = 2 double barrière) n
        #
        0.05 // Delai (appelé aussi fenetre dans le rapport), 
        #       commun aux deux barrières (en année) d
        0 //  type barrière 1 (type_1 = 0 barrière cst, 
        #                      type_1 = 1 barrière variable affine, 
        #                      type_1 = 2 barrière variable exponentielle) type_1
        0 //  type barrière 2 (type_2 = 0 barrière cst, 
        #                      type_2 = 1 barrière variable affine, 
        #                      type_2 = 2 barrière variable exponentielle) type_2
# //////////////////////////////////////////////////////////////////////////
IMPORTANT : Si l'on utilise deux barrières, barrière 1 est celle de niveau le plus faible
IMPORTANT : Si l'on utilise seulement une seule barrière, on remplit barrière 1 uniquement
                Seule barrière 1 sera lue, pas besoin de changer les valeurs de barrière 2
        #
        # //Caractéristiques barrière 1 ------------------------------------
        #
        800// niveau barrière 1 constant (en cash / absolu) b1
        0 //  Coeff a barrière affine/exponentielle a1 (absolu)
        #
        # //Caractéristiques barrière 2 ------------------------------------
        #
        1200// niveau barrière 2 constant (en cash / absolu) b2
        0 //  Coeff a barrière affine/exponentielle a2 (absolu)
        # 
        # //Surface de prix ------------------------------------
# //////////////////////////////////////////////////////////////////////////
IMPORTANT : Ce calcul prend du temps
IMPORTANT : Le grid est centré sur les valeurs de strike/délai et vol/délai
                initialisé plus haut
        0  // spsd Surface de prix strike/délai 
        #   (spsd = 1 si l'on souhaite calculer la surface de prix, 0 sinon)
        0  // spvd Surface de prix volatilité/délai 
        #   (spvd = 1 si l'on souhaite calculer la surface de prix, 0 sinon)
        10 // taille du grid sur lequel on calcule la surface de prix
# //////////////////////////////////////////////////////////////////////////
IMPORTANT : Les données des surfaces de prix sont enregistrées dans des 
                fichiers csv généré automatiquement
----------------------------------------------------------------------------------
Rôle des différents fichiers du programme

On peut distinguer 5 catégories de fichier:

 - fichier de donnée d'input: input.dat

 - fichier de travail sur donnée: io.cpp
    Ce fichier contient toute les fonctions permettant de lire et d'écrire dans 
    des fichiers de donnée

 - fichier pour le pricing par transformée de Laplace:
    error_func.cpp : ce fichier permet de calculer la fonction d'erreur comlexe à
                     partir d'une décomposition en série entière de la fonction de
                     Faddeeva. Il est adapté de la librairie PNL (référence dans 
                     le rapport)

    hat_func.cpp : ce fichier calcule toutes les tranformées de Laplace du prix de 
                   chacune des options. Il contient également la fonction 
                   permettant d'inverser la transformée de Laplace

    utility.cpp : ce fichier calcule des fonctions basiques utilisées à de 
                  nombreuses reprises

    option.cpp : ce fichier déclare la classe option et les fonctions membres 
                 associées

    price_single.cpp : ce fichier contient les fonctions permettant de pricer
                       les options simples barrières

    special_esp.cpp : ce fichier calcule les transformées de Laplace des temps d'
                      arrêt nécessaires pour pricer les options doubles barrières

    price_double.cpp : ce fichier contient les fonctions permettant de pricer
                       les options doubles barrières

    wrapper.cpp : ce fichier combine les fonctions de price_single.cpp et 
                  price_double.cpp en une seule fonction. Il applique également
                  des tests aux données d'input pour éviter les cas dégénérés
                  et pour détecter cetaines erreurs de frappes dans le fichier 
                  input.dat

 - fichier pour le pricing par Méthode de Monte Carlo:
    compose.hpp : ce fichier permet d'effectuer des compositions de fonctions 
    définies dans les structures 

    schemes.hpp : ce fichier définit les structures de schémas de discrétisation 
    qui sont utilisées pour effectuer les itérations dans la procédure de simulation

    schemes_BS.hpp : ce fichier est composé des structures de schéma de 
    discrétisation exacte BS dans les deux mesures de probabilite Q et Qtilde

    Density_stopping_time.hpp : ce fichier permet de simuler les temps aléatoires 
    de dernier passage dans la méthode de Monte Carlo

    mc.hpp : ce fichier définit les structures utilisées dans la méthode de 
    Monte Carlo, notamment le calcul d'estimateurs et d'intervalles de confiance

    pricing_mc.hpp : ce fichier est composé de structures de définition de l'option, 
    des structures de barrières time-dependent et permet de lancer la procédure de
    simulation d'une réalisation du payoff

 - fichier main.cpp
 ----------------------------------------------------------------------------------
