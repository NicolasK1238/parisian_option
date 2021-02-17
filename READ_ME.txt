/////////////////////////////////////
///Projet Option Parisienne
/////////////////////////////////////
//
//Nicolas Klaeylé - Timothée Fabre

----------------------------------------------------------------------------------
Fonctionement du code: 

 - Compilation
Le code doit etre compilé une seule fois avec la commande 
"g++ -o main main.cpp"

Cette commande construit l'executable main qui est appelé par la suite.

- Exécution
Tout les paramètres d'input sont spécifiés dans le fichier input.dat
Après avoir compilé le programme, il n'est plus nécessaire de le recompiler.
Pour chaque nouvelle simulation, il suffit de mettre à jour les paramètres 
d'input dans le fichier input.dat, puis de sauvegarder le fichier, et enfin 
d'exécuter le programme avec la commande "./main"

Cette architecture permet de n'avoir qu'à compiler une seule fois le programme, 
pour ensuite l'exécuter avec des paramètres différents chaque fois que l'on 
veut pricer une option.

 - NB Le programme calcule également des surfaces de prix.
Ce processus est assez long, on peut donc l'éteindre depuis le fichier input.dat

----------------------------------------------------------------------------------
Organisation du fichier input.dat

Le fichier d'input contient, dans l'ordre les variables suivantes:

