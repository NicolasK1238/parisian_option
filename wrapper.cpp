/////////////////////////////////////
///Projet Option Parisienne
/////////////////////////////////////
//
//Nicolas Klaeylé - Timothée Fabre

#pragma once
// Module
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <functional>
#include <cmath>
#include <numeric>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// Fichier de code
#include "option.cpp"
#include "price_single.cpp"
#include "price_double.cpp"


using namespace std;

// Ce fichier contient les registres pour pricer chacune des 12 types d'options


double price_laplace(option par, int type_ud)

    {
        std::string option_tag = par.TAG1;
        std::string type1 = par.TAG3;
        std::string type2 = par.TAG4;
        std::string ref_1 = "constante";
        std::string ref_2 = "-";       

        // Le code Laplace ne price que les options a barrières constantes

        if (type1.compare(ref_1) != 0 ) 
        {   
            cout << "Le code Laplace ne price que les options à barrières constantes"<<"\n";
            cout << "Voir la première barrière"<<"\n";
            return 0;
        }

        if (type2.compare(ref_1) != 0 && type2.compare(ref_2) != 0) 
        {   
            cout << "Le code Laplace ne price que les options à barrières constantes"<<"\n";
            cout << "Voir la deuxième barrière"<<"\n";
            return 0;
        }

        // Cas dégénéré
        if ( par.TAG2.compare("double barrière") == 0 && par.b1 == par.b2)
        {
            cout << "Cas dégénéré: les deux barrières ne peuvent pas etre égale" <<  "\n";
            return 0;
        }

        // On cherche à présent à appeler une subroutine de pricing pour chacun des types d'options
        // On appelle les fonction tag par tag

        std::string ref1 = "SPUIC";
        std::string ref2 = "SPUOC";
        std::string ref3 = "SPDIC";
        std::string ref4 = "SPDOC";

        std::string ref5 = "SPUIP";
        std::string ref6 = "SPUOP";
        std::string ref7 = "SPDIP";
        std::string ref8 = "SPDOP";

        std::string ref9 = "DPOC";
        std::string ref10 = "DPIC";

        std::string ref11 = "DPOP";
        std::string ref12 = "DPIP";

        double res = -1;

        if (option_tag.compare(ref1) == 0 ) // SPUIC
            {              
                res = price_spuic(par,1);
            }
        else if (option_tag.compare(ref2) == 0 ) // SPUOC
            {
                res = price_spuoc(par,1);
            }
        else if (option_tag.compare(ref3) == 0 ) // SPDIC
            {
                res = price_spdic(par,1);
            }
        else if (option_tag.compare(ref4) == 0 ) // SPDOC
            {
                res = price_spdoc(par,1);
            }
//
        else if (option_tag.compare(ref5) == 0 ) // SPUIP
            {            
                res = price_spuip(par,1);
            }
        else if (option_tag.compare(ref6) == 0 ) // SPUOP
            {
                res = price_spuop(par,1);
            }
        else if (option_tag.compare(ref7) == 0 ) // SPDIP
            {
                res = price_spdip(par,1);
            }
        else if (option_tag.compare(ref8) == 0 ) // SPDOP
            {
                res = price_spdop(par,1);
            }
//
        else if (option_tag.compare(ref9) == 0 ) // DPOC
            {            
                res = price_dpoc(par);
            }
        else if (option_tag.compare(ref10) == 0 ) // DPIC
            {
                res = price_dpic(par);
            }
        else if (option_tag.compare(ref11) == 0 ) // DPOP
            {
                res = price_dpop(par);
            }
        else if (option_tag.compare(ref12) == 0 ) // DPIP
            {
                res = price_dpip(par);
            }
        else 
            {
                cout << "Erreur: le type de l'option est mal spécifié" << "\n";
            }

        return res;    

    };




