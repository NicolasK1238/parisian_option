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
#include <numeric>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <complex> 
#include <cfloat>
#include <cmath>
#include <limits>
// Fichier de code
#include "option.cpp"
#include "error_func.cpp"

using namespace std;


// Ce fichier contient l'implémentation de fct permettant de réaliser de petits calculs simples
// Cela permettra dans d'autre bloc de rendre le code plus clair
// Dans la mesure du possible, on gardera les meme notation que dans l'article et le rapport


double normale_cdf(double u)  // n'existe pas directement dans c++
    {
        return erfc(-u / sqrt(2))/2; // on la calcule  avec la fonction erfc
    };

std::complex<double> normale_cdf(std::complex<double> z)  
    {      
        return erfc(-z / sqrt(2) )/2.;
    };


double psi(double x)
    {
        double PI  =3.141592653589793238463;
        return 1 + x * sqrt(2*PI) * exp( pow(x,2)/2 ) * normale_cdf(x);
    };

std::complex<double>  psi(std::complex<double> z)
    {
        double PI  =3.141592653589793238463;
        return 1. + z * sqrt(2*PI) * exp( pow(z,2)/2. ) * normale_cdf(z);
    };

std::complex<double> theta_(std::complex<double> z)  
    {
        return sqrt(2. * z);
    };

int fact(int n) 
    {
    if (n == 0 || n == 1)
        return 1;
    else
        return n * fact(n - 1);
    };

int combinaison(int p, int q)
    {
        return fact(q) / ( fact(p) * fact(q-p) );
    };

// Cette fonction vérifie que le calcul en Laplace est licite
// La transformée de laplace d'un PDIC n'est défini que pour L <= x
// La transformée de laplace d'un PUIC n'est défini que pour x <= L
// De part le fonctionement du programme, cette limitation se transmet à l'évalutation de toute les fct

void existence_check(option par, int type_ud)
{
    std::string option_tag = par.TAG2;

    if (option_tag.compare("mono barrière")==0)
    {
        if (type_ud == 1 && par.x > par.b1) 
            {
                cout << "Le pricing en Laplace n'est pas défini pour cette combinaison (barrière - valeure initiale)"<<"\n";
                cout << "Voir fichier READ_ME" << "\n";
            } 
        else if (type_ud == -1 && par.x < par.b1)        
            {
                cout << "Le pricing en Laplace n'est pas défini pour cette combinaison (barrière - valeure initiale)"<<"\n";
                cout << "Voir fichier READ_ME" << "\n";
            }
    }
    if (option_tag.compare("double barrière")==0)
    {
        if (par.x < par.b1 || par.x > par.b2) 
            {
                cout << "Le pricing en Laplace n'est pas défini pour cette combinaison ";
                cout << "(barrière basse - valeure initiale - barrrière haute)"<<"\n";
                cout << "Voir fichier READ_ME" << "\n";
            } 
    }
    return;

}



// Les fcts ci dessous servent uniquement à harmonisé les notations entre le code MC et le code Laplace

char make_inout(int type_io)
{
    if (type_io == 1)
    {
        return 'I';
    }
    else
    {
        return 'O';
    }
    
}

