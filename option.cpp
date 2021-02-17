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

using namespace std;

// Ce fichier contient les fonctions nécessaires pour distinguer les différents types d'options à pricer
//et pour créer les classes correspondantes

// Cette fonction identifie le type d'option 

void dispatch (int type_pc, int type_ud, int type_io,  int n,  \
                int type_1, int type_2, string& tag1, string& tag2, string& tag3, string& tag4) 
{

    // Tag pour le type de fonction
    if (type_pc==1 && type_ud == 1 && type_io == 1 && n==1) {
        tag1 = "SPUIC";
    } else if (type_pc==1 && type_ud == 1 && type_io == -1 && n==1) {
        tag1 = "SPUOC";
    } else if (type_pc==1 && type_ud == -1 && type_io == 1 && n==1) {
        tag1 = "SPDIC";
    } else if (type_pc==1 && type_ud == -1 && type_io == -1 && n==1) {
        tag1 = "SPDOC";
    } else if (type_pc==-1 && type_ud == 1 && type_io == 1 && n==1) {
        tag1 = "SPUIP";
    } else if (type_pc==-1 && type_ud == 1 && type_io == -1 && n==1) {
        tag1 = "SPUOP";
    } else if (type_pc==-1 && type_ud == -1 && type_io == 1 && n==1) {
        tag1 = "SPDIP";
    } else if (type_pc==-1 && type_ud == -1 && type_io == -1 && n==1) {
        tag1 = "SPDOP";

    } else if (type_pc==1 && type_io == 1 && n==2) {
        tag1 = "DPIC";
    } else if (type_pc==1 && type_io == -1 && n==2) {
        tag1 = "DPOC";
    } else if (type_pc==-1 && type_io == 1 && n==2) {
        tag1 = "DPIP";
    } else if (type_pc==-1 && type_io == -1 && n==2) {
        tag1 = "DPOP";

    } else {
        cout << "Error dans input.dat - voir fichier README" << "\n";    
    };

    // Tag pour le nombre de barrière 
    if (n==1) {
        tag2 = "mono barrière";
    } else if (n==2){
        tag2 = "double barrière";
    } else {
        cout << "Error dans input.dat - voir fichier README" << "\n";
    };

    // Tag pour le type de barrière - 1
    if (type_1==0) {
        tag3 = "constante";
    } else if (type_1==1) {
        tag3 = "affine";
    } else if (type_1==2) {
        tag3 = "exponentielle";
    };

    // Tag pour le type de barrière - 2
    if (n==1) {
        tag4 = "-";
    } else {
        if (type_2==0) {
            tag4 = "constante";
        } else if (type_2==1) {
            tag4 = "affine";
        } else if (type_2==2) {
            tag4 = "exponentielle";
        };
    }

    return ;


    
};




// Déclaration de la classe fonction

class option
    {
        public:

            double r, delta,sigma, x,T,K,time;
            string TAG1,TAG2,TAG3,TAG4;
            double d;
            double a1,b1;
            double a2,b2;

            option (double x1, double x2, double x3, double x4, double x5, double x6, double x7,\
            string x8, string x9, string x10, string x22,  \
            double x11, double x12, double x13, double x14, double x15)

                {
                    r = x1;   //taux d'intéret
                    delta = x2;  // taux de div
                    sigma = x3;   // vol
                    x = x4;   // valeur initiale
                    T = x5;   // maturité
                    K = x6;   // strike
                    time = x7;   // moment ou l'on price, égale à zero tout le temps

                    TAG1 = x8;  // type d'option
                    TAG2 = x9;  // type de barrière 
                    TAG3 = x10; // type de barrière - barrière 1
                    TAG4 = x22; // type de barrière - barrière 2

                    d = x11;   // delai 

                    a1 = x12;   // coef affine barrière 1
                    b1 = x13;   // coeff constant barrière 1

                    a2 = x14;   // coef affine barrière 2
                    b2 = x15;   // coeff constant barrière 2

                };

    // Fonction membre 

    double m_()  ;
    double L_(int which_one);
    double k_();
    double b_(int which_one) ; 
    double d_(int which_one) ; 
    void change_kd(double new_strike, double new_delay);
    void change_sd(double new_sigma, double new_delay);

        private:



    };

    // Fonction membre 

    double option::m_()   
        {
            return (r - delta - pow(sigma,2)/2) / sigma;
        };

    double option::L_(int which_one)
        {
            double b;
            if (which_one == 1) 
                {
                    b = b1;
                }
            else if (which_one == 2) 
                {
                    b = b2;
                }

            return b; 
        };

    double option::k_()
        {
            return log(K / x) / sigma;
        };

    double option::b_(int which_one)  
        {
            double res;
            if (which_one == 1) 
                {
                    res = b1;
                }
            else if (which_one == 2) 
                {
                    res = b2;
                }

            return log( res / x ) / sigma ;
        };

    double option::d_(int which_one)
        {
            return ( b_(which_one) - k_() ) / sqrt( d );
        };


    void option::change_kd(double new_strike, double new_delay)
    {
        K = new_strike;
        d = new_delay;
    }

    void option::change_sd(double new_sigma, double new_delay)
    {
        sigma = new_sigma;
        d = new_delay;
    }
































