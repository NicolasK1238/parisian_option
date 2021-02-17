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
#include<fstream>
#include <iomanip>
#include <stdlib.h>

// Fichier de code
#include "option.cpp"

using namespace std;


// Ce fichier contient toutes les fonctions pour lire/écrire dans un fichier de donnée
// On mettra egalement ici les fonctions affichant des résultats dans la console


// Cette fonction permet de lire le fichier d'input

void read_input (double& r, double& delta, double& sigma, double& x, double& T, double& K, \
                int& type_pc, int& type_ud, int& type_io, double& d, int& n,  \
                int& type_1, int& type_2, double& b1, double& a1, double& b2,  \
                double& a2, int& spsd, int& spvd, int& nbr) 
{
    double val[23];
    int i=0;

    ifstream inputFile("input.dat");  // ouvre le fichier
    std::string str;

    if (inputFile.is_open())  // récupère les données
    {  
        while ( std::getline(inputFile, str) )
        {
            if (str[0] == '#') continue;
            val[i] = atof(str.c_str());
            i += 1;
        }
    }
    // affecte les données

    r = val[0];
    delta = val[1];
    sigma = val[2];
    x = val[3];
    T = val[4];
    K = val[5];
    type_pc = val[6];
    type_ud = val[7];
    type_io = val[8];
    n = val[9];
    d = val[10];
    type_1 = val[11];
    type_2 = val[12];
    b1 = val[13];
    a1 = val[14];
    b2 = val[15];
    a2 = val[16];
    spsd = val[17];
    spvd = val[18];
    nbr = val[19];

    return ;

}


// Cette fonction affiche les information sur la simulation en cours

void display_compilator(option par)
    {
        cout << "Option Parisienne " << "\n";   
        cout << "\n";
        cout << "Caractéristique de l'option -------------------------------" << "\n";
        cout << par.TAG1 << "\n";
        cout << "Strike = " << par.K << " | Maturité = " << par.T << " | S0 = " << par.x << "\n";
        cout << "Taux = " << par.r << " | Vol = " << par.sigma << " | Div = " << par.delta << "\n";  
        cout << "\n";
        cout << "Caractéristiques des barrières ----------------------------" << "\n";
        cout << "\n";
        cout << "Option "<< par.TAG2 << "\n";
        cout << "Delai = " << par.d <<  "\n";
        cout << "\n";

        string tag1 = "constante";
        string tag2 = "affine";
        string tag3 = "exponentielle";
        string tag4 =  "double barrière";


        cout << "Barrière 1 : barrière "<< par.TAG3<< "\n";
        if (tag1.compare(par.TAG3) == 0) 
            {
                cout << "Barrière 1 = "<< par.b1<< "\n";
            }
        else if (tag2.compare(par.TAG3)== 0) 
            {
                cout << "Barrière 1 = "<< par.a1 << " * t + " << par.b1<< "\n";
            }
        else
            {
                cout << "Barrière 1 = exp("<< par.a1 << " * t) * " << par.b1<< "\n";
            }
        
        cout << "\n";

        if (tag4.compare(par.TAG2)== 0) 
            {
                cout << "Barrière 2 : barrière "<< par.TAG4<< "\n";
                if (tag1.compare(par.TAG4)== 0) 
                    {
                        cout << "Barrière 2 = "<< par.b2<< "\n";
                    }
                else if (tag2.compare(par.TAG4)== 0) 
                    {
                        cout << "Barrière 2 = "<< par.a2 << " * t + " << par.b2<< "\n";
                    }
                else
                    {
                        cout << "Barrière 2 = exp("<< par.a2 << " * t) * " << par.b2<< "\n";
                    }
            }
    cout << "\n";
        return;
        

    }



void write_csv(string res[],string name, int n)
{
    std::ofstream myfile;
    myfile.open (name);


    for (int i = 0; i < n; i++) 
        { 
            myfile << res[i] ; 
        }    

    myfile.close();

    return ;

}