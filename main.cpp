/////////////////////////////////////
///Projet Option Parisienne
/////////////////////////////////////
//
//Nicolas Klaeylé - Timothée Fabre
// Fichier main

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
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <functional>
#include <chrono> 

// Fichier de code
#include "io.cpp"
#include "utility.cpp"
#include "option.cpp"
#include "wrapper.cpp"
#include "schemes_BS.hpp"
#include "compose.hpp"
#include "mc.hpp"
#include "schemes.hpp"
#include "Density_stopping_time.h"
#include "pricing_mc.hpp"

using namespace std;
using namespace std::chrono;

int main() {

    // Initialisation des inputs

    // Caractéristiques de l'option
    double r(0.1), delta(0.1),sigma(0.1), x(1),T(1),K(100);
    int type_pc(1), type_ud(1), type_io(1);

    // Caractéristiques des barrières
    double d(0.1);
    int n(1),type_1(0),type_2(-1);

    // Caractéristiques barrière 1 et 2
    double b1(0.3),a1(0);
    double b2(0),a2(0);

    // Variable d'affichage
    double p_laplace = -1;

    // Interrupteur et donnée calcul surface de prix
    int spsd(0);
    int spvd(0);
    int nbr(1);


    // Lecture des inputs dans le fichier input.dat
    read_input (r,delta,sigma,x,T,K,type_pc,type_ud,type_io,d,n, \
                type_1, type_2, b1,a1, b2,a2,spsd,spvd,nbr);


    // Initialisation de la classe option - Calcul en Laplace
    string tag1,tag2,tag3,tag4;  // type d'option,type de barrière - 1,type de barrière - 2
    dispatch (type_pc,type_ud,type_io,n,type_1,type_2,tag1,tag2,tag3,tag4);
    option par_lap = option(r,delta,sigma,x,T,K,0,tag1,tag2,tag3,tag4,d,a1,b1,a2,b2);

    // Initialisation simulation MC
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);
    unsigned nb_pas(1000);

    double p(0);
    double h(T/nb_pas);
    double p_type( make_inout(type_io) );

    barrier barrierDown(0,0,"Z"), barrierUp(0,0,"Z");
    make_barrier(barrierUp,barrierDown, a1, b1, a2, b2,n, type_ud, tag3, tag4);

    point_Payoff pPayoff(*payoff_Put);
    if (type_pc == 1)
        {
            pPayoff = *payoff_Call;
        }


    // Initialisation de la structure option parisienne - Calcul en MC
    parisian_Option par_mc(T, r, delta, sigma, K, x, barrierUp, barrierDown, d, p, p_type, pPayoff);

    // Affichage dans la console de la simulation en cours
    cout << "-----------------------------------------------------------" << "\n";
    cout << " " << "\n";
    display_compilator(par_lap);

// -----------------------------------------------------------------------------------------
    // Pricing avec Transformée de Laplace
    cout << "Début du calcul: prix Laplace -----------------------------" << "\n";
    cout << "|  "<< "\n";

    auto start_la = high_resolution_clock::now();
    p_laplace = price_laplace(par_lap,type_ud);
    auto stop_la = high_resolution_clock::now();
    auto duration_la = duration_cast<microseconds>(stop_la - start_la);

    cout << "|  "<<"Prix - Laplace = "<< p_laplace << "\n";
    cout << "|  "<<"Temps de calcul : "<< duration_la.count() << " microseconds" << "\n";

    cout << "|  "<<"\n";
    cout << "Fin du calcul: prix Laplace" << "\n";
    cout << "\n";

// -----------------------------------------------------------------------------------------
    // Pricing Monte Carlo
    cout << "Début du calcul: prix Monte Carlo -------------------------" << "\n";
    cout << "|  "<<"\n";
    unsigned nb_simulations_batch = 1000; // Nombre de simulations par batch
    double precision = 10; // Pr�cision de la simulation

    auto start_mc = high_resolution_clock::now();
    normal_distribution<> G;
    auto X_exact = make_random_scheme(exact_bs_Qtild(sigma, x, h), G, nb_pas);
    auto struct_Parisian_Realisation = make_Parisian_Realisation(X_exact, par_mc, x, 0);
    auto struct_mc = monte_carlo(struct_Parisian_Realisation, gen, nb_simulations_batch, precision);
    par_mc.price = struct_mc.mean(); // Mise � jour de l'objet
    auto stop_mc = high_resolution_clock::now();
    auto duration_mc = duration_cast<microseconds>(stop_mc - start_mc);

    cout << "|  "<<"Prix - Monte Carlo = "<< par_mc.price << "\n";
    cout << "|  "<<"Intervalle de confiance : [" << max(par_mc.price - struct_mc.ic_size(),0.) <<
    ", " << par_mc.price + struct_mc.ic_size() << "]" << endl;
    cout << "|  "<<"Temps de calcul : "<< duration_mc.count() << " microseconds" << "\n";  

    cout << "|  "<<"\n";
    cout <<"Fin du calcul: prix Monte Carlo" << "\n";
    cout << "\n";

// -----------------------------------------------------------------------------------------

    // Si interrupteur = 1, calcul surface de prix

    string name_sd = par_lap.TAG1 + "_K-D.csv"; // nom des fichiers de sortie
    string name_vd = par_lap.TAG1 + "_Vol-D.csv"; 

    string res[ (nbr+1) * (nbr+1) + 1]; // tableaux qui contient les données
    // +1 pour le header

    if (spsd==1) 
        {
            cout << "Début du calcul: surface de prix strike/delai -------------" << "\n";
            cout << "|" << "\n"; 

            //Calcul du pas
            double step_K = (par_lap.K * 1.4 -  par_lap.K * 0.6) / nbr;
            double step_D = (par_lap.d * 1.5 -  par_lap.d * 0.5) / nbr;

            //Calcul du point de départ du grid
            double start_K = par_lap.K * 0.6;
            double start_D = par_lap.d * 0.6; 

            // Initialise le tableau   
            res[0] = "K,D,Prix_Laplace,Prix_MC,Borne_inf_MC,Borne_sup_MC\n"  ;     

            option par_lap_tilde = option(r,delta,sigma,x,T,K,0,tag1,tag2,tag3,tag4,d,a1,b1,a2,b2);
            parisian_Option par_mc_tilde(T, r, delta, sigma, K, x, barrierUp, barrierDown, d, p, p_type, pPayoff);
            int index = 1;

            for (int i= 0. ; i <= nbr ; i++) // grid K
              {
                  for (int j= 0. ; j <= nbr ; j++) // grid D
                    {
                        par_lap_tilde.change_kd(start_K + i * step_K, start_D + j * step_D);  // Laplace

                        par_mc_tilde.change_kd(start_K + i * step_K, start_D + j * step_D);  // MC
                        auto X_exact = make_random_scheme(exact_bs_Qtild(sigma, x, h), G, nb_pas);
                        auto struct_Parisian_Realisation = make_Parisian_Realisation(X_exact, par_mc_tilde, x, 0);
                        auto struct_mc = monte_carlo(struct_Parisian_Realisation, gen, nb_simulations_batch, precision);
                        par_mc_tilde.price = struct_mc.mean();

                        res[index ] = to_string(0.6*K + i * step_K) + "," +
                        to_string(0.5*d + j * step_D) + ","+ to_string(price_laplace(par_lap_tilde,type_ud)) // laplace
                        +","+to_string(par_mc_tilde.price)+","+to_string(max(par_mc_tilde.price - struct_mc.ic_size(),0.))+
                        ","+to_string(par_mc_tilde.price + struct_mc.ic_size())+"\n";
                        index += 1;


                    }
              }
            cout << "| Enregistrement... " << "\n";

            write_csv(res,name_sd,(nbr+1) * (nbr+1) + 1);

            cout << "Fin du calcul: surface de prix strike/delai" << "\n";
            cout << "\n";       
        };

    if (spvd==1) 
        {
            cout << "Début du calcul: surface de prix vol/delai -------" << "\n";  
            cout << "|" << "\n";

            //Calcul du pas
            double step_sigma = (par_lap.sigma * 1.4 -  par_lap.sigma * 0.6) / nbr;
            double step_D = (par_lap.d * 1.5 -  par_lap.d * 0.5) / nbr;

            //Calcul du point de départ du grid
            double start_sigma = par_lap.sigma * 0.6;
            double start_D = par_lap.d * 0.6; 

            // Initialise le tableau   
            res[0] = "sigma,D,Prix_Laplace,Prix_MC,Borne_inf_MC,Borne_sup_MC\n"  ;     
            option par_lap_tilde = option(r,delta,sigma,x,T,K,0,tag1,tag2,tag3,tag4,d,a1,b1,a2,b2);
            parisian_Option par_mc_tilde(T, r, delta, sigma, K, x, barrierUp, barrierDown, d, p, p_type, pPayoff);
            int index = 1;

            for (int i= 0. ; i <= nbr ; i++) // grid K
              {
                  for (int j= 0. ; j <= nbr ; j++) // grid D
                    {
                        par_lap_tilde.change_sd(start_sigma + i * step_sigma, start_D + j * step_D);  // Laplace

                        par_mc_tilde.change_sd(start_sigma + i * step_sigma, start_D + j * step_D);  // MC
                        auto X_exact = make_random_scheme(exact_bs_Qtild(start_sigma + i * step_sigma, x, h), G, nb_pas);
                        auto struct_Parisian_Realisation = make_Parisian_Realisation(X_exact, par_mc_tilde, x, 0);
                        auto struct_mc = monte_carlo(struct_Parisian_Realisation, gen, nb_simulations_batch, precision);
                        par_mc_tilde.price = struct_mc.mean();

                        res[index ] = to_string(0.6*sigma + i * step_sigma) + "," +
                        to_string(0.5*d + j * step_D) + ","+ to_string(price_laplace(par_lap_tilde,type_ud)) // laplace
                        +","+to_string(par_mc_tilde.price)+","+to_string(max(par_mc_tilde.price - struct_mc.ic_size(),0.))+
                        ","+to_string(par_mc_tilde.price + struct_mc.ic_size())+"\n";
                        index += 1;


                    }
              }
            cout << "| Enregistrement... " << "\n";

            write_csv(res,name_vd,(nbr+1) * (nbr+1) + 1);

            cout << "Fin du calcul: surface de prix vol/delai" << "\n";  
            cout << "\n";         
        };


// -----------------------------------------------------------------------------------------

    // Finition
    cout << "Fin des calculs" << "\n";
    cout << "" << "\n";   
    cout << "-----------------------------------------------------------" << "\n";

    return 0;
}
