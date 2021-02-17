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
#include "utility.cpp"
#include "hat_func.cpp"


using namespace std;

// Ici on price toutes les options mono barrières

double price_dpoc(option par)
    {
        double m = par.m_();
        double fa = exp( - (par.r + 0.5*m*m) * par.T );

        double res = -1;

        res= laplace(par,0,dpoc_hat) * fa;


        return res ;
    }

double price_dpic(option par)
    {
        double m = par.m_();
        double fa = exp( - (par.r + 0.5*m*m) * par.T );

        double res = -1;

        res= laplace(par,0,dpic_hat) * fa;


        return res ;        
    }

// On price ici les put par parités

double price_dpop(option par)
    {
        double x = par.x;
        double K = par.K;
        double L1 = par.L_(1);
        double L2 = par.L_(2);
        double res = -1;

        option par_tilde = option(par.r,par.delta,par.sigma,1/x,par.T,1/K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1 / L2,par.a2,1 / L1);

        res = price_dpoc(par_tilde) * x * K;

        return res;
    }

double price_dpip(option par)
    {
        double x = par.x;
        double K = par.K;
        double L1 = par.L_(1);
        double L2 = par.L_(2);
        double res = -1;

        option par_tilde = option(par.r,par.delta,par.sigma,1/x,par.T,1/K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1 / L2,par.a2,1 / L1);

        res = price_dpic(par_tilde) * x * K;        
        return res;
    }