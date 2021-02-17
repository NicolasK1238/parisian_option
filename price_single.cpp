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

double price_spdic(option par,int which_one)
    {
        double b = par.b_(which_one);
        double m = par.m_();
        double fa = exp( - (par.r + 0.5*m*m) * par.T );
        double L = par.L_(which_one);

        double res = -1;

        if (b <= 0.) 
            {
                res= laplace(par,which_one,pdic_hat_inf) * fa;
            } 
        else
            {
                res = laplace(par,which_one,pdic_hat_sup) * fa;
            } 

        return res;
    };

double price_spdoc(option par,int which_one)
    {
        double b = par.b_(which_one);
        double m = par.m_();
        double fa = exp( - (par.r + 0.5*m*m) * par.T );
        double L = par.L_(which_one);

        double res = -1;


        if (b <= 0.) 
            {   
                res = fa *  laplace(par,which_one,pdoc_hat_inf) ;
            } 
        else
            {
                res = fa * laplace(par,which_one,pdoc_hat_sup);
            } 
        return res;
    };

double price_spuic(option par,int which_one)
    {
        double b = par.b_(which_one);
        double m = par.m_();
        double fa = exp( - (par.r + 0.5*m*m) * par.T );
        double L = par.L_(which_one);

        double res = -1;

        if (b >= 0.) 
            {
                res = laplace(par,which_one,puic_hat_sup) * fa;
            } 
        else
            {
                res = laplace(par,which_one,puic_hat_inf) * fa;
            } 
        return res;
    };




double price_spuoc(option par,int which_one)
    {
        double b = par.b_(which_one);
        double m = par.m_();
        double fa = exp( - (par.r + 0.5*m*m) * par.T );
        double L = par.L_(which_one);

        double res = -1;

        if (b >= 0.) 
            {
                res =  laplace(par,which_one,puoc_hat_sup) * fa;
            } 
        else
            {
                res = laplace(par,which_one,puoc_hat_inf) * fa ;
            } 

        return res;
    };


// On price ici les put par parités

double price_spdop(option par,int which_one)
    {
        double x = par.x;
        double K = par.K;
        double L = par.L_(which_one);
        double res = -1;

        option par_tilde = option(par.r,par.delta,par.sigma,1/x,par.T,1/par.K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1 / L,par.a2,1 / L);

        res = price_spuoc(par_tilde,which_one) * x * K;

        return res;
    };

double price_spuop(option par,int which_one)
    {
        double x = par.x;
        double K = par.K;
        double L = par.L_(which_one);
        double res = -1;

        option par_tilde = option(par.r,par.delta,par.sigma,1/x,par.T,1/par.K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1 / L,par.a2,1 / L);

        res = price_spdoc(par_tilde,which_one) * x * K;
        return res;
    };

double price_spuip(option par,int which_one)
    {
        double x = par.x;
        double K = par.K;
        double L = par.L_(which_one);
        double res = -1;

        option par_tilde = option(par.r,par.delta,par.sigma,1/x,par.T,1/par.K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1 / L,par.a2,1 / L);

        res = price_spdic(par_tilde,which_one) * x * K;
        return res ;

    };


double price_spdip(option par,int which_one)
    {
        double x = par.x;
        double K = par.K;
        double L = par.L_(which_one);
        double res = -1;

        option par_tilde = option(par.r,par.delta,par.sigma,1/x,par.T,1/par.K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1 / L,par.a2,1 / L);
        
        res = price_spuic(par_tilde,which_one) * x * K;
        return res ;
    };
