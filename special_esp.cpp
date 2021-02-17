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
#include <complex>  
#include <functional>
// Fichier de code
#include "option.cpp"
#include "error_func.cpp"
#include "utility.cpp"

using namespace std;

// Ce fichier calcule toute les espérances nécessaire à la valuation des options double barrières

std::complex<double> inter_minus(option par, int which_one,std::complex<double> lambda )
    //  E( exp(- lambda * Tb-) )
    {
        std::complex<double> theta = theta_(lambda);
        double b = par.b_(which_one);
        double D = par.d;

        std::complex<double> res;

        if (b < 0)
            {
                res = exp(theta * b) / psi(theta * sqrt(D));
            }
        else
            {
                res = exp(-lambda*D) * (1 - 2 * normale_cdf( -b/sqrt(D) ) ) +
                ( exp(-theta * b) * normale_cdf( theta * sqrt(D) - b/sqrt(D) ) +
                exp(theta * b) * normale_cdf( -theta * sqrt(D) - b/sqrt(D) ) ) 
                /  psi(theta * sqrt(D));
            }
        return res;
    }

std::complex<double> inter_plus(option par, int which_one,std::complex<double> lambda )
    //  E( exp(- lambda * Tb+) )
    {
        std::complex<double> theta = theta_(lambda);
        double b = par.b_(which_one);
        double D = par.d;

        std::complex<double> res;

        if (b < 0)
            {
                res = exp(-theta * b) / psi(theta * sqrt(D));
            }
        else
            {
                res = exp(-lambda*D) * (1 - 2 * normale_cdf( b/sqrt(D) ) ) +
                ( exp(theta * b) * normale_cdf( theta * sqrt(D) + b/sqrt(D) ) +
                exp(-theta * b) * normale_cdf( -theta * sqrt(D) + b/sqrt(D) ) ) 
                /  psi(theta * sqrt(D));
            }
        return res;
    }

std::complex<double> second_plus(option par, int which_one,std::complex<double> lambda )
    //  E( exp(+ theta * ZTb-) )
    {
        std::complex<double> theta = theta_(lambda);
        double b = par.b_(which_one);
        double D = par.d;

        std::complex<double> res;

        theta = - theta;

        if (b < 0)
            {
                res = exp(-theta * b) * psi( theta*sqrt(D) );
            }
        else
            {
                res = 2. * psi(theta*sqrt(D)) * normale_cdf(-b/sqrt(D)) * exp(-theta * b) +
                exp(lambda*D) * ( normale_cdf((b/sqrt(D))+theta*sqrt(D))-exp(-2.*theta*b) *    
                normale_cdf( -(b/sqrt(D))+theta*sqrt(D))  );
            }
        return res;
    }

std::complex<double> second_minus(option par, int which_one,std::complex<double> lambda )
    //  E( exp(- theta * ZTb+) )
    {
        std::complex<double> theta = theta_(lambda);
        double b = par.b_(which_one);
        double D = par.d;

        std::complex<double> res;

        if (b < 0)
            {
                res = exp(-theta * b) * psi( -theta*sqrt(D) );
            }
        else
            {
                res = 2. *   psi(-theta*sqrt(D)) * normale_cdf(b/sqrt(D)) * exp(-theta * b) +
                exp(lambda*D) * ( normale_cdf((-b/sqrt(D))-theta*sqrt(D))-exp(-2.*theta*b) *    
                normale_cdf( (b/sqrt(D))-theta*sqrt(D))  );
            }
        return res;
    }

std::complex<double> first_minus(option par, std::complex<double> lambda )
// E[ exp(-lambda*Tb1-) * 1(Tb1- < Tb2+)  ]
    {
        std::complex<double> theta = theta_(lambda);
        double D = par.d;
        double b1 = par.b_(1);
        double b2 = par.b_(2);      


        std::complex<double> a1 = exp(theta*b1) * second_minus(par,2,lambda) / psi(theta*sqrt(D));
        std::complex<double> a2 = exp(-theta*b2) * second_plus(par,1,lambda) / psi(theta*sqrt(D));

        std::complex<double> res = ( inter_minus(par,1,lambda) 
        - a1 * inter_plus(par,2,lambda)) / (1. - a1*a2);

        return res;
    }

std::complex<double> first_plus(option par, std::complex<double> lambda )
// E[ exp(-lambda*Tb2+) * 1(Tb2+ < Tb1-)  ]
    {
        std::complex<double> theta = theta_(lambda);
        double D = par.d;
        double b1 = par.b_(1);
        double b2 = par.b_(2);      


        std::complex<double> a1 = exp(theta*b1) * second_minus(par,2,lambda) / psi(theta*sqrt(D));
        std::complex<double> a2 = exp(-theta*b2) * second_plus(par,1,lambda) / psi(theta*sqrt(D));

        std::complex<double> res = ( inter_plus(par,2,lambda) 
        - a2 * inter_minus(par,1,lambda)) / (1. - a1*a2);

        return res;
    }


    