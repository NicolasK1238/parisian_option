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
#include "utility.cpp"
#include "special_esp.cpp"

using namespace std;

// Ce fichier contient l'implémentation des fonctions transformé de laplace. 
// C'est egalement ici que l'on calcule la transformée de laplace inverse

// BS simple ///////////////////////////////////////////////////////////////////////////////////////////////////
std::complex<double> sc_hat(option par, std::complex<double> lambda, int which_one)
    {  // which one n'a pas de role ici, c'est juste pour maintenir un format cst sur tt les fonction de ce fichier
        // Calcule des variables internes
        std::complex<double> theta = theta_(lambda);
        double K = par.K;
        double m =par.m_();
        double k = par.k_();
        double sigma = par.sigma;
        double x = par.x;

        // Variable de retour
        std::complex<double> res;

        // Variable intermediaire
        std::complex<double> term1;
        std::complex<double> term2;
        std::complex<double> term3;   

        if ( K >= x) 
            {
                term1 = K * exp( (m - theta) * k ) / theta;
                term2 = 1./(m-theta);
                term3 = 1./(m+sigma-theta);

                res = term1 * (term2 - term3);

            }
        else 
            {
                term1 = (2*K) / ( m * m - theta * theta );
                term2 = (2*x) / ( (m+sigma) * (m+sigma) - theta * theta );
                term3 = ( K * exp( (m + theta) * k ) / theta ) * ( (1./(m+theta)) - (1./(m+theta+sigma)) );
                
                res = term1 - term2 + term3;
            };

        return res;
    }

// PDIC - PDOC ////////////////////////////////////////////////////////////////////////////////////////////////////

std::complex<double> pdic_hat_inf(option par, std::complex<double> lambda,int which_one)
    {
        // Calcule des variables internes
        std::complex<double> theta = theta_(lambda);
        double K = par.K;
        double m = par.m_();
        double k = par.k_();
        double sigma = par.sigma;
        double b = par.b_(which_one);
        double D = par.d;        
        double L = par.L_(which_one);
        double d = par.d_(which_one);

        std::complex<double> plus = psi(+ theta * sqrt(D) );
        std::complex<double> minus = psi(- theta * sqrt(D) );      
        double PI  =3.141592653589793238463;  

        // Variable de retour
        std::complex<double> res = 0;

        // Variable intermediaire
        std::complex<double> term1;
        std::complex<double> term2;
        std::complex<double> term3;  
        std::complex<double> term4;
        std::complex<double> term5;  
        std::complex<double> term6;  

        if(real(lambda) <= (m+sigma)*(m+sigma)*0.5)
        {
            cout << "Erreur: Re(Lambda) trop faible" << "\n";
        }


        if ( K > L) 
            {
                term1 = (minus/plus) * (K/theta) * exp( 2*b*theta + (m-theta)*k );
                term2 = 1./(m-theta);
                term3 = 1./(m+sigma-theta);

                res = term1 * (term2 -term3);
            }
        else 
            {
                term1 = exp( b*(m + theta) ) / plus;
                term2 = (2*K) / (m*m - theta*theta);
                term3 = psi(m * sqrt(D)) - m*sqrt(2*PI*D)*exp(D*m*m*0.5)*normale_cdf(m*sqrt(D) + d);
                term4 = (2*L) / ( (m+sigma)*(m+sigma) - theta*theta );
                term5 = psi( (m+sigma) * sqrt(D)) - (m+sigma)*sqrt(2*PI*D)*
                exp(D*(m+sigma)*(m+sigma)*0.5)*normale_cdf( (m+sigma)*sqrt(D) + d );

                res = term1 * (term2 * term3 - term4 * term5 );

                term1 = K * exp( k*(m + theta) ) / (plus * theta);
                term2 = ( 1./(m+theta) ) - ( 1./(m+theta+sigma) );
                term3 = plus - theta * sqrt(2*PI*D) * exp(lambda*D) * normale_cdf(theta * sqrt(D) - d);
                term4 = K * sqrt(2*PI*D) * exp(2*b*theta + lambda*D + (m-theta)*k ) 
                * normale_cdf(-theta * sqrt(D) - d) / plus;
                term5 = (1. /(m+sigma-theta)) - (1. /(m-theta));

                res = res + term1 * term2 * term3 + term4*term5;

            };   
        return res;    
    }


        
std::complex<double> pdoc_hat_inf(option par, std::complex<double> lambda,int which_one)
    {
        return sc_hat(par,lambda,which_one) - pdic_hat_inf(par,lambda,which_one);
    }

std::complex<double> pdoc_hat_sup(option par, std::complex<double> lambda,int which_one)
    {
 
        double m = par.m_();
        double b = par.b_(which_one);
        double D = par.d;        
        double L = par.L_(which_one);
        std::complex<double> theta = theta_(lambda);
        std::complex<double> res = 0;

        double term1 = exp(m*b) * L ;
        std::complex<double> term2 = exp( - theta*b) * normale_cdf(theta * sqrt(D) - b / sqrt(D) ) + 
        exp( + theta*b) * normale_cdf(-theta * sqrt(D) - b / sqrt(D) );
       
        option par_tilde = option(par.r,par.delta,par.sigma,1,par.T,par.K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1,par.a2,1);

        res = pdoc_hat_inf(par_tilde,lambda,which_one)*term1*term2;     

        return res * term1 * term2;
    }

std::complex<double> pdic_hat_sup(option par, std::complex<double> lambda,int which_one)
    {
        return sc_hat(par,lambda,which_one) - pdoc_hat_sup(par,lambda,which_one);
    }


// PUIC - PUOC ////////////////////////////////////////////////////////////////////////////////////////////////////

std::complex<double> puic_hat_sup(option par, std::complex<double> lambda,int which_one)
    {
        // Calcule des variables internes
        std::complex<double> theta = theta_(lambda);
        double K = par.K;
        double m = par.m_();
        double k = par.k_();
        double sigma = par.sigma;
        double b = par.b_(which_one);
        double D = par.d;        
        double L = par.L_(which_one);
        double d = par.d_(which_one);

        std::complex<double> plus = psi(+ theta * sqrt(D) );
        std::complex<double> minus = psi(- theta * sqrt(D) );      
        double PI  =3.141592653589793238463;  

        // Variable de retour
        std::complex<double> res = 0;

        // Variable intermediaire
        std::complex<double> term1;
        std::complex<double> term2;
        std::complex<double> term3;  
        std::complex<double> term4;
        std::complex<double> term5;  

        if(real(lambda) <= (m+sigma)*(m+sigma)*0.5)
        {
            cout << "Erreur: Re(Lambda) trop faible" << "\n";
        }

        if ( K > L) 
            {

                term1 = 2*sqrt(2*PI*D)*exp( (m-theta)*b )/plus;
                term2 = ( (K*m)/(m*m-theta*theta) ) * exp(D * m * m * 0.5) * normale_cdf(m * sqrt(D) + d);
                term3 = ( L*(m+sigma)/( (m+sigma)*(m+sigma) - theta*theta ) ) * exp( D * 0.5 *(m+sigma)*(m+sigma) ) *
                normale_cdf( (m+sigma) * sqrt(D) + d );
                term4 = K * sqrt(2*PI*D) * (1./plus) * exp(-2*b*theta + lambda*D + (m+theta)*k ) * 
                normale_cdf(d - theta * sqrt(D));                
                term5 = (1./(m+sigma+theta))  - (1./(m+theta)) ;

                res = term1 * (term2 - term3) + term4 * term5;

                term1 = (K/theta) * (1./plus) * exp( (m-theta)*k );
                term2 = (1./(m-theta))  - (1./(m+sigma-theta)) ;
                term3 = plus - theta * sqrt(2*PI*D) * exp(lambda*D) * normale_cdf(d + theta * sqrt(D));

                res = res + term1 * term2 * term3;
            } 
        else
            {
                term1 = 2. * exp( (m-theta)*b ) / plus;
                term2 = K * psi( m * sqrt(D) ) / (m*m - theta*theta);
                term3 = L * psi( (m+sigma) * sqrt(D) ) / ( (m+sigma)*(m+sigma) - theta*theta );
                term4 = (minus/plus) * (K/theta) * exp( -2.*b*theta + (m+theta)*k );
                term5 = (1./(m+theta)) - (1./(m+theta+sigma));

                res = term1 * (term2 - term3) + term4 * term5;

            }

        return res;
    };

std::complex<double> puoc_hat_sup(option par, std::complex<double> lambda,int which_one)
    {
        return sc_hat(par,lambda,which_one) - puic_hat_sup(par,lambda,which_one);
    }

std::complex<double> puoc_hat_inf(option par, std::complex<double> lambda,int which_one)
    {
        double m = par.m_();
        double b = par.b_(which_one);
        double D = par.d;        
        double L = par.L_(which_one);
        std::complex<double> res = 0;
        std::complex<double> theta = theta_(lambda);

        double term1 = exp(m*b)*L;
        std::complex<double> term2 = exp(-theta*b) * normale_cdf( theta*sqrt(D) - ( b/sqrt(D) ) ) +
         exp(theta*b) * normale_cdf( -theta*sqrt(D) - ( b/sqrt(D) ) );

        option par_tilde = option(par.r,par.delta,par.sigma,1,par.T,par.K,0,
        par.TAG1,par.TAG2,par.TAG3,par.TAG4,par.d,par.a1,1,par.a2,1);

        res = pdoc_hat_inf(par_tilde,lambda,which_one)*term1*term2; 

        return res * term1 * term2;
    }

std::complex<double> puic_hat_inf(option par, std::complex<double> lambda,int which_one)
    {
        return sc_hat(par,lambda,which_one) - puoc_hat_inf(par,lambda,which_one);
    }

// Gestion DPOC-DPIC ////////////////////////////////////////////////////////////////////////////////////////////////////

std::complex<double> pdic_hat(option par, std::complex<double> lambda,int which_one)
    {
        double b = par.b_(which_one);
        if (b<=0)
            {
                return pdic_hat_inf(par,lambda,which_one);
            }
        else
            {
                return pdic_hat_sup(par,lambda,which_one);                
            }
    }

std::complex<double> puic_hat(option par, std::complex<double> lambda,int which_one)
    {
        double b = par.b_(which_one);
        if (b>=0)
            {
                return puic_hat_sup(par,lambda,which_one);
            }
        else
            {
                return puic_hat_inf(par,lambda,which_one);                
            }
    }

std::complex<double> a_hat(option par, std::complex<double> lambda)
    {
        std::complex<double> res;

        res = first_minus(par,lambda)*second_plus(par,1,lambda)*puic_hat_sup(par,lambda,2) +
        first_plus(par,lambda)*second_minus(par,2,lambda)*pdic_hat_inf(par,lambda,1);
        return res;
    }

std::complex<double> dpoc_hat(option par, std::complex<double> lambda, int which_one)
    {
        std::complex<double> res;
        res = sc_hat(par,lambda,0) - pdic_hat(par,lambda,1)
        -puic_hat(par,lambda,2) + a_hat(par,lambda);

        return res;
    }

std::complex<double> dpic_hat(option par, std::complex<double> lambda, int which_one)
    {
        std::complex<double> res;
        res = pdic_hat(par,lambda,1)+puic_hat(par,lambda,2) - a_hat(par,lambda);

        return res;
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Ces fonctions inversent la transformee de laplace

template<class T>
double somme(double t, std::complex<double> alpha, T func,option par,int which_one, int indice)
    {
        std::complex<double> I (0,1);
        std::complex<double> res = 0;
        double PI  = 3.141592653589793238463;

        res = exp(alpha * t) * func(par,alpha,which_one) / (2*t);

        for (double k = 1; k<= indice; k++)
            {
                res = res+ exp(alpha * t) * pow(-1,k) * real( func(par,alpha + I * PI * k / t,which_one) ) / t ;
            } 

        return real(res);
    }

template<class T>
double laplace(option par, int which_one, T func)
    {
        double alpha = 9.2 / par.T ;
        double res=0;
        double t = par.T;

        double m = par.m_();
        double sigma = par.sigma;

        if (alpha <= (m+sigma)*(m+sigma)*0.5  )
            {
                cout << "Alpha n'est pas assez grand";
                return 0;
            };

        res = somme(t,alpha,func,par,which_one,2000);
        // 2000 premiers termes suffisant pour avoir la convergence

        return res;
    }

