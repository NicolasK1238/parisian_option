/////////////////////////////////////
///Projet Option Parisienne
/////////////////////////////////////
//
//Nicolas Klaeylé - Timothée Fabre

#pragma once

// Fonctions "barri�re"
double barrier_exp(double a, double u_0, double t){
    return u_0 * exp(a * t);
}
double barrier_lin(double a, double u_0, double t){
    return a * t + u_0;
}

// Structure barri�re
struct barrier {
    barrier(double a, double u_0, string type_): coeff_barrier(a), initial_value(u_0), type_barrier(type_) {};

    double operator()(double time) 
        {
            if (type_barrier == "constante" || type_barrier == "affine" || 
            type_barrier == "-") {
                return barrier_lin(coeff_barrier ,initial_value, time);
            } else if (type_barrier == "exponentielle") {
                return barrier_exp(coeff_barrier ,initial_value, time);
            } else {
                //std::cout << "Problème avec l'initialisation des barrières, voir input.dat" << std::endl;
            return 0;}
        };


    friend void make_barrier(barrier &bup, barrier &bdown, double a1, double b1, double a2, double b2, 
    int n, int type_ud, string tag1, string tag2)
    {
        // on initialie les barrières en fct des tags
        // Un petit peu laid comme implémentation, mais pas moyen de faire autrement sans refaire tt le code mc
  
        if (n==1)  // mono barrière
        {
            if (type_ud==1) //Up
                {
                    bup.initial_value = b1;
                    bup.coeff_barrier = a1;
                    bup.type_barrier = tag1; 
                    bdown.initial_value = -std::numeric_limits<double>::infinity();
                    bdown.coeff_barrier = 0;
                    bdown.type_barrier = tag2; 
                } 
            if (type_ud==-1) //Down
                {
                    bdown.initial_value = b1;
                    bdown.coeff_barrier = a1;
                    bdown.type_barrier = tag1; 
                    bup.initial_value = +std::numeric_limits<double>::infinity();
                    bup.coeff_barrier = 0;
                    bup.type_barrier = tag2; 
                } 
        }
        else if (n==2) // Double barrière
        // Par convention dans le fichier d'input, en cas de double barrière, la plus petite est b1
            {
                bup.initial_value = b2;
                bup.coeff_barrier = a2;
                bup.type_barrier = tag2; 
                bdown.initial_value = b1;
                bdown.coeff_barrier = a1;
                bdown.type_barrier = tag1; 
            }
        else
            {
                std::cout << "Problème avec l'initialisation des barrières, voir input.dat" << std::endl;
            }
        
    }


    public:
        double coeff_barrier, initial_value;
        string type_barrier;
};

// Fonctions payoff
double payoff_Call(double prix_Sj, double strike){
    return((prix_Sj-strike>0) ? prix_Sj-strike : 0);
}

double payoff_Put(double prix_Sj, double strike){
    return((strike-prix_Sj>0) ? strike-prix_Sj : 0);
}

// Structure option parisienne
typedef double (* point_Payoff)(double prix_Sj, double strike);
struct parisian_Option {
    parisian_Option() = default;
    parisian_Option(double T, double r, double delta, double sigma, double K, double s_0, barrier bu, barrier bd, double delay, double p, char p_type, point_Payoff id_fp): maturity(T), interest_Rate(r),
    dividend_Rate(delta), implied_Vol(sigma), strike(K), S_0(s_0), barrier_Up(bu), barrier_Down(bd), delay(delay), price(p), parisian_Type(p_type), f_Payoff(id_fp) {};
    double operator()(double const& S) {
        return((*f_Payoff)(S,strike));
    }

    void caracteristics() {
        std::cout << "Implied volatility: " << 100*implied_Vol << " %;  Spot price: " << S_0 << ";" << std::endl;
        if (strike == S_0){
            std::cout << "Strike: ATM;  Maturity: " << maturity << " years;" << std::endl;
        } else {
            std::cout << "Strike: " << strike << "; Maturity: " << maturity << " years;" << std::endl;
        }
        std::cout << "Interest rate: " << 100*interest_Rate << " %;  Dividend rate: " << 100*dividend_Rate << " %;" << std::endl;
        std::cout << "Initial corridor: ]" << barrier_Down(0) << ", " << barrier_Up(0) << "[; " << std::endl;
        std::cout << "Final corridor: ]" << barrier_Down(maturity) << ", " << barrier_Up(maturity) << "[; " << std::endl;
        std::cout << "Delay: " << delay << " years; " << std::endl;
        if (parisian_Type=='I') {
            std::cout << "Type: IN;" << std::endl;
        } else if (parisian_Type=='O') {
            std::cout << "Type: OUT;" << std::endl;
        } else {std::cout << "Please initialize the type with 'O' or 'I'! " << std::endl;}
        std::cout << "Option price: " << price << "." << std::endl;
    }

    void change_kd(double new_strike, double new_delay) 
    {
        strike = new_strike;
        delay = new_delay;        
    }

    void change_sd(double new_sigma, double new_delay) 
    {
        implied_Vol = new_sigma;
        delay = new_delay;        
    }


    public:
        double maturity, interest_Rate, dividend_Rate, strike, delay;
        double S_0, price, implied_Vol;
        char parisian_Type;
        barrier barrier_Up, barrier_Down;
        point_Payoff f_Payoff;
};

// Structure de calcul du pay-off de l'option parisienne (utilisant la structure initialis�e "Parisian Option") pour une simulation
template <typename TRealScheme>
struct parisian_Realisation {
    parisian_Realisation(TRealScheme S, parisian_Option PO, double final_Re,  bool activ_Pay = 0) : scheme(S), Par_Opt(PO), activation_Payoff(activ_Pay), final_Realisation(final_Re) {};

    template <typename TGen>
    double operator()(TGen & gen){ // D�claration de l'op�rateur, pour it�rer dans la m�thode Monte Carlo: Renvoie 0 ou 1 suivant si le payoff est nul ou pas
        double g_k(0), tau(0);
        scheme.init(); // Initialisation du sch�ma
        if (Par_Opt.parisian_Type == 'I'){
            activation_Payoff = 0;
        } else if (Par_Opt.parisian_Type == 'O'){
            activation_Payoff = 1;
        }
        double S_k(scheme.algo.state);
        while (scheme.not_end()) {
            double S_k_1(scheme.algo.state); // Valeur de S(h * (k-1))
            double hk_1(scheme.algo.h * scheme.n);
            auto z = scheme.Z(gen);   // Simulation de la v.a. (suivant la loi caract�ristique de l'attribut "Z" de l'objet sch�ma "scheme"
            scheme = scheme.next(z); // Mise � jour du sch�ma (en mettant � jour l'algo, par ex. BS, par le biais de next) => m�j de l'attribut "state" de l'algo
            S_k = scheme.algo.state; // Valeur de S(h * k)
            double hk(scheme.algo.h * scheme.n);
            if (S_k > Par_Opt.barrier_Down(hk) && S_k < Par_Opt.barrier_Up(hk)){ // Si le sous-jacent n'a toujours pas activ� la barri�re
                g_k = hk; // On va chercher le pas d'incr�mentation dans l'algo, par ex. BS
            } else if (S_k <= Par_Opt.barrier_Down(hk)) {
                if (S_k_1 > Par_Opt.barrier_Down(hk_1)){ // Si le sous-jacent n'�tait pas dans la barri�re au pas pr�c�dent, simulation du dernier temps de passage
                    sup_stoppingtime_density dens_st(log(Par_Opt.barrier_Down(hk)/S_k)/Par_Opt.implied_Vol, log(S_k_1/S_k)/Par_Opt.implied_Vol, scheme.algo.h);
                    auto VA_st = acceptance_Rejection_Unif([dens_st](double x) { return dens_st.derivative(x); }, 0, scheme.algo.h, dens_st.max(), gen, 10); // M�thode de rejet
                    g_k = hk - VA_st;
                    //std::cout << "DERNIER TEMPS DE PASSAGE ------> [" << hk_1 << ", " << g_k << ", " << hk << "]" << std::endl;
                } else {
                    double p_k = exp(- 2 * (log(Par_Opt.barrier_Down(hk_1)) - log(S_k_1)) * (log(Par_Opt.barrier_Down(hk_1)) - log(S_k)) / (scheme.algo.h * std::pow(Par_Opt.implied_Vol,2))); // Calcul de la probabilit� de franchir la barri�re
                    std::bernoulli_distribution bool_pk(p_k);
                    bool reach_Barrier(bool_pk(gen)); // Simulation du bool�en
                    if (reach_Barrier) {
                        sup_stoppingtime_density dens_st(log(Par_Opt.barrier_Down(hk)/S_k)/Par_Opt.implied_Vol, log(S_k_1/S_k)/Par_Opt.implied_Vol, scheme.algo.h);
                        auto VA_st = acceptance_Rejection_Unif([dens_st](double x) { return dens_st.derivative(x); }, 0, scheme.algo.h, dens_st.max(), gen, 10); // M�thode de rejet
                        g_k = hk - VA_st;
                        //std::cout << "DERNIER TEMPS DE PASSAGE ------> [" << hk_1 << ", " << g_k << ", " << hk << "]" << std::endl;
                    }
                    }
            } else if (S_k >= Par_Opt.barrier_Up(hk)) { // Si le sous-jacent vient d'activer la barri�re
                if (S_k_1 < Par_Opt.barrier_Up(hk_1)){ // Si le sous-jacent n'�tait pas dans la barri�re au pas pr�c�dent, simulation du dernier temps de passage
                    sup_stoppingtime_density dens_st(log(S_k/Par_Opt.barrier_Up(hk))/Par_Opt.implied_Vol, log(S_k/S_k_1)/Par_Opt.implied_Vol, scheme.algo.h);
                    auto VA_st = acceptance_Rejection_Unif([dens_st](double x) { return dens_st.derivative(x); }, 0, scheme.algo.h, dens_st.max(), gen, 10); // M�thode de rejet
                    g_k = hk - VA_st;
                    //std::cout << "DERNIER TEMPS DE PASSAGE ------> [" << hk_1 << ", " << g_k << ", " << hk << "]" << std::endl;
                } else {
                    double p_k = exp(- 2 * (log(S_k_1) - log(Par_Opt.barrier_Down(hk_1))) * (log(S_k) - log(Par_Opt.barrier_Down(hk_1))) / (scheme.algo.h * std::pow(Par_Opt.implied_Vol,2))); // Calcul de la probabilit� de franchir la barri�re
                    std::bernoulli_distribution bool_pk(p_k);
                    bool reach_Barrier(bool_pk(gen)); // Simulation du bool�en
                    if (reach_Barrier) {
                        sup_stoppingtime_density dens_st(log(S_k/Par_Opt.barrier_Up(hk))/Par_Opt.implied_Vol, log(S_k/S_k_1)/Par_Opt.implied_Vol, scheme.algo.h);
                        auto VA_st = acceptance_Rejection_Unif([dens_st](double x) { return dens_st.derivative(x); }, 0, scheme.algo.h, dens_st.max(), gen, 10); // M�thode de rejet
                        g_k = hk - VA_st;
                        //std::cout << "DERNIER TEMPS DE PASSAGE ------> [" << hk_1 << ", " << g_k << ", " << hk << "]" << std::endl;
                    }
                    }
                }
            tau = hk - g_k; // Calcul du temps pass� dans la zone d'activation
            if (tau > Par_Opt.delay) { // V�rification si le pay-off est activ� (ou d�sactiv�), suivant si la barri�re est activante ou d�sactivante
                //std::cout<< "Barriere activee!!" << std::endl;
                if (Par_Opt.parisian_Type == 'I') {
                    activation_Payoff = 1;
                } else if (Par_Opt.parisian_Type == 'O') {
                    activation_Payoff = 0;
                }
                break;// On sort du "while" si activation (ou d�sactivation)
            }
            }
        //if (activation_Payoff==0) {
        //        std::cout<< "Barriere desactivee!!" << std::endl;
        //}
        //std::cout << scheme.n << "#######" << scheme.n_max << "#######" << scheme.algo.h << "#######" << std::endl;
        auto z_last = scheme.Z(gen);
        unsigned time_to_maturity((scheme.n_max - scheme.n) * scheme.algo.h);
        auto final_Re(scheme.go_to_last(z_last));
        double mu = -(Par_Opt.interest_Rate-Par_Opt.dividend_Rate-pow(Par_Opt.implied_Vol,2)/2)/Par_Opt.implied_Vol;
        final_Realisation = activation_Payoff * exp(- Par_Opt.interest_Rate * Par_Opt.maturity) * Par_Opt.f_Payoff(final_Re.algo.state,Par_Opt.strike)
                            * std::pow(Par_Opt.S_0/S_k, mu/Par_Opt.implied_Vol) * std::exp(-mu * std::sqrt(time_to_maturity) * z_last + std::pow(mu, 2) * Par_Opt.maturity); // Valeur du payoff
        //std::cout<<"Valeur du dernier temps de passage de la barriere: "<< g_k <<std::endl;
        //std::cout<<"Valeur du temps passe dans la barriere: "<< tau <<std::endl;
        //std::cout<<"Valeur de S_T: "<< final_Re.algo.state <<std::endl;
        //std::cout<<"Valeur du booleen: "<< activation_Payoff <<std::endl;
        //std::cout<<"Payoff en T: "<< activation_Payoff * Par_Opt(final_Re.algo.state) <<std::endl;
        return final_Realisation; // On renvoie la derni�re r�alisation
    }

    public:
        TRealScheme scheme;
        parisian_Option Par_Opt;
        bool activation_Payoff;
        double final_Realisation;
};

template <typename TRealScheme>
parisian_Realisation<TRealScheme>
make_Parisian_Realisation(TRealScheme S, parisian_Option PO, double final_Re, bool activ_Payoff = 0) {
    return { S, PO, final_Re, activ_Payoff };
}



