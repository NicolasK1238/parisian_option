/////////////////////////////////////
///Projet Option Parisienne
/////////////////////////////////////
//
//Nicolas Klaeylé - Timothée Fabre

#pragma once

double pi(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651);

// M�thode de la s�cante pour la recherche de z�ros
template <typename TRealFunction>
double secanteR(TRealFunction const & f,
                  double y = 0, double l = 0, double r = 1,
                  double eps = 2*std::numeric_limits<double>::epsilon())
{
    double x_0 = (r - l)/2;
    double x_n_1 = x_0 - (r - l)/100;
    double x_n = x_0;
    do {
        double keep_x_n_1 = x_n_1;
        x_n_1 = x_n;
        x_n -= (f(x_n) - y) * (x_n - keep_x_n_1) / ((f(x_n) - y) - f(keep_x_n_1));
    } while (std::fabs(x_n - x_n_1) > eps);
    return x_n;
};


// Densit� du dernier temps de passage de la barri�re (premier temps de passage d'un pont brownien)
struct sup_stoppingtime_density {
    sup_stoppingtime_density(double a, double y, double T) : a_(a), y_(y), T_(T) {}
    double operator()(double t) const {
        return std::sqrt(T_/8*pi)*std::pow((T_-t)*t, -1.5)*(std::exp(-2*a_*(a_-y_)/T_)*(a_*T_-y_*t)*std::exp(-std::pow(-a_*T_+(2*a_-y_)*t, 2)/(2*T_*(T_-t)*t))
                                                                +(a_*T_-(2*a_-y_)*t)*std::exp(-std::pow(-a_*T_+y_*t, 2)/(2*T_*(T_-t)*t)));
    }
    double value(double t) const {
        return std::sqrt(T_/8*pi)*std::pow((T_-t)*t, -1.5)*(std::exp(-2*a_*(a_-y_)/T_)*(a_*T_-y_*t)*std::exp(-std::pow(-a_*T_+(2*a_-y_)*t, 2)/(2*T_*(T_-t)*t))
                                                                +(a_*T_-(2*a_-y_)*t)*std::exp(-std::pow(-a_*T_+y_*t, 2)/(2*T_*(T_-t)*t)));
    }
    double derivative(double t) const {
        return -std::sqrt(T_/32*pi)*(T_-2*t)*std::pow((T_-t)*t, -2.5)*(std::exp(-2*a_*(a_-y_)/T_)*(a_*T_-y_*t)*std::exp(-std::pow(-a_*T_+(2*a_-y_)*t, 2)/(2*T_*(T_-t)*t))+(a_*T_-(2*a_-y_)*t)*std::exp(-std::pow(-a_*T_+y_*t, 2)/(2*T_*(T_-t)*t)))
                + std::sqrt(T_/8*pi)*std::pow((T_-t)*t, -1.5)*(std::exp(-2*a_*(a_-y_)/T_)*(-y_+(a_*T_-y_*t)*(2*y_*(-a_*T_+(2*a_-y_)*t)*T_*(T_-t)*t+std::pow(-a_*T_+(2*a_-y_)*t, 2)*(T_-2*t)*T_)/(2*std::pow(T_*(T_-t)*t, 2)))*std::exp(-std::pow(-a_*T_+(2*a_-y_)*t, 2)/(2*T_*(T_-t)*t))
                                                         + (y_+(a_*T_-(2*a_-y_)*t)*(std::pow(-a_*T_+y_*t, 2)*(T_-2*t)*T_+2*y_*(a_*T_-y_*t)*T_*(T_-t)*t)/(2*std::pow(T_*(T_-t)*t, 2)))*std::exp(-std::pow(-a_*T_+y_*t, 2)/(2*T_*(T_-t)*t)));
    }
    double mode() const { return secanteR([this](double x) { return derivative(x); }, 0, 0, T_); }
    double max() const { return value(secanteR([this](double x) { return derivative(x); }, 0, 0, T_)); }
    private:
        double a_, y_, T_;
};


// M�thode "Acceptance-rejection" pour la simulation de la densit�, via la simulation d'une loi uniforme
template <typename TRealFunction, typename TGen>
double acceptance_Rejection_Unif(TRealFunction const & f, double a, double b, double borne_sup, TGen & gen , double const_ = 1) {
    std::uniform_real_distribution<double> U_1(0.0, 1.0), U_2(a, b);
    double realisation_Test(U_1(gen));
    double realisation_Interv(U_2(gen));
    while (f(realisation_Interv)/(const_ * borne_sup) < realisation_Test) {
        realisation_Test = U_1(gen);
        realisation_Interv = U_2(gen);
    }
    return realisation_Interv;
}



