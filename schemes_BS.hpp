/////////////////////////////////////
///Projet Option Parisienne
/////////////////////////////////////
//
//Nicolas Klaeylé - Timothée Fabre

#pragma once

struct exact_bs {
    using result_type = double;
    operator result_type() const { return state; }
    template <typename TAlgo, typename TRandom> friend struct random_scheme;

    exact_bs() = default;
    exact_bs(double r, double d, double sigma, double x0, double h)
        : r(r), delta(d), sigma(sigma), state(x0), h(h) {}

    template <typename TWhiteNoise>
    double operator()(TWhiteNoise const & z, unsigned last_n) {
        return state *= exp((r - delta - 0.5*sigma*sigma) * last_n * h + sigma * sqrt(last_n * h) * z);
    }

    template <typename TWhiteNoise>
    double operator()(TWhiteNoise const & z) {
        return state *= exp((r - delta - 0.5*sigma*sigma) * h + sigma * sqrt(h) * z);
    }
public:
    double r, delta, sigma;
    double state;
    double h;
};

struct exact_bs_Qtild {
    using result_type = double;
    operator result_type() const { return state; }
    template <typename TAlgo, typename TRandom> friend struct random_scheme;

    exact_bs_Qtild() = default;
    exact_bs_Qtild(double sigma, double x0, double h)
        : sigma(sigma), state(x0), h(h) {}

    template <typename TWhiteNoise>
    double operator()(TWhiteNoise const & z, unsigned last_n) {
        return state *= exp(sigma * sqrt(last_n * h) * z);
    }

    template <typename TWhiteNoise>
    double operator()(TWhiteNoise const & z) {
        return state *= exp(sigma * sqrt(h) * z);
    }
public:
    double sigma;
    double state;
    double h;
};
