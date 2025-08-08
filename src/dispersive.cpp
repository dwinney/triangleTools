// Evaluate triangles with dispersion relations of a single variable
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "dispersive.hpp"

namespace triangleTools
{
    complex dispersive::evaluate(const arguments & args)
    {
        // Save all the masses
        save_args(args);
        if (is_zero(_M3)) { return 0.; }

        // Integration is with respect to the _M3 variable
        // i.e. cut along the m1-m2 intermediate state
        double low  = std::norm(csqrt(_m1) + csqrt(_m2));
        double pth  = std::norm(csqrt(_M1) - csqrt(_M2));
        double high = std::numeric_limits<double>::infinity();

        if (pth <= low) return disperse(low, high);
        return disperse(low, pth) + disperse(pth, high);
    };

    complex dispersive::disperse(double low, double high)
    {
        using namespace boost::math::quadrature;
        double inf = std::numeric_limits<double>::infinity();

        // Subtraction for Cauchy trick
        complex subtraction = (real(_M3)>=low) ? rho(_M3)*discontinuity(_M3) : 0.;

        complex log_term = (high == inf) ?                 - log(1-_M3/low)
                                         : log(1-_M3/high) - log(1-_M3/low);
        log_term *= subtraction;

        // Remaining principal value integral
        auto integrand = [&](double x)
        {
            complex num = rho(x)*discontinuity(x) - subtraction;
            return num*(_M3/x)/(x-_M3);
        };
        complex integral = gauss_kronrod<double,61>::integrate(integrand, low, high, _depth, 1.E-9, NULL);
        
        return (integral+log_term);
    };
    
    complex dispersive::discontinuity(complex x)
    {
        switch (_id)
        {
            case id::convergent:    return Q0(x);
            case id::log_divergent:
            {
                complex a, b;
                a = 2;
                b = (x*x-x*(_M1+_M2+_m1+_m2)+(_m1-_m2)*(_M1-_M2))/x;

                complex term_1 = (a*a*Q2(x)+2*a*b*Q1(x)+b*b*Q0(x))/norm(p(x));
                complex term_2 = norm(q(x))*Q0(x);
                return term_2 - term_1;
            };
            default: return NaN<complex>();
        };
    };
};