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
        using namespace boost::math::quadrature;

        // Save all the masses
        save_args(args);
        if (is_zero(_M3)) { return 0.; }

        // Integration is with respect to the _M3 variable
        // i.e. cut along the m1-m2 intermediate state
        double low  = std::norm(csqrt(_m1) + csqrt(_m2));
        double high = std::numeric_limits<double>::infinity();

        // Subtraction for Cauchy trick
        complex subtraction = (real(_M3)>=low) ? rho(_M3)*discontinuity(_M3) : 0.;
        complex log_term    = subtraction*log(1-_M3/low);

        // Remaining principal value integral
        auto integrand = [&](double x)
        {
            complex num = rho(x)*discontinuity(x) - subtraction;
            return num*(_M3/x)/(x-_M3);
        };
        complex integral = gauss_kronrod<double, 61>::integrate(integrand, low, high, _depth, 1.E-9, NULL);
        
        return (integral-log_term);
    };

    complex dispersive::discontinuity(complex x)
    {
        switch (_id)
        {
            case id::convergent: return Q0(x);
            default: return NaN<complex>();
        };
    };
};