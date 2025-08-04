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
        if (is_zero(_s)) { return 0.; }

        // Integration is with respect to the _M3 variable
        // i.e. cut along the m1-m2 intermediate state
        double low  = std::norm(sqrt(real(_m1)) + sqrt(real(_m2)));
        double high = std::numeric_limits<double>::infinity();

        // Subtraction for Cauchy trick
        complex subtraction = (_s >= low) ? rho(_s)*discontinuity(_s) : 0.;
        complex log_term    = subtraction*log(1 - _M3/low);

        // Remaining principal value integral
        auto integrand = [&](double x)
        {
            complex num = rho(x)*discontinuity(x) - subtraction;
            return num*(_s/x)/(x - _M3);
        };
        complex integral = gauss_kronrod<double, 61>::integrate(integrand, low, high, _depth, 1.E-9, NULL);
        
        return (integral-log_term);
    };

    complex dispersive::discontinuity(complex s)
    {
        switch (_id)
        {
            case id::convergent: return Q0();
            default: return NaN<complex>();
        };
    };
};