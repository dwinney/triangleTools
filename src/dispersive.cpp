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
        double low  = std::norm(sqrt(real(_m2)) + sqrt(real(_m3)));
        double high = std::numeric_limits<double>::infinity();

        // Subtraction for Cauchy trick
        complex subtraction = rho(_s)*discontinuity(_s);
        complex log_term    = -subtraction*log(1 - _M3/low);

        // Remaining principal value integral
        auto integrand = [&](double x)
        {
            complex num = rho(x)*discontinuity(x)*(_s/x) - subtraction;
            return num/(x - _M3);
        };
        complex integral = gauss_kronrod<double, 61>::integrate(integrand, low, high, _depth, 1.E-9, NULL);
        
        return (integral+log_term) / PI;
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