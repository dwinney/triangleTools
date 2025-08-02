// Evaluate triangles with dispersion relations of a single variable
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "dispersive.hpp"

namespace triangleTools
{
    complex dispersive::evaluate(const args & args)
    {
        using namespace boost::math::quadrature;

        // Save all the masses
        _M1  = args._external[0], _M2 = args._external[1], _M3 = args._external[2];
        _m1  = args._internal[0], _m2 = args._internal[1], _m3 = args._internal[2];

        // Integration is with respect to the _M3 variable
        // i.e. cut along the m1-m2 intermediate state
        double low  = std::norm(sqrt(real(_m2)) + sqrt(real(_m3)));
        double high = std::numeric_limits<double>::infinity();

        // Subtraction for Cauchy trick
        complex subtraction = rho(_M3) * discontinuity(_M3);
        complex log_term    = - subtraction * log(1. - _M3/low);

        // Remaining principal value integral
        auto integrand = [&](double x)
        {
            complex num = rho(x) * discontinuity(x) - subtraction;
            return num / x / (x - _M3);
        };
        complex integral = gauss_kronrod<double, 61>::integrate(integrand, low, high, _depth, 1.E-9, NULL);
        
        return (integral + log_term) / PI;
    };
};