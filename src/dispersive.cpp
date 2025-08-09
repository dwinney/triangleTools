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
        using boost::math::quadrature::gauss_kronrod;

        // Save all the masses
        save_args(args);
        if (is_zero(_M3)) { return 0.; }

        // Integration is with respect to the _M3 variable
        // i.e. cut along the m1-m2 intermediate state
        double low  = std::norm(csqrt(_m1) + csqrt(_m2));
        double high = std::numeric_limits<double>::infinity();

        // Subtraction for Cauchy trick
        complex subtraction = (real(_M3)>=low) ? discontinuity(_M3) : 0.;
        complex log_term    = subtraction*log(1-_M3/low);

        // Remaining principal value integral
        auto    integrand = [&](double x){ return (discontinuity(x)-subtraction)*(_M3/x)/(x-_M3); };
        complex integral  = gauss_kronrod<double,61>::integrate(integrand, low, high, _depth, 1.E-9, NULL);
        
        return (integral-log_term);
    };
    
    complex dispersive::discontinuity(complex x)
    {
        switch (_id)
        {
            case id::convergent: return rho(x)*Q0(x);
            case id::log_divergent:
            {
                complex px = p(x), qx = q(x);
                complex b = x-_M1-_M2-_m1-_m2+(_m1-_m2)*(_M1-_M2)/x;
                complex term_1 = qx*qx*Q0(x);
                complex term_2 = (Q2(x)+b*Q1(x)+b*b/4*Q0(x))/4/px/px;
                return rho(x)*(term_1-term_2);
            };
            default: return NaN<complex>();
        };
    };
};