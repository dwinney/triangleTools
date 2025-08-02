// Generic structure for a triangle graph. The primary motivation of is 
// to test and match different ways to evaluate the same triangle.
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "feynman.hpp"

namespace triangleTools
{
    // Use cubature to do the 2D feynman parameter integral with 
    // brute force quadrature
    complex feynman::evaluate(const arguments & args)
    {
        // Desination for the result and assosiated errors
        double val[2], err[2];

        // Integrate both x and y from 0 to 1
        double min[2] = {0., 0.};
        double max[2] = {1., 1.};

        // Pass args to integrand
        _integrand.save_args(args);

        // Integrate over x and y
        hcubature(2, wrapped_integrand, &_integrand, 2, min, max, _maxcalls, 0, 1e-4, ERROR_INDIVIDUAL, val, err);

        // Assemble the result as a complex double
        complex result(val[0], val[1]);
        // result /= pow(4.*PI, 2.); // Rest of left-over factors from covariant loop normalization

        return result;
    };

    // Interface with cubature
    int feynman::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
    {
        integrand* feyn_integrand = (integrand*) fdata;

        // Feynman parameters
        double x = in[0] * in[1];
        double y = in[0] * (1. - in[1]);
        double z = 1. - x - y;

        complex result = in[0]*feyn_integrand->evaluate(x, y, z);

        // Split up the real and imaginary parts to get them out
        fval[0] = real(result);
        fval[1] = imag(result);

        return 0;
    };

    // Combination of denominators 
    complex feynman::integrand::D(double x1, double x2, double x3)
    {
        return + x1   *_m1 + x2   *_m2 + x3   *_m3 
               - x2*x3*_M1 - x1*x3*_M2 - x1*x2*_M3;
    };

    // Filter integrands by _id
    complex feynman::integrand::evaluate(double x1, double x2, double x3)
    {
        switch (_id)
        {
            case id::convergent:    return 1./D(x1,x2,x3);
            case id::log_divergent: return - log(D(x1,x2,x3));
            default: return NaN<double>();
        };

        return NaN<double>();
    };
};