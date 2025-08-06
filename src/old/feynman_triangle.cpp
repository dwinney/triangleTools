// This object defines a Triangle amplitude, i.e. the rescattering
// diagram associated with an intermediate t-channel exchange.
//
// Evaluated specifically with the Feynman method
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "feynman_triangle.hpp"

// ---------------------------------------------------------------------------
// Evaluate the triangle assuming a fixed mass exchange with mass t
std::complex<double> triangleTools::feynman_triangle::eval(double s, double t)
{
    // Desination for the result and assosiated errors
    double val[2], err[2];

    // Integrate both x and y from 0 to 1
    double min[2] = {0., 0.};
    double max[2] = {1., 1.};

    // Fix the "masses" s and t
    integrand.set_energies(s, t);

    // TODO: Set relative errors and max calls to actual good values
    // Integrate over x and y
    hcubature(2, wrapped_integrand, &integrand, 2, min, max, 0, 0, 1e-3, ERROR_INDIVIDUAL, val, err);

    // Assemble the result as a complex double
    std::complex<double> result = val[0] + xi * val[1];
    result *= 2.; // Factor of 2 from the normalization of dF_3 integration measure

    return result;
};

// ---------------------------------------------------------------------------
// Wrapper for the feynman parameter integrands to fit into hcubature
int triangleTools::feynman_triangle::wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval)
{
  dF3_integrand* integrand = (dF3_integrand *) fdata;

  // Feynman parameters
  double x = in[0] * in[1];
  double y = in[0] * (1. - in[1]);
  double z = 1. - x - y;

  std::complex<double> result = in[0] * integrand->eval(x, y, z);

  // Split up the real andi imaginary parts to get them out
  fval[0] = std::real(result);
  fval[1] = std::imag(result);

  return 0.;
};