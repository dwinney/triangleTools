// Spin projection functions, Q_j(s,t)
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "projection_function.hpp"

// ---------------------------------------------------------------------------
// Evaluate the cross-channel exchange projected amplitude
// Q_{jjp}(s,t)
std::complex<double> triangleTools::projection_function::eval(double s, double t)
{
    set_energies(s, t);

    std::complex<double> result = 0.;

    auto error = [&] ()
    {
        std::cout << "\nError! projection_function:";
        std::cout << " j = " << std::to_string(qns->j);
        std::cout << " and j' = " << std::to_string(qns->jp);
        std::cout << " (code " << std::to_string(qns->id()) << ")";
        std::cout << " combination not available. Quitting... \n";
        exit(1);
    };

    // debug(mDec2);
    switch (qns->id())
    {

        // s wave, scalar exchange
        case 0:
        {
            result = Q(qns->l);
            break;
        }

        // s wave, vector exchange
        case 1:
        {
            result  = Q(qns->l+1);
            result += (2.*s - mDec2 - 3.*M_PION*M_PION) * Q(qns->l);
            break;
        }

        // p - wave, scalar exchangze
        case 10:
        {
            result = 2. * Q(qns->l+1);
            result += (s - mDec2 - 3.*M_PION*M_PION) * Q(qns->l);
            result /= psqr();
            break;
        }

        // p - wave, vector exchange
        case 11:
        {
            result = 2. * Q(qns->l+2);
            result += (5.*s + - 3. * mDec2 - 9. * M_PION*M_PION) * Q(qns->l+1);
            result += (2.*s*s - 3.*mDec2*s - 9.*M_PION*M_PION*s + mDec2*mDec2 + 6.*mDec2*M_PION*M_PION + 9.*M_PION*M_PION*M_PION*M_PION) * Q(qns->l);
            result /= psqr();
            break;
        }

        // d-wave scalar exchange
        case 20:
        {
            std::complex<double> term1, term2;

            // z^2
            term1  = 4.* Q(qns->l+2);
            term1 += (4. * s - 4. * mDec2 - 12. * M_PION*M_PION) * Q(qns->l+1);
            term1 += (s*s - 2.*mDec2*s + mDec2*mDec2 - 6.*M_PION*M_PION*s + 6.*M_PION*M_PION*mDec2 + 9.*M_PION*M_PION*M_PION*M_PION) * Q(qns->l);
            term1 /= psqr() * psqr();

            // 1
            term2  = Q(qns->l);
            term2 *= qsqr() / psqr();

            result = 3. * term1;
            result -= term2;

            break;
        }

        // a1 lam = 0 lamp = 0, s-wave, scalar exchange
        case 10000:
        {
            result  = s * Q(qns->l+1);
            result += (mDec2 - M_PION*M_PION) * Q(qns->l+1);
            result += (s - mDec2 - M_PION*M_PION) * (mDec2 - M_PION*M_PION) * Q(qns->l);

            break;
        };
            
        // Omega case
        case -11111:
        {
            std::complex<double> term1, term2;

            // q^2 * z^2
            term1  = 0.;
            term1 += 4.* Q(qns->l+2);
            term1 += (4. * s - 4. * mDec2 - 12. * M_PION*M_PION) * Q(qns->l+1);
            term1 += (s*s - 2.*mDec2*s + mDec2*mDec2 - 6.*M_PION*M_PION*s + 6.*M_PION*M_PION*mDec2 + 9.*M_PION*M_PION*M_PION*M_PION) * Q(qns->l);
            term1 *= 1. / psqr();
            // term1 /= s;

            // q^2 * 1
            term2  = Q(qns->l);
            term2 *= qsqr();

            result = (term2 - term1)/4;
            break;
        };

        case 11010:
        {
            std::complex<double> term1, term2;

            // q^2 * z^2
            term1  = 4.* Q(qns->l+2);
            term1 += (4. * s - 4. * mDec2 - 12. * M_PION*M_PION) * Q(qns->l+1);
            term1 += (s*s - 2.*mDec2*s + mDec2*mDec2 - 6.*M_PION*M_PION*s + 6.*M_PION*M_PION*mDec2 + 9.*M_PION*M_PION*M_PION*M_PION) * Q(qns->l);
            term1 *= 1. / psqr();

            // q^2 * 1
            term2  = Q(qns->l);
            term2 *= qsqr();

            result = term2 - term1;
            result *= - sqrt(mDec2);
            break;
        };

        case 11011:
        {
            auto F = [&](int i)
            {
                std::complex<double> term1, term2;
                // q^2 * z^2
                term1  = 4.* Q(qns->l+2 + i);
                term1 += (4. * s - 4. * mDec2 - 12. * M_PION*M_PION) * Q(qns->l+1 + i);
                term1 += (s*s - 2.*mDec2*s + mDec2*mDec2 - 6.*M_PION*M_PION*s + 6.*M_PION*M_PION*mDec2 + 9.*M_PION*M_PION*M_PION*M_PION) * Q(qns->l + i);
                term1 *= 1. / psqr();

                // q^2 * 1
                term2  = Q(qns->l+i);
                term2 *= qsqr();

                return term2 - term1; 
            }; 

            result = sqrt(mDec2) * (F(1) + (2.*s -mDec2-3.*M_PION*M_PION)*F(0)); 

            break;
        };

        case -420:
        {
            std::complex<double> term1, term2;

            qns->l = 1;

            // q^2 * z^2
            term1  = 4.* Q(qns->l+2);
            term1 += (4. * s - 4. * mDec2 - 12. * M_PION*M_PION) * Q(qns->l+1);
            term1 += (s*s - 2.*mDec2*s + mDec2*mDec2 - 6.*M_PION*M_PION*s + 6.*M_PION*M_PION*mDec2 + 9.*M_PION*M_PION*M_PION*M_PION) * Q(qns->l);
            term1 *= 1. / psqr();

            // q^2 * 1
            term2  = Q(qns->l);
            term2 *= qsqr();

            result = term2 - term1;
            break;
        };

        default: error();
    }

    // result /= pow(t, double(qns->l));

    return result;
};

// ---------------------------------------------------------------------------
// Angular projection Q kernel functions
// These are of the form:
// Q_n(s,t) = 1/Kacser(s) * \int_{t_minus}^{t_plus}dx  x^n / (x - t - ieps)
std::complex<double> triangleTools::projection_function::Q_0()
{
  std::complex<double> result;
  result  = log(t - ieps - t_minus());
  result -= log(t - ieps - t_plus());

  result /= Kacser();
  return result;
};

std::complex<double> triangleTools::projection_function::Q(int k)
{
  switch (k)
  {
    case 0:
    {
      return Q_0();
    }
    case 1:
    {
      return t * Q_0() - 1.;
    }
    case 2:
    {
      return t*t * Q_0() - t - 0.5 * (pow(t_plus(), 2.) - pow(t_minus(), 2.)) / Kacser();
    }
    case 3:
    {
      return t*t*t * Q_0() - t*t - t*(pow(t_plus(), 2.) - pow(t_minus(), 2.)) / (2.*Kacser()) - (pow(t_plus(), 3.) - pow(t_minus(), 3.)) / (3.*Kacser());
    }
    default:
    {
     std::cout << "\nNot enough Q's!!!\n";
     exit(1);
    }
  }
};

// ---------------------------------------------------------------------------
// Usual Kallen triangle function
std::complex<double> triangleTools::Kallen(std::complex<double> x, std::complex<double> y, std::complex<double> z)
{
  return x * x + y * y + z * z - 2. * (x * z + y * z + x * y);
};

// ---------------------------------------------------------------------------
// Kacser function which includes the correct analytic structure of
// product of breakup momenta, p(s) * q(s)
std::complex<double> triangleTools::projection_function::psqr()
{
    std::complex<double> result;
    result = pow(sqrt(s) + M_PION, 2.) - mDec2 - ieps;
    result *= pow(sqrt(s) - M_PION, 2.) - mDec2 - ieps;
    result /= s;

    return result;
};

std::complex<double> triangleTools::projection_function::qsqr()
{
    std::complex<double> result;
    result  = Kallen(s, M_PION*M_PION, M_PION*M_PION);
    result /= s;

    return result;
};

std::complex<double> triangleTools::projection_function::Kacser()
{
  std::complex<double> result;

  result  = sqrt(pow(sqrt(s) + M_PION, 2.) - mDec2 - ieps);
  result *= sqrt(pow(sqrt(s) - M_PION, 2.) - mDec2 - ieps);
  result *= sqrt(Kallen(s, M_PION*M_PION, M_PION*M_PION)) / s;

  return result;
};

// Ratio of agular momentum barrier factors that are removed when partial wave projecting
// 1 / p^2(s)
std::complex<double> triangleTools::projection_function::barrier_ratio(int ell)
{
  if (ell == 0)
  {
    return 1.;
  }
  else
  {
    std::complex<double> result;
    result  = s;
    result /= pow(sqrt(s) + M_PION, 2.) - mDec2 - ieps;
    result /= pow(sqrt(s) - M_PION, 2.) - mDec2 - ieps;

    return pow(result, xr * double(ell));
  }
};

// ---------------------------------------------------------------------------
// complex Bounds of integration
std::complex<double> triangleTools::projection_function::t_minus()
{
  return (mDec2 + ieps) + M_PION*M_PION - (s + mDec2 + ieps - M_PION*M_PION) / 2. - Kacser() / 2.;
};

std::complex<double> triangleTools::projection_function::t_plus()
{
  return (mDec2 + ieps) + M_PION*M_PION - (s + mDec2 + ieps - M_PION*M_PION) / 2. + Kacser() / 2.;
};