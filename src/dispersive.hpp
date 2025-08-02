// Evaluate triangles with dispersion relations of a single variable
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef DISPERSIVE_HPP
#define DISPERSIVE_HPP

#include "utilities.hpp"
#include "constants.hpp"
#include "interface.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace triangleTools
{
    class dispersive
    {
        public: 

        // Use all default settings
        dispersive(){};

        // Change integration settings
        dispersive(double depth) : _depth(depth) {};

        // Evaluate the dispersion integral along the m1-m2 cut
        complex evaluate(const args & args); 
        inline complex operator()(const args & ms){ return evaluate(ms); };
        
        // private:

        // Masses squared so as to not need to pass around
        complex _M1, _M2, _M3, _m1, _m2, _m3;

        // Integration options 
        double _depth = 0;
        
        //--------------------------------------------------------------
        // Momenta and things of that sort

        // Triangle function for all the momenta we'll need
        inline complex kallen(complex x, complex y, complex z)
        {
            return x*x + y*y + z*z - 2*(x*y + y*z + x*z);
        };

        // Two-body phase in terms of masses m2 and m3
        inline complex rho(complex s)
        {
            return csqrt(kallen(s, _m1, _m2))/s;
        };

        // Incoming momenta squared (M1 - M2 cm-momentum)
        inline complex psqr()
        {
            complex result;
            result  = pow(csqrt(_M3)+csqrt(_M2), 2) - _M1;
            result *= pow(csqrt(_M3)-csqrt(_M2), 2) - _M1; 
            result /= _M3;
            return result;
        };

        // Outgoing momenta squared (m1 - m2 cm-momentum)
        inline complex qsqr()
        {
            complex result;
            result  = pow(csqrt(_M3)+csqrt(_m2), 2) - _m1;
            result *= pow(csqrt(_M3)-csqrt(_m2), 2) - _m1; 
            result /= _M3;
            return result;
        };

        // Product of momenta 
        inline complex kacser(){ return csqrt(psqr() * qsqr()); };

        // Bounds of momentum transfer 
        inline complex t(double z)
        {
            return _M1 + _M2 - (_M3 + _M1 - _M2)*(_M3 + _m1 - _m2)/2 - z*kacser()/2;
        };

        // Discontinuity across the cut
        inline complex discontinuity(complex s){ return 0.; };
    };
};

#endif