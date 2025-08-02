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
        complex evaluate(const arguments & args); 
        inline complex operator()(const arguments & args){ return evaluate(args); };
        
        // private:
        
        // Masses squared so as to not need to pass around
        complex _M1, _M2, _M3, _m1, _m2, _m3;
        double  _s, _ieps;
        id _id;
        
        // save masses 
        inline void save_args(const arguments & args)
        {
            _M1  = args._external[0], _M2 = args._external[1], _M3 = args._external[2];
            _m1  = args._internal[0], _m2 = args._internal[1], _m3 = args._internal[2];
            _s = real(_M3), _ieps = imag(_M3);
            _id = args._id;
        };
        
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
        inline complex rho(complex s){ return csqrt(kallen(s, _m1, _m2))/s; };
        // Incoming momenta squared (M1 - M2 cm-momentum)
        inline complex p(){ return csqrt(kallen(_s,_M1,_M2)+_ieps)/2./sqrt(_s); };
        // Outgoing momenta squared (m1 - m2 cm-momentum)
        inline complex q(){ return csqrt(kallen(_s,_m1,_m2)+_ieps)/2./sqrt(_s); };
        // Product of momenta 
        inline complex kacser(){ return 4.*p()*q(); };

        // Bounds of momentum transfer 
        inline complex t(double z)
        {
            return _M1+_M2-(_M3+_M1-_M2)*(_M3+_m1-_m2)/2/_M3+z*kacser()/2.;
        };

        //--------------------------------------------------------------
        // Legendres of the second kind (put in terms of t not z_t)

        inline complex Q0(){ return (log(_m3 - t(-1)) - log(_m3 - t(+1)))/kacser(); };
        inline complex Q1(){ return _m3*Q0() - 1; };

        //--------------------------------------------------------------
        // Finally the disc across the cut of our triangle

        complex discontinuity(complex s);
    };
};

#endif