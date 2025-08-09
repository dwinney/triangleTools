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
        
        // save masses 
        inline void save_args(const arguments & args)
        {
            _M1  = args._external[0], _M2 = args._external[1], _M3 = args._external[2];
            _m1  = args._internal[0], _m2 = args._internal[1], _m3 = args._internal[2];
            _id  = args._id;
        };

        // Bounds of momentum transfer 
        inline complex t(complex x, double z)
        {
            return _M1+_m1-(x+_M1-_M2)*(x+_m1-_m2)/x/2.+z*kacser(x)/2.;
        };

        inline complex discontinuity(const arguments & args)
        {
            save_args(args);
            return discontinuity(_M3); 
        };

        private:
        
        // Masses squared so as to not need to pass around
        complex _M1, _M2, _M3, _m1, _m2, _m3;
        id _id;
        
        
        // Integration options 
        double _depth = 0;
        
        //--------------------------------------------------------------
        // Momenta and things of that sort

        // Triangle function for all the momenta we'll need
        inline complex kallen(complex x, complex y, complex z)
        {
            return x*x + y*y + z*z - 2*(x*y + y*z + x*z);
        };

        // Incoming momenta squared (M1 - M2 cm-momentum)
        inline complex p     (complex x){ return csqrt(kallen(x,_M1,_M2))/2/csqrt(x); };
        // Outgoing momenta squared (m1 - m2 cm-momentum)
        inline complex q     (complex x){ return csqrt(kallen(x,_m1,_m2))/2/csqrt(x); };
        // Two-body phase in terms of masses m1 and m2
        inline complex rho   (complex x){ return 2*q(x)/csqrt(x); };
        // Product of momenta 
        inline complex kacser(complex x){ return 4*p(x)*q(x); };

        //--------------------------------------------------------------
        // Legendres of the second kind (put in terms of t not z_t)

        inline complex Q0(complex x){ return log((_m3-t(x,-1))/(_m3-t(x,+1)))/kacser(x); };
        inline complex Q1(complex x){ return _m3*Q0(x) - 1; };
        inline complex Q2(complex x){ return _m3*Q1(x)-t(x,0); }

        //--------------------------------------------------------------
        // Finally the disc across the cut of our triangle

        complex discontinuity(complex s);
    };
};

#endif