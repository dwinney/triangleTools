// Evaluate triangles with a brute force integration of Feynman parameters
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef FEYNMAN_TRIANGLE_HPP
#define FEYNMAN_TRIANGLE_HPP

#include <array>
#include "args.hpp"
#include "constants.hpp"
#include "cubature.h"

namespace triangleTools
{
    // 
    class feynman_triangle
    {
        public:

        // Default triangle, dont really need to specify anything
        feynman_triangle(){};

        // Brute force feynman parameters integration
        complex evaluate(const args & ms); 
        inline complex operator()(const args & ms){ return evaluate(ms); };

        private:
            
        // Interface with cubature
        static int wrapped_integrand(unsigned ndim, const double *in, void *fdata, unsigned fdim, double *fval);
        
        // Pass specific integrand desired into the wrapper above
        struct integrand
        {
            // Save the arguments originally passed through feynman_triangle::evaluate
            void pass_args(args args)
            {
                _id = args._id;
                _M1 = args._external[0], _M2 = args._external[1], _M3 = args._external[2]; 
                _m1 = args._internal[0], _m2 = args._internal[1], _m3 = args._internal[2]; 
            };

            // Masses SQUARED of all particles involved
            complex _M1, _M2, _M3; // external 
            complex _m1, _m2, _m3; // internal 

            id _id;

            // Combination of denominators 
            complex D(double x1, double x2, double x3);

            // Filter _id 
            complex evaluate(double x1, double x2, double x3);
        };

        integrand _integrand;
    };

};

#endif