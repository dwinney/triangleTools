// Evaluate triangles with by simply calling LoopTools functions
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "looptools.hpp"

namespace triangleTools
{
    complex looptools::evaluate(const arguments & args)
    {
        complex result = 0;
        complex M1 = args._external[0], M2 = args._external[1], M3 = args._external[2];
        complex m1 = args._internal[1], m2 = args._internal[0], m3 = args._internal[2];

        ltini();
        setwarndigits(60);
        switch (args._id)
        {
            case id::convergent: 
            {
                result = - C0iC(cc0, M1, M2, M3, m2, m3, m1);
                 break;
            };
            case id::log_divergent: 
            {
                result = 2*C0iC(cc00, M1, M2, M3, m2, m3, m1);
                 break;
            };
            default: { result = NaN<complex>(); break; };
        }; 
        ltexi();

        return result;
    };
};