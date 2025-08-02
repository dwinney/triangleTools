// Evaluate triangles with by simply calling LoopTools functions
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef LOOPTOOLS_TRIANGLE_HPP
#define LOOPTOOLS_TRIANGLE_HPP

#include "clooptools.h"
#include <array>
#include "triangle_args.hpp"
#include "constants.hpp"
#include "utilities.hpp"

namespace triangleTools
{
    class looptools_triangle
    {
        public: 

        looptools_triangle(){};

        // Brute force feynman parameters integration
        complex evaluate(const triangle_args & ms); 
        inline complex operator()(const triangle_args & ms){ return evaluate(ms); };
    };
};

#endif