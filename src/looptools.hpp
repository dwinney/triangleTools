// Evaluate triangles with by simply calling LoopTools functions
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef LOOPTOOLS_HPP
#define LOOPTOOLS_HPP

#include "clooptools.h"
#include <array>
#include "interface.hpp"
#include "constants.hpp"
#include "utilities.hpp"

namespace triangleTools
{
    class looptools
    {
        public: 

        looptools(){};

        // Structure is simple because we only need  
        // to call the appropriate C0iC function
        complex evaluate(const args & args); 
        inline complex operator()(const args & args){ return evaluate(args); };
    };
};

#endif