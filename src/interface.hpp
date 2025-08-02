// Interface to pass around masses and stuff easily
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <array>
#include "constants.hpp"

namespace triangleTools
{
    enum class id : unsigned int
    {
        convergent, log_divergent
    };

    struct args
    {
        args(id xid) : _id(xid)
        {};

        // Which topology we will evaluate
        id _id; 

        // Masses involved
        std::array<complex,3> _internal, _external;
    };
};

#endif