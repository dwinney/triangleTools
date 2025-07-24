// Pass all masses in a general struct
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include "constants.hpp"

namespace triangleDeck 
{
    struct masses
    {
        std::array<complex,3> _internal;
        std::array<complex,3> _external;
    };
};