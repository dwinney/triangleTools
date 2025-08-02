#include "constants.hpp"
#include "utilities.hpp"
#include "feynman_triangle.hpp"

void feynman_plot()
{
    using namespace triangleTools;
    
    // Masses
    double mu2   = norm(M_PION), sig = norm(0.770);
    double m3pi2 = norm(1.4),     t  = -0.1;

    triangle_args args(id::convergent);
    args._external = {mu2, m3pi2 + I*0.1, t};
    args._internal = {mu2, mu2, sig};
    
    feynman_triangle T;
    print(T(args));
};