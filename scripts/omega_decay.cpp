// Test of loop function which appears in ω -> 3π P-wave decay
//
// Author:       Daniel Winney (2025)
// Email:        daniel.winney@gmail.com
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "utilities.hpp"
#include "feynman.hpp"
#include "looptools.hpp"
#include "dispersive.hpp"
#include "plotter.hpp"
#include "plot.hpp"

void omega_decay()
{
    using namespace triangleTools;
    using triangleTools::complex;

    // Masses
    double  mu2   = norm(M_PION);
    double  M2    = norm(0.780),  t  = norm(0.770);
    complex ieps  = I*1E-5;

    arguments args(id::log_divergent);
    args._external = {M2,  mu2, 0};
    args._internal = {mu2, mu2, t};

    feynman    fT(1E8);
    dispersive dT(20);

    int N = 100;
    double min = EPS, max = 1.0;
    
    complex fsub = fT(args), dsub = dT(args);
    std:vector<double> x, fre, fim, dre, dim;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {M2, mu2, xi+ieps};

        complex fTxi = fT(args) - fsub;
        complex dTxi = dT(args) - dsub;

        print(xi, fTxi);
        
        x.push_back(xi);
        fre.push_back(real(fTxi)); fim.push_back(imag(fTxi));
        dre.push_back(real(dTxi)); dim.push_back(imag(dTxi));
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.add_curve(x, fim, solid(jpacColor::DarkGrey));
    p.add_curve(x, fre, solid(jpacColor::DarkGrey));
    p.add_curve(x, dre, dashed(jpacColor::Blue));
    p.add_curve(x, dim, dashed(jpacColor::Red));
    p.set_labels("#it{m}_{3#pi}^{2}  [GeV^{2}]", "#it{T}_{1}");
    p.save("T1_decay.pdf");
};