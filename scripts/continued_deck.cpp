// Test of Deck-type loop in pi1 production with realistic spins
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

void continued_deck()
{
    using namespace triangleTools;
    using triangleTools::complex;

    // Masses
    double  mu2   = norm(M_PION), sig = norm(0.770);
    double  m3pi2 = norm(1.4),     t  = -0.1;
    complex ieps  = I*1E-5;

    arguments args(id::log_divergent);
    args._external = {mu2, t,   m3pi2+ieps};
    args._internal = {mu2, sig, mu2};

    feynman    fT(1E8);
    dispersive dT(20);

    int N = 500;
    double min = -0.5, max = 3.25;
    
    std:vector<double> x, fim, dim, dre;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._internal = {mu2, xi, mu2};

        complex fTxi = fT(args);
        complex dTxi = dT.discontinuity(args);

        x.push_back(xi);
        fim.push_back(imag(fTxi));
        dre.push_back(real(dTxi));
        dim.push_back(imag(dTxi));
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_legend(0.32, 0.7);
    p.add_header("#it{t} = #minus 0.1, #it{m}_{3#pi}^{2} = (1.4)^{2}");
    p.add_curve(x, dre, solid(jpacColor::Blue, "Real"));
    p.add_curve(x, dim, solid(jpacColor::Red, "Imaginary"));
    p.add_curve(x, fim, dashed(jpacColor::DarkGrey, "Im #it{T}(#it{t}, #it{m}^{2}_{3#pi} ; #sigma)"));
    p.set_labels("#sigma  [GeV^{2}]", "#Delta(#it{t}, #it{m}_{3#pi}^{2} #; #sigma)");
    p.shade_region({4*mu2, norm(sqrt(m3pi2)-sqrt(mu2))});
    p.add_vertical(norm(sqrt(m3pi2)+sqrt(mu2)));
    p.save("T1_deck.pdf");
};