// Test of Deck-type loop in pi1 production with all particles as scalars
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

void scalar_deck()
{
    using namespace triangleTools;
    using triangleTools::complex;

    // Masses
    double  mu2   = norm(M_PION), sig = norm(0.770);
    double  m3pi2 = norm(1.4),     t  = -0.1;
    complex ieps  = I*1E-5;

    arguments args(id::convergent);
    args._external = {mu2, t,   EPS};
    args._internal = {mu2, sig, mu2};

    looptools  lT;
    feynman    fT(1E8);
    dispersive dT(20);

    int N = 100;
    double min = EPS, max = 2.4;
    
    complex lsub = lT(args), fsub = fT(args), dsub = dT(args);
    std:vector<double> x, lre, lim, fre, fim, dre, dim;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {mu2, t, xi+ieps};

        complex lTxi = lT(args) - lsub;
        complex fTxi = fT(args) - fsub;
        complex dTxi = dT(args) - dsub;

        x.push_back(xi);
        lre.push_back(real(lTxi)); lim.push_back(imag(lTxi));
        fre.push_back(real(fTxi)); fim.push_back(imag(fTxi));
        dre.push_back(real(dTxi)); dim.push_back(imag(dTxi));
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_legend(0.7, 0.45);
    p.add_header("#minus #it{t} = 0.1, #sigma = #it{m}_{#rho}^{2}");
    p.add_curve(x, lim, solid(jpacColor::Green, "LoopTools"));
    p.add_curve(x, lre, solid(jpacColor::Green));
    p.add_curve(x, fim, dotted(jpacColor::Red, "Feynman"));
    p.add_curve(x, fre, dotted(jpacColor::Red));
    p.add_curve(x, dre, dashed(jpacColor::Blue, "Dispersive"));
    p.add_curve(x, dim, dashed(jpacColor::Blue));
    p.set_labels("#it{m}_{3#pi}^{2}  [GeV^{2}]", "#it{T}_{0}(#it{t}, #it{m}^{2}_{3#pi} #; #sigma)");
    p.save("T0_deck.pdf");
};