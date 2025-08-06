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
    args._external = {mu2, t,   0. };
    args._internal = {mu2, sig, mu2};

    feynman    fT(1E8);
    dispersive dT(20);

    int N = 10;
    double min = EPS, max = 2.4;
    
    complex fsub = fT(args), dsub = dT(args);
    std:vector<double> x, fre, fim, dre, dim;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {mu2, t, xi+ieps};

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
    p.set_labels("#it{m}_{3#pi}^{2}  [GeV^{2}]", "#it{T}_{0}");
    p.save("t0.pdf");
};