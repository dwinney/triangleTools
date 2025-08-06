#include "constants.hpp"
#include "utilities.hpp"
#include "feynman.hpp"
#include "looptools.hpp"
#include "dispersive.hpp"
#include "plotter.hpp"
#include "plot.hpp"

void eta_decay()
{
    using namespace triangleTools;
    using triangleTools::complex;

    // Masses
    double  mu2   = norm(M_PION);
    double  M2    = norm(0.540),  t  = norm(0.770);
    complex ieps  = I*1E-5;

    arguments args(id::convergent);
    args._external = {M2,  mu2, 0};
    args._internal = {mu2, mu2, t};

    feynman    fT(1E8);
    dispersive dT(20);

    int N = 10;
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
    p.set_labels("#it{m}_{3#pi}^{2}  [GeV^{2}]", "#it{T}_{0}");
    p.save("t0.pdf");
};