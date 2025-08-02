#include "constants.hpp"
#include "utilities.hpp"
#include "feynman.hpp"
#include "looptools.hpp"
#include "plotter.hpp"
#include "plot.hpp"

void plot_logdivergent()
{
    using namespace triangleTools;
    using triangleTools::complex;
    
    // Masses
    double  mu2   = norm(M_PION), sig = norm(0.770);
    double  m3pi2 = norm(1.4),     t  = -0.1;
    complex ieps  = I*1E-5;

    arguments args(id::log_divergent);
    args._external = {mu2, 0.,  t  };
    args._internal = {mu2, mu2, sig};

    feynman   fT;
    looptools lT;

    int N = 50;
    double min = 0, max = 2.5;
    
    complex fsub = fT(args), lsub = lT(args);
    std:vector<double> x, lre, lim, fre, fim;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {mu2, norm(xi) + ieps, t};

        complex lTxi = lT(args) - lsub;
        complex fTxi = fT(args) - fsub;
        
        print (i, xi, lTxi, fTxi);
        x.push_back(xi);
        fre.push_back(real(fTxi)); lre.push_back(real(lTxi)); 
        fim.push_back(imag(fTxi)); lim.push_back(imag(lTxi));
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.add_curve(x, fre, solid(jpacColor::Blue, "Real"));
    p.add_curve(x, fim, solid(jpacColor::Red,  "Imag"));
    p.add_curve(x, lre, dashed(jpacColor::Blue));
    p.add_curve(x, lim, dashed(jpacColor::Red));
    p.set_labels("#it{m}_{3#pi}^{2}  [GeV^{2}]", "#it{T}_{1}");
    p.save("t1.pdf");
};