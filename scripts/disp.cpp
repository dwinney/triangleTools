#include "constants.hpp"
#include "utilities.hpp"
#include "dispersive.hpp"
#include "plotter.hpp"
#include "plot.hpp"

void disp()
{
    using namespace triangleTools;
    using triangleTools::complex;
    
    // Masses 
    double  mu2   = norm(M_PION);
    double  M2    = norm(0.780),  t  = norm(0.770);

    arguments args(id::convergent);
    args._external = {M2,  mu2, 0};
    args._internal = {mu2, mu2, t};

    dispersive dT;

    int N = 50;
    double min = 0, max = 1.0;
    
    std:vector<double> x, lre, lim, fre, fim, dre, dim;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {M2, mu2, xi+IEPS};

        complex dTxi = dT(args);
        print(xi, dTxi);        
        x.push_back(xi);
        dre.push_back(real(dTxi)); dim.push_back(imag(dTxi));
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.add_curve(x, dre, dotted(jpacColor::Blue));
    p.add_curve(x, dim, dotted(jpacColor::Red));
    p.set_labels("#it{m}_{3#pi}^{2}  [GeV^{2}]", "#it{T}_{0}");
    p.add_vertical({norm(2*M_PION), norm(0.780-M_PION), norm(0.780+M_PION)});
    p.save("t0.pdf");
};