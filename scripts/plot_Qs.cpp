#include "constants.hpp"
#include "utilities.hpp"
#include "dispersive.hpp"
#include "plotter.hpp"
#include "plot.hpp"

void plot_Qs()
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
        
    int N = 200;
    double min = norm(2*M_PION), max = 1.1;
    std:vector<double> x, rq0, iq0, rq1, iq1;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {M2, mu2, xi+IEPS};

        dT.save_masses(args);
        complex q0 = dT.Q0(), q1 = dT.Q1();
        
        x.push_back(xi);
        rq0.push_back(real(q0)); iq0.push_back(imag(q0)); 
        rq1.push_back(real(q1)); iq1.push_back(imag(q1)); 
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_legend(0.7, 0.7);
    p.add_curve(x, rq0, solid (jpacColor::Blue, "#it{Q}_{0}"));
    p.add_curve(x, iq0, dashed(jpacColor::Blue));
    p.add_curve(x, rq1, solid (jpacColor::Red, "#it{Q}_{1}"));
    p.add_curve(x, iq1, dashed(jpacColor::Red));
    p.set_labels("#it{s}  [GeV^{2}]", "#it{Q}");
    p.save("qs.pdf");
};