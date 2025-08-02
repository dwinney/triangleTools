#include "constants.hpp"
#include "utilities.hpp"
#include "dispersive.hpp"
#include "plotter.hpp"
#include "plot.hpp"

void plot_bounds()
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
    std:vector<double> x, rtp, itp, rtm, itm;
    for (int i = 0; i <= N; i++)
    {
        double xi = min + i*(max - min)/N;
        args._external = {M2, mu2, xi+IEPS};

        dT.save_masses(args);
        complex tp = dT.t(+1), tm = dT.t(-1);
        
        x.push_back(xi);
        rtp.push_back(real(tp)); itp.push_back(imag(tp)); 
        rtm.push_back(real(tm)); itm.push_back(imag(tm)); 
    };

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_legend(0.3, 0.2);
    p.add_curve(x, rtp, solid (jpacColor::Blue, "#it{t}_{#plus}"));
    p.add_curve(x, itp, dashed(jpacColor::Blue));
    p.add_curve(x, rtm, solid (jpacColor::Red, "#it{t}_{#minus}"));
    p.add_curve(x, itm, dashed(jpacColor::Red));
    p.set_labels("#it{s}  [GeV^{2}]", "#it{t}  [GeV^{2}]");
    p.save("tpm.pdf");
};