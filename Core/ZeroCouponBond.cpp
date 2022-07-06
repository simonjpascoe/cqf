#include "stdafx.h"
#include "ZeroCouponBond.h"

namespace Instruments {

ZeroCouponBond::ZeroCouponBond(double maturity) : Instrument(maturity)
{
}

double ZeroCouponBond::HeathJarrowMortonEvaluate(double dt, const Matrix& simulation, const Matrix& tenors) const
{
    // The first column of the path simulations corresponds to the short term rate
    // thus V(tau) = \exp (-dt * \sum_0^tau r_t)

    auto finalRow = (size_t)ceil(_maturity/dt + 0.5);
    if (finalRow > simulation.Rows()) {
        throw "HJM simulation is not over a sufficiently long time horizon";
    }

    double value = 0;
    for (auto i=0; i<finalRow; i++)
    {
        value += simulation(i,0);
    }
    return exp(-value * dt);
}

}