#include "stdafx.h"
#include "Cap.h"
#include "ZeroCouponBond.h"

namespace Instruments {

Cap::Cap(double maturity, double strike) : Instrument(maturity), _strike(strike)
{
}

double Cap::HeathJarrowMortonEvaluate(double dt, const Matrix& simulation, const Matrix& tenors) const
{
	// Calculate and sum up all the discounted caplet prices
	auto finalRow = (size_t)ceil(_maturity/dt + 0.5);
	if (finalRow > simulation.Rows()) {
		throw "HJM simulation is not over a sufficiently long time horizon";
	}

	double value = 0;
	for (auto i=0; i<finalRow; i++)
	{
		auto caplet = calculateCaplet(i, simulation, tenors) * dt;
		ZeroCouponBond zcb(i*dt);
		value += caplet * zcb.HeathJarrowMortonEvaluate(dt, simulation, tenors);
	}

	return value;
}

double Cap::calculateCaplet(size_t ti, const Matrix& simulation, const Matrix& tenors) const
{
	// 1. Calculate forward libor at time represented by row ti using simple integration
	auto libor = 0.0;
	for (auto i = 0; i< tenors.Cols()-1; i++)
	{
		libor += (tenors(0,i+1) - tenors(0,i)) * (simulation(ti,i+1) + simulation(ti,i));
	}
	libor = libor * 0.5/tenors(0,tenors.Cols()-1);

	// 2. calculate the caplet value
	return max(libor - _strike, 0);
}

}