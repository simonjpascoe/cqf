#pragma once
#include "Instruments.h"
#include "HeathJarrowMortonPriceable.h"

using namespace HJM;

namespace Instruments {

class Cap : public Instrument, public HeathJarrowMortonPriceable
{
public:
	Cap(double maturity, double strike);
	virtual ~Cap(void) {
	}

	virtual double HeathJarrowMortonEvaluate(double dt, const Matrix& simulation, const Matrix& tenors) const;

private:
	double calculateCaplet(size_t ti, const Matrix& simulation, const Matrix& tenors) const;
	const double _strike;
};

}