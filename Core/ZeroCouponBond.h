#pragma once
#include "Instruments.h"
#include "HeathJarrowMortonPriceable.h"

using namespace HJM;

namespace Instruments {

class ZeroCouponBond : public Instrument, public HeathJarrowMortonPriceable
{
public:
	ZeroCouponBond(double maturity);
	virtual ~ZeroCouponBond(void) {
	}

	virtual double HeathJarrowMortonEvaluate(double dt, const Matrix& simulation, const Matrix& tenors) const;

};

};

