#pragma once
#include "Matrix.h"

using namespace Common;

namespace HJM {
class HeathJarrowMortonPriceable
{
public:
	HeathJarrowMortonPriceable(void) {}
	virtual ~HeathJarrowMortonPriceable(void) {}

	virtual double HeathJarrowMortonEvaluate(double dt, const Matrix& simulation, const Matrix& tenors) const = 0;
};
}

