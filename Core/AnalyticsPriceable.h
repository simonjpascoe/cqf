#pragma once

#include "Financials.h"

using namespace Financials;

namespace Analytic {
class AnalyticsPriceable
{
public:
	AnalyticsPriceable(void) {}
	virtual ~AnalyticsPriceable(void) {}

	virtual pair<double,double> AnalyticsEvaluate(const MarketData& marketData) const = 0;
};
}

