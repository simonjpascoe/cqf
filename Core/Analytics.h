#pragma once

#include "Financials.h"
#include "AnalyticsPriceable.h"

using namespace Financials;

namespace Analytic {

class Analytics
{
public:
    Analytics(void) {}
    virtual ~Analytics(void) {}

    //pair<double, double> CalculateImpliedVol(const Quote<AnalyticsPriceable*>& quote, const MarketData& marketData) const; 
};

}

