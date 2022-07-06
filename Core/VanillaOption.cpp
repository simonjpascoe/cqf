#include "stdafx.h"
#include "VanillaOption.h"

#include "Mathematics.h"

using namespace Common;
using namespace Enums;

namespace Instruments {

VanillaOption::VanillaOption(double maturity, double strike, Enums::OptionType optionType) 
    : Instrument(maturity), _strike(strike), _optionType(optionType)
{
}

double VanillaOption::Payoff(double spot) const
{
    switch (_optionType)
    {
        case OptionType::Call: return max(0.0, spot - _strike);
        case OptionType::Put: return max(0.0, _strike - spot);
        case OptionType::BinaryCall: return (spot >= _strike ? 1.0 : 0.0);
        case OptionType::BinaryPut: return (spot <= _strike ? 1.0 : 0.0);
        case OptionType::Forward: return spot - _strike;
    }

    throw exception("Unsupported option type");
}

pair<double, double> VanillaOption::AnalyticsEvaluate(const MarketData& marketData) const
{
    switch (_optionType)
    {
        case OptionType::Call: { 
            auto low = call(marketData.Spot, marketData.LowVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            auto high = call(marketData.Spot, marketData.HighVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            return make_pair(low,high);
            };
        case OptionType::Put: {
            auto low = put(marketData.Spot, marketData.LowVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            auto high = put(marketData.Spot, marketData.HighVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            return make_pair(low,high);
            };
        case OptionType::BinaryCall: {
            auto low = binaryCall(marketData.Spot, marketData.LowVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            auto high = binaryCall(marketData.Spot, marketData.HighVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            return make_pair(low,high);
        };
        case OptionType::BinaryPut: {
            auto low = binaryPut(marketData.Spot, marketData.LowVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            auto high = binaryPut(marketData.Spot, marketData.HighVolatility, marketData.InterestRate, Instrument::_maturity, _strike);
            return make_pair(low,high);
        };
        case OptionType::Forward: {
            // no volatility dependence, so it's the same either way
            auto fwd = forward(marketData.Spot, marketData.InterestRate, Instrument::_maturity, _strike);
            return make_pair(fwd, fwd);
        };
    }

    throw exception("Unsupported option type");
}


double VanillaOption::call(double spot, double volatility, double interestRate, double maturity, double strike) const
{
    auto d1 = (log(spot/strike) + (interestRate + (volatility * volatility / 2.0)) * maturity)/(volatility*sqrt(maturity));
    auto d2 = d1 - volatility * sqrt(maturity);
    return spot * StandardCumulativeNormal(d1) - StandardCumulativeNormal(d2) * strike * exp(-interestRate * maturity);
};

double VanillaOption::put(double spot, double volatility, double interestRate, double maturity, double strike) const
{
    auto d1 = (log(spot/strike) + (interestRate + (volatility * volatility / 2.0)) * maturity)/(volatility*sqrt(maturity));
    auto d2 = d1 - volatility * sqrt(maturity);
    return strike * exp(-interestRate * maturity) * StandardCumulativeNormal(-d2) - spot * StandardCumulativeNormal(-d1) ;
};

double VanillaOption::binaryCall(double spot, double volatility, double interestRate, double maturity, double strike) const
{
    auto d2 = (log(spot/strike) + (interestRate - (volatility * volatility / 2.0)) * maturity)/(volatility*sqrt(maturity));
    return exp(-interestRate * maturity) * StandardCumulativeNormal(d2);
};

double VanillaOption::binaryPut(double spot, double volatility, double interestRate, double maturity, double strike) const
{
    auto d2 = (log(spot/strike) + (interestRate - (volatility * volatility / 2.0)) * maturity)/(volatility*sqrt(maturity));
    return exp(-interestRate * maturity) * StandardCumulativeNormal(-d2);
}; 

double VanillaOption::forward(double spot, double interestRate, double maturity, double strike) const
{
    return (spot - strike) * exp(interestRate * maturity);
};

}