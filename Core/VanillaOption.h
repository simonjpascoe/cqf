#pragma once
#include "Instruments.h"
#include "FiniteDifferencePriceable.h"
#include "AnalyticsPriceable.h"

using namespace std;
using namespace FD;
using namespace Analytic;

namespace Enums
{
enum OptionType
{
    Call,
    Put,
    BinaryCall,
    BinaryPut,
    Forward
};
}

namespace Instruments {

class VanillaOption : public Instrument, public FiniteDifferencePriceable, public AnalyticsPriceable
{
public:
    VanillaOption(double maturity, double strike, Enums::OptionType optionType);
    virtual ~VanillaOption(void) {}

    double Payoff(double spotValue) const;
    double Maturity() const { return Instrument::Maturity(); }

    pair<double, double> AnalyticsEvaluate(const MarketData& marketdata) const;

    double call(double spot, double volatility, double interestRate, double maturity, double strike) const;
    double put(double spot, double volatility, double interestRate, double maturity, double strike) const;
    double binaryCall(double spot, double volatility, double interestRate, double maturity, double strike) const;
    double binaryPut(double spot, double volatility, double interestRate, double maturity, double strike) const;
    double forward(double spot, double interestRate, double maturity, double strike) const;

private:    
    const double _strike;
    const Enums::OptionType _optionType;

};

}