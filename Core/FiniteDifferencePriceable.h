#pragma once


namespace FD {
class FiniteDifferencePriceable
{
public:
    FiniteDifferencePriceable(void) {}
    virtual ~FiniteDifferencePriceable(void) {}

    virtual double Payoff(double spotValue) const = 0;
    virtual double Maturity() const = 0;
};
}

