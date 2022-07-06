#pragma once
#include "stdafx.h"
#include <string>
using namespace std;

namespace Enums {

    enum class InstrumentType {
        ZeroCouponBond,
        Cap,
        VanillaOption
};

inline InstrumentType EnumerateInstrument(string name) {
    if (name == "ZeroCouponBond") { return InstrumentType::ZeroCouponBond; }
    if (name == "Cap") { return InstrumentType::Cap; }
    if (name == "VanillaOption") { return InstrumentType::VanillaOption; }
    
    throw exception("Unsupported instrument type");
}

}

namespace Instruments {
class Instrument
{
public:
    Instrument(double maturity) : _maturity(maturity) { };
    virtual ~Instrument() {} ;

    double Maturity() const { return _maturity; }

protected:
    const double _maturity;
};

}