#pragma once
#include "stdafx.h"
#include <string>
using namespace std;

namespace Enums {

enum InstrumentType {
	ZeroCouponBond,
	Cap,
	VanillaOption
};

inline InstrumentType EnumerateInstrument(string name) {
	if (name == "ZeroCouponBond") { return ZeroCouponBond; }
	if (name == "Cap") { return Cap; }
	if (name == "VanillaOption") { return VanillaOption; }
	
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