#pragma once

#include <functional>
#include <memory>
#include <string>
#include "Financials.h"
#include "Matrix.h"
#include "FiniteDifferencePriceable.h"

using namespace Common;
using namespace Financials;

namespace Enums {

enum EvaluationMode {
	Optimize,
	Pricing
};

inline EvaluationMode EnumerateEvaluationMode(string name) {
	if (name == "Optimize") { return Optimize; }
	if (name == "Pricing") { return Pricing; }
	
	throw exception("Unsupported evaluation mode");
}

enum FiniteDifferenceScheme {
	Explicit,
	ExplicitRichardson,
	Implicit
};

inline FiniteDifferenceScheme EnumerateFiniteDifferenceScheme(string name) {
	if (name == "Explicit") { return Explicit; }
	if (name == "ExplicitRichardson") { return ExplicitRichardson; }
	if (name == "Implicit") { return Implicit; }
	throw exception("Unsupported finite difference scheme");
}

}

namespace FD {

class FiniteDifferenceParameters
{
public:
	unsigned int SpaceSteps;
	double MaxTimestepLength;
	Enums::FiniteDifferenceScheme Scheme;
	double GridWidth;
};

class FiniteDifferenceSolution
{
public:
	double Value;
	Matrix Payoff;
	Matrix Curve;
	Matrix ResidualPayoff;
	Matrix HedgeWeights;
};

class FiniteDifferenceResult
{
public:
	Matrix SpotLevels;
	pair<double, FiniteDifferenceSolution> Bid;
	pair<double, FiniteDifferenceSolution> Offer;
};


class FiniteDifferencePricer
{
public:
	FiniteDifferencePricer() { };
	virtual ~FiniteDifferencePricer() { };

	FiniteDifferenceResult Evaluate(
		const Enums::EvaluationMode evaluationMode,
		const vector<pair<FiniteDifferencePriceable*,double>>& exotic,
		const Quotes<FiniteDifferencePriceable*>& hedgeInstruments,
		const FiniteDifferenceParameters& parameters,
		const MarketData& marketData) const;

	FiniteDifferenceSolution EvaluateHedge(
		const FiniteDifferenceParameters& parameters,
		const vector<pair<FiniteDifferencePriceable*,double>>& exotic, 
		const MarketData& marketData,
		const Quotes<FiniteDifferencePriceable*>& hedgeInstruments,
		const vector<double>& hedgeWeights) const;

	const vector<double> FiniteDifferencePricer::CalculateSpaceVector(
		const double maturity,
		const MarketData& marketData,
		const FiniteDifferenceParameters& parameters) const;

private:
	pair<double, FiniteDifferenceSolution> OptimizeWeights(
		const FiniteDifferenceParameters& parameters,
		const vector<pair<FiniteDifferencePriceable*,double>>& exotic, 
		const MarketData& marketData,
		const Quotes<FiniteDifferencePriceable*>& hedgeInstruments,
		const bool solveForBid) const;

	vector<double> Kernel(
		const FiniteDifferenceParameters& parameters,
		const vector<pair<FiniteDifferencePriceable*, double>>& instruments,
		const MarketData& marketData) const;
	
	vector<double> RichardsonKernel(
		const FiniteDifferenceParameters& parameters,
		const vector<pair<FiniteDifferencePriceable*, double>>& instruments,
		const MarketData& marketData,
		const double factor) const;

	vector<double> KernelImplicit(
		const FiniteDifferenceParameters& parameters,
		const vector<pair<FiniteDifferencePriceable*, double>>& instruments,
		const MarketData& marketData) const;

	vector<pair<FiniteDifferencePriceable*, double>> invertPosition( 
		const vector<pair<FiniteDifferencePriceable*, double>>& instruments) const;

	double calculateMaturity(const vector<pair<FiniteDifferencePriceable*, double>>& instruments) const;
	double calculateMaturity(const vector<FiniteDifferencePriceable*>& instruments) const;
	double calculateMaturity(const Quotes<FiniteDifferencePriceable*>& quotes) const;
};

}