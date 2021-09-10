#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

#include "ExcelHelpers.h"
#include "XDictionary.h"
#include "KeyStore.h"

#include "Analytics.h"
#include "AnalyticsPriceable.h"
#include "Financials.h"
#include "Instruments.h"
#include "VanillaOption.h"

using namespace xlw;
using namespace std;

using namespace Analytic;
using namespace Financials;
using namespace Common;
using namespace Instruments;

extern "C" {
	LPXLFOPER EXCEL_EXPORT xlAnalyticsEvaluate(LPXLFOPER dictionary)
	{
		EXCEL_BEGIN;
		NO_FUNC_WIZARD;

		XlfOper xDictionary(dictionary);
		auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

		auto xInstrument = dict.Lookup("instrument").StringValue();
		auto xMarketData = dict.Lookup("market_data").StringValue();

		auto instrument = ObjectCache<shared_ptr<Instruments::Instrument>>::Instance().Lookup(xInstrument);
		auto ae = dynamic_cast<AnalyticsPriceable*>(instrument.get());

		auto marketData = ObjectCache<MarketData>::Instance().Lookup(xMarketData);

		auto result = ae->AnalyticsEvaluate(marketData);

		XDictionary output;
		output.Add("low", CellMatrix(result.first));
		output.Add("high", CellMatrix(result.second));

		// ownership management
		auto key = cache_item(output);
		return XlfOper(key);

		EXCEL_END;
	}

	LPXLFOPER EXCEL_EXPORT xlBlackScholesCalculate(
		LPXLFOPER iType, 
		LPXLFOPER iSpot,
		LPXLFOPER iLowVolatility,
		LPXLFOPER iHighVolatility,
		LPXLFOPER iInterestRate,
		LPXLFOPER iMaturity,
		LPXLFOPER iStrike)
	{
		EXCEL_BEGIN;

		auto type  = string(XlfOper(iType).AsString());
		Enums::OptionType ty;
		if (type == "Call") {ty = Enums::Call; }
		else if (type == "Put") {ty = Enums::Put;}
		else if (type == "BinaryCall") { ty = Enums::BinaryCall; }
		else if (type == "BinaryPut") { ty = Enums::BinaryPut;}
		else if (type == "Forward") { ty = Enums::Forward;}
		else THROW_XLW("Unrecognised option type");

		VanillaOption option(XlfOper(iMaturity).AsDouble(), XlfOper(iStrike).AsDouble(), ty);
		
		MarketData m;
		m.LowVolatility = XlfOper(iLowVolatility).AsDouble();
		m.HighVolatility = XlfOper(iHighVolatility).AsDouble();
		m.Spot = XlfOper(iSpot).AsDouble();
		m.InterestRate = XlfOper(iInterestRate).AsDouble();

		auto result = option.AnalyticsEvaluate(m);

		XDictionary output;
		output.Add("low", CellMatrix(result.first));
		output.Add("high", CellMatrix(result.second));

		// ownership management
		auto key = cache_item(output);
		return XlfOper(key);

		EXCEL_END;
	}
}

namespace {
	XLRegistration::Arg argsXlDictionary[] = {
        { "dictionary", "Input dictionary", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__analytics_evaluate(
        "xlAnalyticsEvaluate", "analytics_evaluate", "Price using analytics",
        "CQF", argsXlDictionary, 1);

	XLRegistration::Arg argsXlBlackScholesCalculate[] = {
        { "type", "Type", "XLF_OPER" },
		{ "spot", "Spot", "XLF_OPER" },
		{ "low_volatility", "Low Volatility", "XLF_OPER" },
		{ "high_volatility", "High Volatility", "XLF_OPER" },
		{ "interest_rate", "Interest Rate", "XLF_OPER" },
		{ "maturity", "Maturity", "XLF_OPER" },
		{ "strike", "Strike", "XLF_OPER" }
    };

	XLRegistration::XLFunctionRegistrationHelper register__blackScholes_evaluate(
        "xlBlackScholesCalculate", "blackScholes_evaluate", "Price using BS formulae",
        "CQF", argsXlBlackScholesCalculate, 7);
}