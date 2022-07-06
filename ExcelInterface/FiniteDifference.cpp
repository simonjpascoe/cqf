#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

#include "ExcelHelpers.h"
#include "XDictionary.h"
#include "KeyStore.h"
#include "Matrix.h"
#include "ObjectCache.h"
#include "ObjectOwner.h"
#include "Instruments.h"
#include "Portfolio.h"
#include "Financials.h"

#include "FiniteDifferencePricer.h"
#include "FiniteDifferencePriceable.h"


using namespace xlw;
using namespace std;
using namespace Common;
using namespace Financials;
using namespace FD;
using namespace Instruments;

extern "C" {
    LPXLFOPER EXCEL_EXPORT xlFDMEvaluate(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto evaluationMode = Enums::EnumerateEvaluationMode(dict.Lookup("evaluation_mode").StringValue());

        auto xParameters = dict.Lookup("fdm_parameters").StringValue();
        auto xMarketData = dict.Lookup("market_data").StringValue();

        vector<pair<FiniteDifferencePriceable*,double>> exotic;
        auto keyStore = dict.Lookup("portfolio").StringValue();
        auto portfolio = ObjectCache<Portfolio>::Instance().Lookup(keyStore);
        for (auto k : portfolio.Values()) {
            if (k.second == 0) continue; // don't add zero weighted
            auto instrument = ObjectCache<shared_ptr<Instruments::Instrument>>::Instance().Lookup(k.first);
            auto two = dynamic_cast<FiniteDifferencePriceable*>(instrument.get());
            exotic.push_back(make_pair(two, k.second));
        }

        auto cParameters = ObjectCache<FiniteDifferenceParameters>::Instance().Lookup(xParameters);
        auto cMarketData = ObjectCache<MarketData>::Instance().Lookup(xMarketData);

        Quotes<FiniteDifferencePriceable*> cHedgeInstruments;

        if (evaluationMode == Enums::Optimize) {
            auto xQuotes = dict.Lookup("market_quotes").StringValue();
            auto quotes = ObjectCache<Quotes<string>>::Instance().Lookup(xQuotes);
            for (auto q : quotes.All())
            {
                auto instrument = ObjectCache<shared_ptr<Instruments::Instrument>>::Instance().Lookup(q.Instrument);
                auto two = dynamic_cast<FiniteDifferencePriceable*>(instrument.get());
                Quote<FiniteDifferencePriceable*> quote;
                quote.Bid = q.Bid;
                quote.Offer = q.Offer;
                quote.Instrument = two;
                cHedgeInstruments.Add(quote);
            }
        }

        FiniteDifferencePricer fde;
        auto result = fde.Evaluate(evaluationMode, exotic, cHedgeInstruments, cParameters, cMarketData);

        // bid and offer results
        auto createNested = [&evaluationMode] (const pair<double, FiniteDifferenceSolution>& f) {
            XDictionary output;
            vector<pair<string, freeF>> children;
            output.Add("payoff", cache_child(f.second.Payoff,children));
            output.Add("curve", cache_child(f.second.Curve,children));

            if (evaluationMode == Enums::Optimize) {
                output.Add("hedge_weights", cache_child(f.second.HedgeWeights, children));
                output.Add("residual_payoff", cache_child(f.second.ResidualPayoff,children));
                output.Add("fit_quality", CellMatrix(f.first));
            }

            output.Add("value", CellMatrix(f.second.Value));
            return output;
        };

        auto bid_xd = createNested(result.Bid);
        auto offer_xd = createNested(result.Offer);

        // top level
        XDictionary output;
        vector<pair<string, freeF>> children;
        output.Add("spot_levels",cache_child(result.SpotLevels, children));
        output.Add("bid", cache_child(bid_xd,children));
        output.Add("offer", cache_child(offer_xd,children));

        // ownership management
        auto key = cache_dictionary(output, children);

        return XlfOper(key);

        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlFDMEvaluateHedge(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xParameters = dict.Lookup("fdm_parameters").StringValue();
        auto xMarketData = dict.Lookup("market_data").StringValue();
        auto xHedgeWeights = dict.Lookup("hedge_weights").StringValue();

        vector<pair<FiniteDifferencePriceable*,double>> exotic;
        auto maturity = 0.0;
        auto keyStore = dict.Lookup("portfolio").StringValue();
        auto portfolio = ObjectCache<Portfolio>::Instance().Lookup(keyStore);
        for (auto k : portfolio.Values()) {
            if (k.second == 0) continue; // don't add zero weighted
            auto instrument = ObjectCache<shared_ptr<Instruments::Instrument>>::Instance().Lookup(k.first);
            maturity = max(maturity, instrument->Maturity());
            auto two = dynamic_cast<FiniteDifferencePriceable*>(instrument.get());
            exotic.push_back(make_pair(two, k.second));
        }

        auto cParameters = ObjectCache<FiniteDifferenceParameters>::Instance().Lookup(xParameters);
        auto cMarketData = ObjectCache<MarketData>::Instance().Lookup(xMarketData);
        auto cHedgeWeights = ObjectCache<Matrix>::Instance().Lookup(xHedgeWeights);

        auto xQuotes = dict.Lookup("market_quotes").StringValue();
        auto quotes = ObjectCache<Quotes<string>>::Instance().Lookup(xQuotes);
        Quotes<FiniteDifferencePriceable*> cHedgeInstruments;
        for (auto q : quotes.All())
        {
            auto instrument = ObjectCache<shared_ptr<Instruments::Instrument>>::Instance().Lookup(q.Instrument);
            auto two = dynamic_cast<FiniteDifferencePriceable*>(instrument.get());
            maturity = max(maturity, two->Maturity());
            Quote<FiniteDifferencePriceable*> quote;
            quote.Bid = q.Bid;
            quote.Offer = q.Offer;
            quote.Instrument = two;
            cHedgeInstruments.Add(quote);
        }

        FiniteDifferencePricer fde;
        auto result = fde.EvaluateHedge(cParameters, exotic, cMarketData, cHedgeInstruments, cHedgeWeights.ToVector(0));
        auto space = Matrix::FromVector(fde.CalculateSpaceVector(maturity, cMarketData, cParameters));
        XDictionary output;
        vector<pair<string, freeF>> children;
        
        output.Add("spot_levels",cache_child(space, children));
        output.Add("payoff", cache_child(result.Payoff,children));
        output.Add("curve", cache_child(result.Curve,children));
        output.Add("residual_payoff", cache_child(result.ResidualPayoff,children));
        output.Add("value", CellMatrix(result.Value));

        // ownership management
        auto key = cache_dictionary(output, children);

        return XlfOper(key);

        EXCEL_END;
    }
    
    LPXLFOPER EXCEL_EXPORT xlCreateParameters(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xSpaceSteps = dict.Lookup("space_steps").NumericValue();
        auto xGridWidth = dict.Lookup("grid_width").NumericValue();
        auto xScheme = dict.Lookup("scheme").StringValue();
        auto xMaxTimestep = dict.Lookup("max_timestep_length").NumericValue();
    
        FiniteDifferenceParameters fdp;
        fdp.SpaceSteps = (unsigned int)xSpaceSteps;
        fdp.GridWidth = xGridWidth;
        fdp.Scheme = Enums::EnumerateFiniteDifferenceScheme(xScheme);
        fdp.MaxTimestepLength = xMaxTimestep;
    
        auto key = cache_item(fdp);
        return XlfOper(key);

        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlQuotesCreate(LPXLFOPER source)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;
        
        XlfOper xSource(source);

        auto cells = xSource.AsCellMatrix();
        if (cells.ColumnsInStructure() !=3)
            { THROW_XLW("Columns must be three"); }

        Quotes<string> quotes;
        for (size_t i = 0; i<cells.RowsInStructure(); i++)
        {
            Quote<string> quote;
            quote.Instrument = cells(i,0).StringValue();
            quote.Bid = cells(i,1).NumericValue();
            quote.Offer = cells(i,2).NumericValue();
            quotes.Add(quote);
        }

        auto key = cache_item(quotes, "Quotes");
        
        return XlfOper(key);

        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlQuoteCreate(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;
        
        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xInstrument = dict.Lookup("instrument").StringValue();
        auto xBid = dict.Lookup("bid").NumericValue();
        auto xOffer = dict.Lookup("offer").NumericValue();

        Quote<string> quote;
        quote.Instrument = xInstrument;
        quote.Bid = xBid;
        quote.Offer = xOffer;

        auto key = cache_item(quote, "Quote");
        
        return XlfOper(key);

        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlCreateMarketDataContainer(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xHighVol = dict.Lookup("high_volatility").NumericValue();
        auto xLowVol = dict.Lookup("low_volatility").NumericValue();
        auto xInterestRate = dict.Lookup("interest_rate").NumericValue();
        auto xSpot = dict.Lookup("spot").NumericValue();
        
        if (xLowVol > xHighVol) { THROW_XLW("Low Volatility must be lower than High Volatility");}

        MarketData md;
        md.HighVolatility = xHighVol;
        md.LowVolatility = xLowVol;
        md.InterestRate = xInterestRate;
        md.Spot = xSpot;

        auto key = cache_item(md);
        return XlfOper(key);

        EXCEL_END;
    }
}

namespace {
    XLRegistration::Arg argsXlXDictionary[] = {
        { "dictionary", "Dictionary", "XLF_OPER" }
    };

    XLRegistration::XLFunctionRegistrationHelper register__fdm_evalulate(
        "xlFDMEvaluate", "fdm_evaluate", "Evaluate portfolio with fdm",
        "CQF", argsXlXDictionary, 1);

    XLRegistration::XLFunctionRegistrationHelper register__fdm_evalulate_hedge(
        "xlFDMEvaluateHedge", "fdm_evaluate_hedge", "Evaluate a hedged portfolio with fdm",
        "CQF", argsXlXDictionary, 1);

    XLRegistration::XLFunctionRegistrationHelper register__create_parameters(
        "xlCreateParameters", "fdm_parameters_create", "Create engine parameters",
        "CQF", argsXlXDictionary, 1);
    
    XLRegistration::XLFunctionRegistrationHelper register__create_marketdata(
        "xlCreateMarketDataContainer", "fdm_marketdata_create", "Create uncertain vol fdm market data",
        "CQF", argsXlXDictionary, 1);

    XLRegistration::XLFunctionRegistrationHelper register__create_quotes(
        "xlQuotesCreate", "quotes_create", "Create quotes",
        "CQF", argsXlXDictionary, 1);

    XLRegistration::XLFunctionRegistrationHelper register__create_quote(
        "xlQuoteCreate", "quote_create", "Create quote",
        "CQF", argsXlXDictionary, 1);
}