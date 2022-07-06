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

#include "ZeroCouponBond.h"
#include "Cap.h"
#include "VanillaOption.h"

using namespace xlw;
using namespace std;
using namespace Enums;
using namespace Instruments;


extern "C" {
    LPXLFOPER EXCEL_EXPORT xlInstrumentCreate(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xProduct = Enums::EnumerateInstrument(dict.Lookup("class").StringValue());

        string key = "no product created";
        auto xMaturity = dict.Lookup("maturity").NumericValue();

        switch (xProduct) {
        case InstrumentType::ZeroCouponBond:
            {
                shared_ptr<Instrument> zcb(new ZeroCouponBond(xMaturity));
                key = cache_item(zcb,"Instrument");
                break;
            }
        case InstrumentType::Cap:
            {
                auto xStrike = dict.Lookup("strike").NumericValue();
                shared_ptr<Instrument> zcb(new Cap(xMaturity,xStrike));
                key = cache_item(zcb,"Instrument");
                break;
            }    
        case InstrumentType::VanillaOption:
            {
                auto xStrike = dict.Lookup("strike").NumericValue();
                auto xType = dict.Lookup("type").StringValue();
                OptionType ty;
                if (xType == "Call") {ty = OptionType::Call; }
                else if (xType == "Put") {ty = OptionType::Put;}
                else if (xType == "BinaryCall") { ty = OptionType::BinaryCall; }
                else if (xType == "BinaryPut") { ty = OptionType::BinaryPut;}
                else if (xType == "Forward") { ty = OptionType::Forward;}
                else THROW_XLW("Unrecognised option type");
                shared_ptr<Instrument> option(new VanillaOption(xMaturity, xStrike, ty));
                key = cache_item(option,"Instrument");
                break;
            }
        }
        
        return XlfOper(key);
        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlPortfolioCreate(LPXLFOPER source)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xSource(source);

        auto cells = xSource.AsCellMatrix();

        Portfolio portfolio;
        for (size_t j = 0; j<cells.RowsInStructure(); j++)
        {
            portfolio.Add(cells(j,0).StringValue(), cells(j,1).NumericValue());
        }

        auto key = cache_item(portfolio);
        
        return XlfOper(key);
        EXCEL_END;
    }
}
namespace {
    XLRegistration::Arg argsXlXDictionary[] = {
        { "dictionary", "Dictionary", "XLF_OPER" }
    };

    XLRegistration::XLFunctionRegistrationHelper register__instrument_create(
        "xlInstrumentCreate", "instrument_create", "Create a priceable instrument",
        "CQF", argsXlXDictionary, 1);

    XLRegistration::Arg argsXlPortfolio[] = {
        { "instruments", "Mapping", "XLF_OPER" }
    };

    XLRegistration::XLFunctionRegistrationHelper register__portfolio_create(
        "xlPortfolioCreate", "portfolio_create", "Create a portfolio",
        "CQF", argsXlPortfolio, 1);
}